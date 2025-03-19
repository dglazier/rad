#pragma once

/*!
 *  Action class for filling a vector of histograms
 *  Only 1 is filled per event, indexed by "bin" in Exec
 */
#include <TROOT.h>
#include <TH1.h>
#include <ROOT/TSeq.hxx>
#include <THn.h>
#include <ROOT/RDF/RActionImpl.hxx>
#include <ROOT/RDF/RMergeableValue.hxx>
#include <memory>
#include <stdexcept>
#include <vector>
   
namespace rad{
  
  namespace histo{


    /*
     * Template on the data type for intermediate THnT<T>
     * this is only used for processing and final histograms
     * are given as TH1s
     */
    template <typename T >
    class SplitHistoHelper : public ROOT::Detail::RDF::RActionImpl<SplitHistoHelper<T>> {
      
    public:
      /// This is a handy, expressive shortcut.
      using THn_t = THnT<T>;

      /// Dfeine ptrs and vector used
      using Histogram_t = TH1; // Abstract histogram class, so can use TH1D, TH2D,..
      using Histogram_ptr = std::shared_ptr<Histogram_t>; 
      /// This type is a requirement for every helper.
      using Result_t = std::vector<Histogram_ptr>; //result is vector of different histograms

      using Result_ptr = std::shared_ptr<Result_t>;
      using THn_ptr = std::shared_ptr<THn_t>;
      
      
    public:
      /// This constructor takes all the parameters necessary to build the THs.
      /// hists should be a vector of TH1s for the same variable
      template <typename TH >
      SplitHistoHelper(const std::vector<TH>& hists )
      {
	  // std::cout<< "SplitHistoHelper construct" <<std::endl;
	  _Nhistos =hists.size();
	  _resultHists=std::make_shared<Result_t>();
	  
	  //create results ptrs
	  for(auto& hist:hists){
	    _resultHists->emplace_back( std::make_shared<TH>( hist) );
	  }
	  
	  const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
	  _histos.resize(nSlots);
	  for (auto i : ROOT::TSeqU(nSlots)) {
	    //copy histograms into vector shared_ptr, so have nSlots copies of hists
	    std::vector<THn_ptr> slotHists(_Nhistos);
	    for(size_t j=0; j<_Nhistos; ++j){//for each histogram convert from TH1 to THn_t
	      auto& his = hists[j];
	      auto histPtr = std::make_shared<THn_t>(*dynamic_cast<THn_t*>(THn::CreateHn(his.GetName(),his.GetTitle(),&his)));
	      slotHists[j]=( histPtr );// vector of splits
	    }
	    _histos[i]=(slotHists); //vector of slots
	    (void)i;
	  }
	  // std::cout<< "SplitHistoHelper done " <<_histos.size()<<" "<<_resultHists->size()<<" "<<_Nhistos<<std::endl;
	}
      SplitHistoHelper(SplitHistoHelper &&) = default;
      SplitHistoHelper(const SplitHistoHelper &) = delete;
      
      Result_ptr GetResultPtr() const { return  _resultHists; }
      void Initialize() {}
      void InitTask(TTreeReader *, unsigned int) {}

      /// This is a method executed at every entry
      template <typename... ColumnTypes>
      void Exec(unsigned int slot,int bin,  ColumnTypes... values)
	{
	  if(bin<0||bin>=static_cast<int>(_Nhistos)) return; //event out of range
	  // std::cout<< "SplitHistoHelper Exec " <<_Nhistos<<" "<<slot<<" "<<_histos[slot].size()<<" "<<bin<<std::endl;
	  
	  //Since TH* fill expects doubles, we build it passing through elements of  std::array<double>.
	  //array should have 1 entry for 1D hist
	  std::array<double, sizeof...(ColumnTypes)> valuesArr{static_cast<double>(values)...};

	  //Only want to fill 1 histogram corresponding to element bin in _histos
	  _histos[slot].at(bin)->Fill(valuesArr.data()); 
	}
      /// This method is called at the end of the event loop. It is used to merge all the internal THnTs which
      /// were used in each of the data processing slots.
      void Finalize()
      {

	//loop over slots
	for (auto slot : ROOT::TSeqU(0, _histos.size())) {
	  //loop over histograms and sum each slot
	  for (size_t ihist=0; ihist<_Nhistos; ++ihist){
	    TH1* temp = nullptr;
	    if(_histos[slot].at(ihist)->GetNdimensions()==1){
	      temp = _histos[slot].at(ihist)->Projection(0,"OE");
	    }
	    else if(_histos[slot].at(ihist)->GetNdimensions()==2){
	      temp = _histos[slot].at(ihist)->Projection(1,0,"OE");
	    }
	    else if(_histos[slot].at(ihist)->GetNdimensions()==3){
	      temp = _histos[slot].at(ihist)->Projection(2,1,0,"OE");
	    }

	    _resultHists->at(ihist)->Add( (temp) );
	    delete temp;temp=nullptr;
	  }
	}
      }
 
      std::string GetActionName(){
	return "SplitHistoHelper";
      }

         private:
      
      std::vector< std::vector<THn_ptr> > _histos; // one per data processing slot
      size_t _Nhistos=0;
      Result_ptr _resultHists;
 
    };
 
  }
}
