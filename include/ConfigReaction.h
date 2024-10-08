#pragma once

//!  Base configuration class for electron scattering reactions

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
*/
#include "DefineNames.h"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

namespace rad{
  namespace config{

    //! Code simplifications
  
    using ROOT::RVecI;
  
    /*! RDFstep is used to update the dataframe handle after each filter or define step */
    // using RDFstep = ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> ;
    using RDFstep = ROOT::RDF::RNode ;

    /*! RVecIndexMap is gives access to specific particle indices in action functions.
      Use reaction part index getters e.g. ScatEleIdx() to get the index from the map.
      e.g. rdf::RVecIndexMap& react; ... ; auto index = react[ScatEleIdx()]; 
    */
    using RVecIndexMap = ROOT::RVec<ROOT::RVecI>;


    void PrintDefinedColumnNames(RDFstep  df){
      std::cout<<"Print Column Names : ";
      auto cols =  df.GetDefinedColumnNames();
      for(auto& col:cols){
	std::cout<<col<<", ";
      }
      cout<<"\n";
    }

    bool ColumnExists(const string& col,RDFstep  df){
      auto cols =  df.GetDefinedColumnNames();
      if(std::find(cols.begin(),cols.end(),col)==cols.end()) return false;
      return true;
    }
    
    std::string as_string(std::string_view v) { 
      return {v.data(), v.size()}; 
    }
    
    //! Class definition

    class ConfigReaction {


    public:

      //     ConfigReaction(const std::string& treeName, const std::string& fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ) : _orig_df{treeName,{fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data(),fileNameGlob.data()},columns},_curr_df{_orig_df}{
      ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ) : _orig_df{treeName,{fileNameGlob.data()},columns},_curr_df{_orig_df},_base_df{_orig_df}{
      }
      ConfigReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ) : _orig_df{treeName,filenames,columns},_curr_df{_orig_df},_base_df{_orig_df}{
      }
      

      /** 
       * Get the current dataframe to add further actions
       */
      RDFstep CurrFrame(){return _curr_df;}
      /** 
       * Set the current dataframe after adding further actions
       */
      void setCurrFrame(RDFstep  df){ _curr_df = df ;}

      /**
       * Interface to RDataFrame Define
       */
      void Define(const string_view name,const string& expression){
	setCurrFrame(CurrFrame().Define(name,expression));
      }
      //  void Define(const TString& name,const string& expression){
      // 	setCurrFrame(CurrFrame().Define(name,expression));
      // }
      template<typename Lambda>
      void Define(const string_view name,Lambda&& func,const ROOT::RDF::ColumnNames_t&  columns ){
	setCurrFrame(CurrFrame().Define(name,func,columns));
      }

      void DefineForAllTypes(const string& name,const string& expression){
	for(auto &atype:_type_comps){
	  TString type_expr = expression.data();
	  type_expr.ReplaceAll("components_p4",atype.second["components_p4"]);
	  type_expr.ReplaceAll("components_p3",atype.second["components_p3"]);
	  Define(atype.first + name.data(),type_expr.Data());
	}
      }
      template<typename Lambda>
      void DefineForAllTypes(const string& name,Lambda&& func,const ROOT::RDF::ColumnNames_t&  columns ){
	for(auto &atype:_type_comps){
	  ROOT::RDF::ColumnNames_t type_cols;
	  for(auto& acol:columns){
	    type_cols.push_back(atype.first + acol);
	  }
	  Define(atype.first + name.data(),func,type_cols);
	}
      }
 
      /**
       * Interface to RDataFrame Redefine
       */
      void RedefineExpr(const string& name,const string& expression){
	setCurrFrame(CurrFrame().Redefine(name,expression));
      }
      void Redefine(const string& name,const string& expression){
	setCurrFrame(CurrFrame().Redefine(name,expression));
      }
      template<typename Lambda>
      void Redefine(const string& name,Lambda&& func,const ROOT::RDF::ColumnNames_t& columns = {}){
	setCurrFrame(CurrFrame().Redefine(name,func,columns));
      }
      /**
       * Interface to RDataFrame Redefine via any aliases that may be used
       */
      void RedefineViaAlias(const string& alias,const string& expression){
	Redefine(_aliasMap[alias],expression);
      }
      template<typename Lambda>
      void RedefineViaAlias(const string& alias,Lambda&& func,const ROOT::RDF::ColumnNames_t& columns ){
	Redefine(_aliasMap[alias],func,columns);
      }
      /** 
       * Add an alias for a branch and update the current frame to the aliased one
       */
      void setBranchAlias(const string& old_name,const string& new_name){
	_aliasMap[new_name] = old_name;
	setCurrFrame(CurrFrame().Alias(new_name,old_name));
      }
      /**
       * Interface to RDataFrame Filter
       */
      template<typename Lambda>
      void Filter(Lambda&& func, const ROOT::RDF::ColumnNames_t& columns = {},std::string name = "" ){
	setCurrFrame(CurrFrame().Filter(func,columns,name));
      }
      void Filter(std::string& expression,std::string 	name = "" ){
	setCurrFrame(CurrFrame().Filter(expression,name));
      }
     
      /**
       * Make a snapshot of newly defined columns
       */
      void Snapshot(const string& filename){
	auto cols = CurrFrame().GetDefinedColumnNames();
	RemoveSnapshotColumns(cols);
	CurrFrame().Snapshot("rad_tree",filename, cols );
      }
      virtual void RemoveSnapshotColumns(std::vector<string>& cols){
	cols.erase(std::remove(cols.begin(), cols.end(), names::ReactionMap() ), cols.end());
      }
      /** 
       * Set constant index in collection for particle
       * This assumes constant position in collection (e.g in some HepMC3 files)
       * and update the current frame to the aliased one
       */
      void setParticleIndex(const string& particle, const int idx, int pdg=0 ){
	setCurrFrame(CurrFrame().Define(particle,[idx](){return idx;}));
	if(pdg!=0) Define(particle+"_OK",Form("rec_pid[%s]==%d",particle.data(),pdg));
      }
      
      /** 
       * Set function, func,  which defines variable index in collection for particle
       * The given function func is responsible for giving the correct position of particle in the collection
       * Names of required branches or identifiers given as vector of column names : columns
       * and update the current frame to the aliased one
       */
      template<typename Lambda>
      void setParticleIndex(const string& particle, Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {}, int pdg=0 ){
	setCurrFrame(CurrFrame().Define(particle,func,columns));
	if(pdg!=0) Define(string(particle)+"_OK",Form("rec_pid[%s]==%d",particle.data(),pdg));
      }

      /**
       * Make map that links particle names to indices in user functions
       * in C++ functions you can use the RVecIndexMap object indexed by 
       * name of the reaction component you need
       *
       * This function must be implmented by a derived class
       */
      virtual void makeParticleMap() = 0;
     // /**
     //   * Make map that links beam names to indices in user functions
     //   * in C++ functions you can use the RVecIndexMap object indexed by 
     //   * name of the reaction component you need
     //   *
     //   * This function must be implmented by a derived class
     //   */
     //  virtual void makeBeamMap() = 0;


       /**
       * Collect constant indices for final state mesons and baryons
       */
      void setMesonIndices(const RVecI& indices){
	Define(as_string(names::Mesons()),[indices](){return indices;},{});
      } 
      void setBaryonIndices(const RVecI& indices){
	Define(as_string(names::Baryons()),[indices](){return indices;},{});
      }

       /**
       * Collect variable indices for final state mesons and baryons
       */
      void setMesonParticles(const  ROOT::RDF::ColumnNames_t& particles){
	setGroupParticles(as_string(names::Mesons()),particles);
      }
      void setBaryonParticles(const  ROOT::RDF::ColumnNames_t& particles){
	setGroupParticles(as_string(names::Baryons()),particles);
      }
     
      /**
       * Generic function to collect variable indices for groups of particle
       * name = identifier of the indice collection, 
       * particles = list of particle names known to ConfigureReaction
       */

      /*
      //Would like the code below to work. Compiles OK, but the dataframe infers 0 arguments, rather than particles.size()
      template <typename... Args>
      void setGroupParticles(const string& name,const  ROOT::RDF::ColumnNames_t& particles){
	auto ncols=particles.size();
	//create lambda to return current indices of particles
	//args will be the data in the columns defined by particles
	auto func = [ncols](Args&&... args) {
	  RVecI vec;
	  vec.reserve(ncols);
	  (vec.emplace_back(std::forward<Args>(args)), ...);
	  return vec;};

	//need to test which of these works fastest
	// auto func = [ncols](Args&&... args) {
	//   RVecI vec;
	//   vec.reserve(ncols);
	//   for (auto i = 0UL; i < size; ++i) {
	//     ret.emplace_back(args[i]...);
	//   }
	  
	// auto func = [ncols](Args&&... args) {
	//   RVecI vec(ncols);
	//   for (auto i = 0UL; i < size; ++i) {
	//     ret[i]=(args[i]...);
	//   }
	//  return vec;};
	
   	setCurrFrame(CurrFrame().Define(names::Baryons(),func,particles));
   
      }
      */
      
      void setGroupParticles(const string& name,const ROOT::RDF::ColumnNames_t &particles){
	//should be able to do this with variadic args, but don't know how to pass as argument to lambda
	switch( particles.size() ){
	case 1 : 
	  Define(name,[](const int p0){return RVecI{p0};},particles);
	  break;
	case 2 : 
	  Define(name,[](const int p0,const int p1){return RVecI{p0,p1};},particles);
	  break;
	case 3 : 
	  Define(name,[](const int p0,const int p1,const int p2){return RVecI{p0,p1,p2};},particles);
	  break;
	case 4 : 
	  Define(name,[](const int p0,const int p1,const int p2,const int p3){return RVecI{p0,p1,p2,p3};},particles);
	  break;
	case 5 : 
	  Define(name,[](const int p0,const int p1,const int p2,const int p3,const int p4){return RVecI{p0,01,p2,p3,p4};},particles);
	  break;
	case 6 : 
	  Define(name,[](const int p0,const int p1,const int p2,const int p3,const int p4,const int p5){return RVecI{p0,01,p2,p3,p4,p5};},particles);
	  break;
	default:
	  std::cerr<<"setGroupParticles only defined up to 6 particles! "<<std::endl;
	  exit(0);
	  break;
	}
	     
      }
      /**
       * Set all Pid (aka PDG) values to -1 so particles ignored
       */
      template<typename Lambda>
      void  ApplyPidMask(const string& mask_name,const string& type,Lambda&& func,const ROOT::RDF::ColumnNames_t  columns = {}){
	//Define a column for the mask variable
	Define(mask_name,func,columns);
	//Apply mask to PID column, so entries indexed by -1 will not be used
	RedefineViaAlias(type+"_"+"pid",[](const RVecI& pid,const RVecI& mask){
	  return ROOT::VecOps::Where(mask,pid,-1);},
	  {type+"_"+"pid",mask_name.data()});
      }
      /**
       * create shortcut string for 3 and 4 momentum components
       */
      void AddType(const string& atype){
	if(_primary_type.empty()==true) _primary_type=atype;
	_type_comps[atype]["components_p4"] = Form("%spx,%spy,%spz,%sm",atype.data(),atype.data(),atype.data(),atype.data());
	_type_comps[atype]["components_p3"] = Form("%spx,%spy,%spz",atype.data(),atype.data(),atype.data());
     }
      /**
       *  change name of this first type column to same name without type prefix
       */
      void AliasToPrimaryType(const string& name){
	if(_primary_type.empty()==true) return;
	std::cout<<"AliasToPrimaryType "<<_primary_type<<" "<<name<<std::endl;
	setBranchAlias(_primary_type+name,name);
      }
     /**
      *  copy name of this first type column to same name without type prefix
       */
      void CopyToPrimaryType(const string& name){
	if(_primary_type.empty()==true) return;
	std::cout<<"CopyToPrimaryType "<<_primary_type<<" "<<name<<std::endl;
	Define(name,Form("return %s;",(_primary_type+name).data() ) );
      }
      enum class ColType{Undef,Int,UInt,Float,Double,Short,Bool,Long};
      
      ColType DeduceColumnVectorType(const string& name){
	TString col_type = CurrFrame().GetColumnType(name);
	if(col_type.Contains("Int_t")) return ColType::Int;
	if(col_type.Contains("UInt_t")) return ColType::UInt;
 	if(col_type.Contains("Float_t")) return ColType::Float;
 	if(col_type.Contains("Double_t")) return ColType::Double;
 	if(col_type.Contains("Short_t")) return ColType::Short;
 	if(col_type.Contains("Bool_t")) return ColType::Bool;
 	if(col_type.Contains("Long_t")) return ColType::Long;
	return ColType::Undef;
      }

      std::map<string, std::map<string,string>> GetTypes() const {return _type_comps;}

      /**
       * get or set the base df
       * this should include the configuration of the input
       * data to rad format (arrays of px,px,pz,m)
       * sharing this allows parallisation between reactions
       */
      void setBaseFrame(RDFstep step){
	_base_df = step;
      }
      void setMyBaseFrame(){
	_base_df = CurrFrame();
      }
      RDFstep getBaseFrame() const{
	return _base_df;
      }
    protected:

      const std::map<string,string>& AliasMap() const {return _aliasMap;}
      
    private :
    
      /**
       * _bare_df
       * Base dataframe constructed prior to adding actions
       * It should not be necessary to use this object for anything
       */
      ROOT::RDataFrame _orig_df;
      /** 
       * _curr_df 
       * Handle for the current dataframe state. 
       * This includes all added defines and filters so far.
       * Additional actions will be applied to this.
       */
      RDFstep _curr_df;
      /** 
       * _base_df 
       * Handle for the basic dataframe state. 
       * This will have columns organised for processing.
       * Typically after a SetAlias call.
       */
      RDFstep _base_df;

      /**
       * Store name of any aliased branches
       */
      std::map<string,string> _aliasMap;
      /**
       * type of data linked to component names e.g._type_comps["tru"]["components_p4"] = "tru_px,tru_py,tru_pz,tru_m"
       */
      std::map<string, std::map<string,string>> _type_comps;
      /*
       * First/Primary type name
       */
      std::string _primary_type;
      
    };//ConfigReaction

  }//config
}//rad
