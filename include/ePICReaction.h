#pragma once

//!  Derived class to configure ePIC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for ePIC files with fixed particle order
*/
#include "ElectroIonReaction.h"
#include "ElectroIonReaction.h"

namespace rad{
  namespace config{

    //! Class definition

    class ePICReaction : public ElectroIonReaction {


    public:

      ePICReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ElectroIonReaction{treeName,fileNameGlob,columns} {

      }
     ePICReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ElectroIonReaction{treeName,filenames,columns} {

      }

      /**
       * Only alias ReconstructedParticles columns
       */ 
       void AliasColumns(Bool_t IsEnd=kTRUE){
	AddType("rec_");
	setBranchAlias("ReconstructedParticles.momentum.x","rec_px");
	setBranchAlias("ReconstructedParticles.momentum.y","rec_py");
	setBranchAlias("ReconstructedParticles.momentum.z","rec_pz");
	setBranchAlias("ReconstructedParticles.mass","rec_m");
	setBranchAlias("ReconstructedParticles.PDG","rec_pid");
	
	 if(IsEnd) AddAdditionalComponents();
       }
      
      /**
       * Alias ReconstructedParticles and MCParticle columns
       */ 
     void AliasColumnsAndMC(Bool_t IsEnd=kTRUE){
	AliasColumns(kFALSE);
	AliasColumnsMC(kFALSE);
	/*
	AddType("rec_");
	setBranchAlias("ReconstructedParticles.momentum.x","rec_px");
	setBranchAlias("ReconstructedParticles.momentum.y","rec_py");
	setBranchAlias("ReconstructedParticles.momentum.z","rec_pz");
	setBranchAlias("ReconstructedParticles.mass","rec_m");
	setBranchAlias("ReconstructedParticles.PDG","rec_pid");
	*/
	/*
	AddType("tru_");
	setBranchAlias("MCParticles.momentum.x","tru_px");
	setBranchAlias("MCParticles.momentum.y","tru_py");
	setBranchAlias("MCParticles.momentum.z","tru_pz");
	setBranchAlias("MCParticles.mass","tru_m");
	setBranchAlias("MCParticles.PDG","tru_pid");

	setBranchAlias("MCParticles.generatorStatus","tru_genStat");
	*/
	
	//remove all but true generated beam+final state particles
	//rec_match_id : 0,1,2,...N=number beam+final particles 
	Define("rec_match_id",[](const ROOT::RVecI& stat){

	  auto filtstat = stat[stat==1||stat==4];
	  auto id = helpers::Enumerate<uint>(filtstat.size());
	  return id;//[filtstat==1];
	},{"tru_genStat"}); //just need tru gen status to get N
	
	//make an mc_match branch cut on actual generated particles (no secondaries)
	//Points rec array to tru array. rec_array has no beam particles, so can ignore
	Define("tru_match_id",[](const ROOT::RVecI& stat,const ROOT::RVecU& simID){
	  
	  ROOT::RVecU match_id;
	  const auto n = simID.size();
	  for(uint i=0;i<n;++i){
	    if( (stat[simID[i]] == 1)  ){
	      match_id.push_back(simID[i]);
	    }
	  }
	  return match_id;
	},{"tru_genStat","ReconstructedParticleAssociations.simID"});//simID points from rec to tru
	
	
	//make an branch with size of number of generator particles (status 1 or 4)
	//used to truncate tru arrays
	Define("tru_n","rad::helpers::Count(tru_genStat,1)+rad::helpers::Count(tru_genStat,4)");
	
	
	if(IsEnd) AddAdditionalComponents();
     }
      /**
       * Only alias MCParticle columns
       */ 
      void AliasColumnsMC(Bool_t IsEnd=kTRUE){
	AddType("tru_");
	setBranchAlias("MCParticles.momentum.x","tru_px");
	setBranchAlias("MCParticles.momentum.y","tru_py");
	setBranchAlias("MCParticles.momentum.z","tru_pz");
	setBranchAlias("MCParticles.mass","tru_m");
	setBranchAlias("MCParticles.PDG","tru_pid");
	setBranchAlias("MCParticles.generatorStatus","tru_genStat");
    

	 if(IsEnd) AddAdditionalComponents();
      }
      /**
       * Alias the columns and rearrange entries 
       * according to ReconstructedParticleAssociations
       * this reorders reconstructed to match mc
       * Note, we should consider changing to the opposite
       * As this is slower and produces larger output trees
       * than matching to reconstructed order.
       * The advantage this way is that the particles are
       * well defined.
       */

      //this function is very ugly at the moment. This is due to template type requiring
      //explicit type when call the Re* functions.
      void AliasColumnsAndMatchWithMC(Bool_t IsEnd=kTRUE){

	/*Remake this function, probably using string define
	  currently OK for rec, but need to truncate mc columns to tru_n entries*/
	//for tru : RedefineViaAlias(alias,"rad::helpers::Truncate(tru_n))
	//for rec : RedefineViaAlias(alias,Form("helpers::Reorder(%s,rec_match_id,tru_match_id,tru_n)",alias.data());
	
	AliasColumnsAndMC(kFALSE);

	for(const auto& col : AliasMap()){
	  cout<<"synch "<<col.first<<endl;
	  const auto& alias = col.first;
	  
	  switch(static_cast<int>(DeduceColumnVectorType(col.second )) ) {
	    
	  case static_cast<int>(ColType::Undef):
	    break;
	  case static_cast<int>(ColType::UInt):
	    RedefineFundamental<UInt_t>(alias);
	    break;
	  case static_cast<int>(ColType::Int):
	    RedefineFundamental<Int_t>(alias);
	    break;
	  case static_cast<int>(ColType::Float):
	    RedefineFundamental<Float_t>(alias);
	    break;
	  case static_cast<int>(ColType::Double):
	    RedefineFundamental<Double_t>(alias);
	    break;
	  case static_cast<int>(ColType::Short):
	    RedefineFundamental<Short_t>(alias);
	    break;
	  case static_cast<int>(ColType::Bool):
	    RedefineFundamental<Bool_t>(alias);
	    break;
	  case static_cast<int>(ColType::Long):
	    RedefineFundamental<Long_t>(alias);
	    break;
	    
	  default:
	    break;
	  }
	}
	return;


	//need to do it here or these calculations are done before synching
	//which will not work as different number of elements
	 if(IsEnd) AddAdditionalComponents();
      }//AliasColumnsAndMatchWithMC


      void AddAdditionalComponents(){
	//and add some additional columns
	DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
	DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
	DefineForAllTypes("eta", Form("rad::ThreeVectorEta(components_p3)"));
	DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
      }
      
      template<typename T> 
	void RedefineFundamental( const string& name ){
	
	auto contains = [](const std::string&  s1,const std::string& s2){
	  return (s1.find(s2) != std::string::npos);
	};
	
	if(contains(name,"rec") ){
	    RedefineViaAlias(name,helpers::Reorder<T>,{name.data(),"rec_match_id","tru_match_id","tru_n"});
	  }
	else if(contains(name,"tru") ){
	  RedefineViaAlias(name,helpers::Truncate<T>,{name.data(),"tru_n"});
	  
	}
	
      }
      

      
      /**
       * calculate the difference in reconsutructed and truth variables
       */
      /**
       * Case Reconstructed and truth synched via AliasColumnsAndMatchWithMC()
       */
      void ResolutionSynched(const string& var){
	Define(string("res_")+var,[](const ROOT::RVecF &rec,const ROOT::RVecF &tru){
	  return (rec - tru);
	},{string("rec_")+var,string("tru_")+var});
      }
      void ResolutionFractionSynched(const string& var){
	Define(string("res_")+var,[](const ROOT::RVecF &rec,const ROOT::RVecF &tru){
	  return (rec - tru)/tru;
	},{string("rec_")+var,string("tru_")+var});
      }
      
      /**
       * Mask tracks that do not have a valid ReconstructedParticleAssociations
       */
      void MaskMCMatch(){
	//define function to create and apply mask for rec or tru
	auto match = [this](const string& type){
	  //remove Pids that are not matched
	  ApplyPidMask(type+"mask_mcmatch",type,
		       [](const ROOT::RVecU& matchid,const ROOT::RVecI& pid){
			 //create a mask array with indices with matchid values==1
			 RVecI mask(pid.size());
			 //Fill mask by looping over matched indexes and checking if
			 //particle exists in recID
			 for(auto& idx: matchid){
			   mask[idx]=1; //this will add a 1 at mask[idx]
			 }
			 return mask;
		       },{type+"_match_id",type+"_pid"} );
	};//end of match lambda
	
	match("rec");
	match("tru");
      }
      
    private:
      
    };//ePICReaction

    
  }//config
}//rad
