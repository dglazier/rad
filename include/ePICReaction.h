#pragma once

//!  Derived class to configure ePIC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for ePIC files with fixed particle order
*/
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
       void AliasColumns(){
	AddType("rec");
	setBranchAlias("ReconstructedParticles.momentum.x","rec_px");
	setBranchAlias("ReconstructedParticles.momentum.y","rec_py");
	setBranchAlias("ReconstructedParticles.momentum.z","rec_pz");
	setBranchAlias("ReconstructedParticles.mass","rec_m");
	setBranchAlias("ReconstructedParticles.PDG","rec_pid");
	AddAdditionalComponents();
       }
      
      /**
       * Alias ReconstructedParticles and MCParticle columns
       */ 
     void AliasColumnsAndMC(){
	AddType("rec");
	setBranchAlias("ReconstructedParticles.momentum.x","rec_px");
	setBranchAlias("ReconstructedParticles.momentum.y","rec_py");
	setBranchAlias("ReconstructedParticles.momentum.z","rec_pz");
	setBranchAlias("ReconstructedParticles.mass","rec_m");
	setBranchAlias("ReconstructedParticles.PDG","rec_pid");
	
	AddType("tru");
	setBranchAlias("MCParticles.momentum.x","tru_px");
	setBranchAlias("MCParticles.momentum.y","tru_py");
	setBranchAlias("MCParticles.momentum.z","tru_pz");
	setBranchAlias("MCParticles.mass","tru_m");
	setBranchAlias("MCParticles.PDG","tru_pid");

	setBranchAlias("MCParticles.generatorStatus","tru_genStat");
	setBranchAlias("ReconstructedParticleAssociations.recID","rec_match_id");
	//make an branch with size of number of generator particles (status 1 or 4)
	//therefor ignoring any additional MCParticles
	Define("tru_n","rad::helpers::Count(tru_genStat,1)+rad::helpers::Count(tru_genStat,4)");
	//make an mc_match branch cut on actual generated particles (no secondaries)
	Define("tru_match_id",[](const ROOT::RVecI& stat,const ROOT::RVecU& simID){
	  ROOT::RVecU match_id;//(ngen);
	  const auto n = simID.size();
	  for(uint i=0;i<n;++i){
	    if( (stat[simID[i]] ==1) || ( stat[simID[i]]==4) ){
	      match_id.push_back(simID[i]);
	    }
	  }
	  return match_id;
	  },{"tru_genStat","ReconstructedParticleAssociations.simID"});
	

	
	if(_isRecTruSynched==false)
	  AddAdditionalComponents();
     }
      /**
       * Only alias MCParticle columns
       */ 
      void AliasColumnsMC(){
	AddType("tru");
	setBranchAlias("MCParticles.momentum.x","tru_px");
	setBranchAlias("MCParticles.momentum.y","tru_py");
	setBranchAlias("MCParticles.momentum.z","tru_pz");
	setBranchAlias("MCParticles.mass","tru_m");
	setBranchAlias("MCParticles.PDG","tru_pid");
	setBranchAlias("MCParticles.generatorStatus","tru_genStat");
    

	//alternative with edm4hep would look like
  	//setCurrFrame(CurrFrame().Define("pmag", getMag<edm4hep::MCParticleData> ,{"MCParticles"}));

	AddAdditionalComponents();
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
      void AliasColumnsAndMatchWithMC(bool with_removal=true){

	/*Remake this function, probably using string define
	  currently OK for rec, but need to truncate mc columns to tru_n entries*/
	//for tru : RedefineViaAlias(alias,"rad::helpers::Truncate(tru_n))
	//for rec : RedefineViaAlias(alias,Form("helpers::Reorder(%s,rec_match_id,tru_match_id,tru_n)",alias.data());
	
	_isRecTruSynched=true;
	AliasColumnsAndMC();
	
	for(const auto& col : AliasMap()){
	  cout<<"synch "<<col.first<<endl;
	  const auto& alias = col.first;
	  if(alias=="rec_match_id") continue;
	  if(alias=="tru_match_id") continue;
	  std::string match_id = TString(alias).Contains("rec") ? "rec_match_id" : "tru_match_id";
	  /*  if(TString(alias).Contains("rec")){
	    RedefineViaAlias(alias,std::string_view(Form("helpers::Reorder(%s,rec_match_id,tru_match_id,tru_n)",alias.data())));
	  }
	  else{
	    RedefineViaAlias(alias,std::string_view(Form("helpers::Truncate(%s,tru_n)",alias.data())));
	    }*/
	  //Note expression function does not work due to some epic branches having
	  // a . in the name, so must do this ugly check all types method
	    switch(static_cast<int>(DeduceColumnVectorType(col.second))) {
	    //    enum class ColType{Undef,Int,UInt,Float,Double,Short,Bool,Long};
	    case static_cast<int>(ColType::Undef):
	    std::cerr<<"Warning ePICReaction::AliasColumnsAndMatchWithMC : cannot deduce type for column "<<col.second<<" aliased to "<<alias<<std::endl;
	      break;
	      case static_cast<int>(ColType::Int):
	      if(with_removal)RedefineViaAlias(alias,helpers::Rearrange<Int_t>,{alias.data(),match_id.data()});
	      else if(match_id=="rec_match_id") RedefineViaAlias(alias,helpers::Reorder<Int_t>,{alias.data(),"rec_match_id","tru_match_id","tru_n"});
	      break;
	      case static_cast<int>(ColType::UInt):
	      if(with_removal)RedefineViaAlias(alias,helpers::Rearrange<UInt_t>,{alias.data(),match_id.data()});
	      else if(match_id=="rec_match_id")  RedefineViaAlias(alias,helpers::Reorder<UInt_t>,{alias.data(),"rec_match_id","tru_match_id","tru_n"});
	      break;
	      case static_cast<int>(ColType::Float):
	      if(with_removal)RedefineViaAlias(alias,helpers::Rearrange<Float_t>,{alias.data(),match_id.data()});
	      else if(match_id=="rec_match_id"){
		if(TString(alias).Contains("rec")){
		  //for rec vectors we need to synchornise with tru particles up to tru_n i.e. actual generated particles
		  RedefineViaAlias(alias,helpers::Reorder<Float_t>,{alias.data(),"rec_match_id","tru_match_id","tru_n"});
		}		
	      }
	      else{
		//for tru particles we need to remove secondaries, i.e. truncate at number of generated
		//due to large number of secondaries files can get very large wihtout truncation
		RedefineViaAlias(alias,helpers::Truncate<Float_t>,{alias.data(),"tru_n"});
	      }
	      
	      
	      break;
	      case static_cast<int>(ColType::Double):
	      if(with_removal)RedefineViaAlias(alias,helpers::Rearrange<Double_t>,{alias.data(),match_id.data()});
	      else if(match_id=="rec_match_id")  RedefineViaAlias(alias,helpers::Reorder<Double_t>,{alias.data(),"rec_match_id","tru_match_id","tru_n"});
	      break;
	      case static_cast<int>(ColType::Short):
	      if(with_removal)RedefineViaAlias(alias,helpers::Rearrange<Short_t>,{alias.data(),match_id.data()});
	      else if(match_id=="rec_match_id")  RedefineViaAlias(alias,helpers::Reorder<Short_t>,{alias.data(),"rec_match_id","tru_match_id","tru_n"});
	      break;
	      case static_cast<int>(ColType::Bool):
	      if(with_removal)RedefineViaAlias(alias,helpers::Rearrange<Bool_t>,{alias.data(),match_id.data()});
	      else if(match_id=="rec_match_id")  RedefineViaAlias(alias,helpers::Reorder<Bool_t>,{alias.data(),"rec_match_id","tru_match_id","tru_n"});
	      break;
	      case static_cast<int>(ColType::Long):
	      if(with_removal)RedefineViaAlias(alias,helpers::Rearrange<Long_t>,{alias.data(),match_id.data()});
	      else if(match_id=="rec_match_id")  RedefineViaAlias(alias,helpers::Reorder<Long_t>,{alias.data(),"rec_match_id","tru_match_id","tru_n"});
	      break;
	    
	      default:
	      break;
	    
	      }
	  
	}//for AliasMap

	//need to do it here or these calculations are done before synching
	//which will not work as different number of elements
	AddAdditionalComponents();
      }//AliasColumnsAndMatchWithMC

      void AddAdditionalComponents(){
	//and add some additional columns
	DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(components_p3)"));
	DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(components_p3)"));
	DefineForAllTypes("eta", Form("rad::ThreeVectorEta(components_p3)"));
	DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(components_p3)"));
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
       * Case Reconstructed and truth not synched and need associations
       */
       void Resolution(const string& var){
	 if(_isRecTruSynched) return ResolutionSynched(var);
	 
	Define(string("res_")+var,[](const ROOT::RVecU &irec,const ROOT::RVecU &itru,const ROOT::RVecF &rec,const ROOT::RVecF &tru){
	  //create array with values of rec, indiced by irec=rec_match_id
	  auto rec_matched = ROOT::VecOps::Take(rec,irec);
	  //create array with values of tru, indiced by isim=sim_match_id
	  auto tru_matched = ROOT::VecOps::Take(tru,itru);
	  return (rec_matched - tru_matched);
	},{"rec_match_id","tru_match_id",string("rec_")+var,string("tru_")+var});
      }
      void ResolutionFraction(const string& var){
	 if(_isRecTruSynched) return ResolutionFractionSynched(var);
	Define(string("res_")+var,[](const ROOT::RVecU &irec,const ROOT::RVecU &itru,const ROOT::RVecF &rec,const ROOT::RVecF &tru){
	  //create array with values of rec, indiced by irec=rec_match_id
	  auto rec_matched = ROOT::VecOps::Take(rec,irec);
	  //create array with values of tru, indiced by isim=sim_match_id
	  auto tru_matched = ROOT::VecOps::Take(tru,itru);
	  return (rec_matched - tru_matched)/tru_matched;
	},{"rec_match_id","tru_match_id",string("rec_")+var,string("tru_")+var});
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
      bool _isRecTruSynched=false; //should remove this and have a derived class...
   
    };//ePICReaction

    
  }//config
}//rad
