#pragma once

//!  Derived class to configure ePIC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for ePIC files with fixed particle order
*/
#include "ElectroIonReaction.h"
#include "ePICUtilities.h"

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
	
	if(IsEnd){
	  /////RedefineFundamentalAliases();

	  rad::epic::UndoAfterBurn undoAB{_p4ion_beam,_p4el_beam,-0.025};
	  auto undoAB_rec=[undoAB](ROOT::RVecF &px,ROOT::RVecF &py,ROOT::RVecF&pz, const ROOT::RVecF &m){
	    undoAB(px,py,pz,m);
	    return px;
	  };
	  
	  
	  //Undo the afterburner procedure
	  //here we just account for crossing angle
	  //just need to redefine 1 component. Other 2 have been updated
	  //need to redefine at least one to make sure this is called before
	  //any of the components
	  RedefineViaAlias("rec_px", undoAB_rec, {"rec_px","rec_py","rec_pz","rec_m"});
	  
	  //undo afterburn on beam components
	  auto beams_px = ROOT::RVecD{_p4el_beam.X(),_p4ion_beam.X()};
	  auto beams_py = ROOT::RVecD{_p4el_beam.Y(),_p4ion_beam.Y()};
	  auto beams_pz = ROOT::RVecD{_p4el_beam.Z(),_p4ion_beam.Z()};
	  auto beams_m = ROOT::RVecD{_p4el_beam.M(),_p4ion_beam.M()};
	  undoAB(beams_px,beams_py,beams_pz,beams_m);
	  _p4el_beam.SetCoordinates(beams_px[0],beams_py[0],beams_pz[0],beams_m[0]);
	  _p4ion_beam.SetCoordinates(beams_px[1],beams_py[1],beams_pz[1],beams_m[1]);

	  cout<<"Undo afterburn head on beam 4-vectors : "<< _p4el_beam<<" "<<_p4ion_beam<<endl;
	  
	  DefineBeamElectron();
	  DefineBeamIon();

	  AddAdditionalComponents();
	}
       }
      
      /**
       * Alias ReconstructedParticles and MCParticle columns
       */ 
     void AliasColumnsAndMC(Bool_t IsEnd=kTRUE){
	AliasColumns(kFALSE);
	AliasColumnsMC(kFALSE);
	//Matching reconstructed to truth :
	// 1) Find final state truth tru_genStat == 1 =>rec_match_id 0,,N
	// 2) Map from old to new : final_match_id
	//    final_match_id = sizeof(MCParticles)
	//             value = order in new arrays or -1 if not included
	// 3) Create new simID : final_match_id[simID[]]
	//    converts to position in new truth array
	// 4) Create new arrays : sizeof(tru_final_state)
	//    tru_ entries = all filled as truth
	//    rec_ entries = filled if that particle was reconstructed
	
	//remove all but true generated beam+final state particles
	//rec_match_id : 0,1,2,...N=number beam+final particles 
	Define("rec_match_id",[](const ROOT::RVecI& stat){

	  auto filtstat = stat[stat==1];
	  auto id = helpers::Enumerate<uint>(filtstat.size());
	  return id;//[filtstat==1];
	},{"tru_genStat"}); //just need tru gen status to get N
	
	//make an mc_match branch cut on actual generated particles (no secondaries)
	//Points rec array to tru array. rec_array has no beam particles, so can ignore
	Define("tru_match_id",[](const ROOT::RVecI& stat,const ROOT::RVecU& simID,const ROOT::RVecU& finalID){
	  const auto n = finalID.size(); //mcparticles stat==1
	  ROOT::RVecI final_match_id(n,-1);
	  for(uint i=0;i<n;++i){
	    if(i>=simID.size())break;
	    //if this truth particle was reconstructed add its new id
	    // if(rad::helpers::Contains(simID,finalID[i]))
	    //final_match_id[finalID[i]]=i;

	    if(rad::helpers::Contains(finalID,simID[i])){
	      //final_match_id[finalID[i]]=simID[i]-2;
	      final_match_id[i]=rad::helpers::findIndex(finalID,simID[i]);
	    }
	  }
	  ROOT::RVecU tru_match_id =final_match_id[final_match_id!=-1]; //Filter valid ids
	  //std::cout<<"tru_match_id "<<simID<<" "<<simID.size()<<" "<<finalID<<" "<<finalID.size()<<" "<<tru_match_id<<" "<<tru_match_id.size()<<std::endl;
	  return tru_match_id;
	  
	},{"tru_genStat","ReconstructedParticleAssociations.simID","tru_final_id"});//simID points from rec to tru
	

	//make an branch with size of number of generator particles (status 1 or 4)
	//used to truncate tru arrays
	//Define("tru_n","rad::helpers::Count(tru_genStat,1)+rad::helpers::Count(tru_genStat,4)");
	Define("tru_n","rad::helpers::Count(tru_genStat,1)");
	
	
	if(IsEnd){
	  RedefineFundamentalAliases();

	  rad::epic::UndoAfterBurn undoAB{ _p4ion_beam,_p4el_beam,-0.025};
	  auto undoAB_tru=[undoAB](ROOT::RVecF &px,ROOT::RVecF &py,ROOT::RVecF&pz, const ROOT::RVecD &m){
	    undoAB(px,py,pz,m);
	    return px;
	  };
	  //Undo the afterburner procedure
	  //here we just account for crossing angle
	  //just need to redefine 1 component. Other 2 have been updated
	  //need to redefine at least one to make sure this is called before
	  //any of the components
	  RedefineViaAlias("tru_px", undoAB_tru, {"tru_px","tru_py","tru_pz","tru_m"});

	  	  //undo afterburn on beam components
	  auto beams_px = ROOT::RVecD{_p4el_beam.X(),_p4ion_beam.X()};
	  auto beams_py = ROOT::RVecD{_p4el_beam.Y(),_p4ion_beam.Y()};
	  auto beams_pz = ROOT::RVecD{_p4el_beam.Z(),_p4ion_beam.Z()};
	  auto beams_m = ROOT::RVecD{_p4el_beam.M(),_p4ion_beam.M()};
	  undoAB(beams_px,beams_py,beams_pz,beams_m);
	  _p4el_beam.SetCoordinates(beams_px[0],beams_py[0],beams_pz[0],beams_m[0]);
	  _p4ion_beam.SetCoordinates(beams_px[1],beams_py[1],beams_pz[1],beams_m[1]);

	  cout<<"Undo afterburn head on beam 4-vectors : "<< _p4el_beam<<" "<<_p4ion_beam<<endl;
	  
	  DefineBeamElectron();
	  DefineBeamIon();

	  //after undo afterburn
	  AddAdditionalComponents();
	  //needed to make sure tru and rec are defined at the same time
	  //and therefore contain same elements.
	  Filter([](const ROOT::RVecF& th,const ROOT::RVecF& pmag,const ROOT::RVecF& ph,const ROOT::RVecF& eta,const ROOT::RVecF& rth,const ROOT::RVecF& rpmag,const ROOT::RVecF& rph,const ROOT::RVecF& reta){return true;},{"tru_pmag","tru_theta","tru_phi","tru_eta","rec_pmag","rec_theta","rec_phi","rec_eta"});
	}
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
    
	Define("tru_final_id",[](const ROOT::RVecI& stat){
	  //map full array to final state only array  
	  auto indices = helpers::Enumerate<uint>(stat.size());
	  return indices[stat==1];

	},{"tru_genStat"});//simID points from rec to tru


	
	if(IsEnd){
	  RedefineFundamentalAliases();

	  rad::epic::UndoAfterBurn undoAB{ _p4ion_beam,_p4el_beam,-0.025};
	  auto undoAB_tru=[undoAB](ROOT::RVecF &px,ROOT::RVecF &py,ROOT::RVecF&pz, const ROOT::RVecD &m){
	    undoAB(px,py,pz,m);
	    return px;
	  };
	  //Undo the afterburner procedure
	  //here we just account for crossing angle
	  //just need to redefine 1 component. Other 2 have been updated
	  //need to redefine at least one to make sure this is called before
	  //any of the components
	  RedefineViaAlias("tru_px", undoAB_tru, {"tru_px","tru_py","tru_pz","tru_m"});
	  //undo afterburn on beam components
	  auto beams_px = ROOT::RVecD{_p4el_beam.X(),_p4ion_beam.X()};
	  auto beams_py = ROOT::RVecD{_p4el_beam.Y(),_p4ion_beam.Y()};
	  auto beams_pz = ROOT::RVecD{_p4el_beam.Z(),_p4ion_beam.Z()};
	  auto beams_m = ROOT::RVecD{_p4el_beam.M(),_p4ion_beam.M()};
	  undoAB(beams_px,beams_py,beams_pz,beams_m);
	  _p4el_beam.SetCoordinates(beams_px[0],beams_py[0],beams_pz[0],beams_m[0]);
	  _p4ion_beam.SetCoordinates(beams_px[1],beams_py[1],beams_pz[1],beams_m[1]);

	  cout<<"Undo afterburn head on beam 4-vectors : "<< _p4el_beam<<" "<<_p4ion_beam<<endl;
	  
	  DefineBeamElectron();
	  DefineBeamIon();

	  //after undo afterburn
	  AddAdditionalComponents();
 
	}
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


	//need to do it here or these calculations are done before synching
	//which will not work as different number of elements
	
	if(IsEnd){
	  
	  RedefineFundamentalAliases();
	  //use true masses for each rec track
	  //i.e. ignore PID
	  //Must truncate to make sure return array is same size as in array
	  //RDF may add beams to tru_m before calling this giving a size mismatch
	  RedefineViaAlias("rec_m",[](const ROOT::RVecF& recm,const ROOT::RVecD& trum){return helpers::Truncate(ROOT::RVecF(trum),recm.size());},{"rec_m","tru_m"});
	  
	  rad::epic::UndoAfterBurn undoAB{_p4ion_beam,_p4el_beam,-0.025};
	  auto undoAB_tru=[undoAB](ROOT::RVecF &px,ROOT::RVecF &py,ROOT::RVecF&pz, const ROOT::RVecD &m){
	    undoAB(px,py,pz,m);
	    return px;
	  };
	  auto undoAB_rec=[undoAB](ROOT::RVecF &px,ROOT::RVecF &py,ROOT::RVecF&pz, const ROOT::RVecF &m){
	    undoAB(px,py,pz,m);

	    return px;
	  };
	  
	  
	  //Undo the afterburner procedure
	  //here we just account for crossing angle
	  //just need to redefine 1 component. Other 2 have been updated
	  //need to redefine at least one to make sure this is called before
	  //any of the components
	  RedefineViaAlias("tru_px", undoAB_tru, {"tru_px","tru_py","tru_pz","tru_m"});
	  RedefineViaAlias("rec_px", undoAB_rec, {"rec_px","rec_py","rec_pz","rec_m"});
	  
	  //undo afterburn on beam components
	  auto beams_px = ROOT::RVecD{_p4el_beam.X(),_p4ion_beam.X()};
	  auto beams_py = ROOT::RVecD{_p4el_beam.Y(),_p4ion_beam.Y()};
	  auto beams_pz = ROOT::RVecD{_p4el_beam.Z(),_p4ion_beam.Z()};
	  auto beams_m = ROOT::RVecD{_p4el_beam.M(),_p4ion_beam.M()};
	  undoAB(beams_px,beams_py,beams_pz,beams_m);
	  _p4el_beam.SetCoordinates(beams_px[0],beams_py[0],beams_pz[0],beams_m[0]);
	  _p4ion_beam.SetCoordinates(beams_px[1],beams_py[1],beams_pz[1],beams_m[1]);

	  cout<<"Undo afterburn head on beam 4-vectors : "<< _p4el_beam<<" "<<_p4ion_beam<<endl;
	  
	  DefineBeamElectron();
	  DefineBeamIon();
	
	  AddAdditionalComponents();
	  //add resolution functions
	  ResolutionFraction<float>("pmag");
	  Resolution("theta");
	  Resolution("phi");
	  Resolution("eta");

	  //needed to make sure tru and rec are defined at the same time
	  //and therefore contain same elements.
	  Filter([](const ROOT::RVecF& th,const ROOT::RVecF& pmag,const ROOT::RVecF& ph,const ROOT::RVecF& eta,const ROOT::RVecF& rth,const ROOT::RVecF& rpmag,const ROOT::RVecF& rph,const ROOT::RVecF& reta){return true;},{"tru_pmag","tru_theta","tru_phi","tru_eta","rec_pmag","rec_theta","rec_phi","rec_eta"});
 
	}
      }//AliasColumnsAndMatchWithMC


      void RedefineFundamentalAliases(){

	for(const auto& col : AliasMap()){
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

    }
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
	  RedefineViaAlias(name,helpers::Rearrange<T>,{name.data(),"tru_final_id"});
	  
	}
	
      }
      

      
      /**
       * calculate the difference in reconsutructed and truth variables
       * Case Reconstructed and truth synched via AliasColumnsAndMatchWithMC()
       */
      //template<typename T>
      void Resolution(const string& var){
	// Define(string("res_")+var,[](const ROOT::RVec<T> &rec,const ROOT::RVec<T> &tru){
	//   return (rec - tru);
	// },{string("rec_")+var,string("tru_")+var});
	Define(string("res_")+var,Form("%s-%s",(string("tru_")+var).data(),(string("rec_")+var).data() ));
      }
      template<typename T>
      void ResolutionFraction(const string& var){
	Define(string("res_")+var,[](const ROOT::RVec<T> &rec,const ROOT::RVec<T> &tru){
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

      void SetBeamsFromMC(Long64_t nrows=100){
	auto nthreads  = ROOT::GetThreadPoolSize();
	//Range only works in single thread mode
	if(nthreads) ROOT::DisableImplicitMT();

	auto tempframe = GetFileNames().size()==0 ?
	  ROOT::RDataFrame{GetTreeName(),GetFileName()} :
	  ROOT::RDataFrame{GetTreeName(),GetFileNames()};
	
	auto beamdf = tempframe.Range(nrows).Define("emean","MCParticles.momentum.z[MCBeamElectrons_objIdx.index[0]]").Define("pzmean","MCParticles.momentum.z[MCBeamProtons_objIdx.index[0]]").Define("pxmean","MCParticles.momentum.x[MCBeamProtons_objIdx.index[0]]");
	auto pze  = beamdf.Mean("emean");
	auto pzp  = beamdf.Mean("pzmean");
	auto pxp  = beamdf.Mean("pxmean");

	setBeamElectron(0,0,*pze);//afterburned number
	setBeamIon(*pxp,0,*pzp);//afterburned number

	if(nthreads) ROOT::EnableImplicitMT(nthreads);
      }
      
    private:
      
    };//ePICReaction

    
  }//config
}//rad
