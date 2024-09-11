#pragma once

//!  Derived class to configure ePIC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for Electron scattering reactions
  Derive from this if your experiment is electron scattering
*/
#include "ConfigReaction.h"
#include "RVecHelpers.h"
#include "BasicKinematics.h"
#include "Beams.h"

namespace rad{
  namespace config{
    

    //! Class definition

    class ElectroIonReaction : public ConfigReaction {


    public:

      ElectroIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ConfigReaction{treeName,fileNameGlob,columns} {

      }
     ElectroIonReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ConfigReaction{treeName,filenames,columns} {

      }

      /**
       * Make map that links particle names to indices in user functions
       * in C++ functions you can use the RVecIndexMap object indexed by 
       * name of the reaction component you need
       */
      void makeParticleMap() override {
	//note, ordering in arguments, map and names must be maintained
  	Define(names::ReactionMap().data(),
	       [](const int& beamion, const int& beamel,
		  const RVecI& baryons,const RVecI& mesons,const int& scatel){
		 return RVecIndexMap{{beamion},{beamel},baryons,mesons,{scatel}};},
	       {names::BeamIon().data(),names::BeamEle().data(),
		names::Baryons().data(),names::Mesons().data(),names::ScatEle().data()});
      }
      // /**
      //  * Make map that links particle names to indices in user functions
      //  * in C++ functions you can use the RVecIndexMap object indexed by 
      //  * name of the reaction component you need
      //  */
      // void makeBeamIndices() override {
      // 	//note, ordering in arguments, map and names must be maintained
      // 	Define(names::ReactionMap().data(),
      // 	       [](const int& beamion, const int& beamel){
      // 		 return ROOT::RVecI{{beamion},{beamel} } };
      // 	       );
      // }
      /** 
       * Set constant index for beam electron, scattered electron and beam ion
       * This assumes constant position in collection (e.g in some HepMC3 files)
       * and update the current frame to the aliased one
       */
      void setBeamElectronIndex(const int idx){
	setParticleIndex(names::BeamEle().data(),idx);
      }
      void setScatElectronIndex(const int idx){
	setParticleIndex(names::ScatEle().data(),idx);
      }
      template<typename Lambda>
      void setScatElectronIndex(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {}){
	 setParticleIndex(names::ScatEle().data(),func,columns);
      }
      void setBeamIonIndex(const int idx){
	setParticleIndex(names::BeamIon().data(),idx);
      }
      /**
       * Collect variable indices for scattered electron
       */
       template<typename Lambda>
      void setScatElectron(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {} ){
	setParticleIndex(names::ScatEle(),func, columns);
      }
      
    private:
      

    };//ElectroIonReaction

  }//config
  namespace electroion{
    /**
     * virtual photon 4-vector
     */
    template<typename Tp, typename Tm>
    PxPyPzMVector PhotoFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      
      auto phot =  beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      SubtractFourVector(phot,react[names::ScatEleIdx()],px,py,pz,m);
      return phot;
      
    }
    const std::string_view  BeamIndices() {return Form("%s,%s",names::BeamIon().data(),names::BeamEle().data()); }//"beam_ele,beam_ion";}

  }
  
}//rad

//Declare we are using this PhotoFourVector in kinematics
using rad::electroion::PhotoFourVector;
using rad::electroion::BeamIndices;
