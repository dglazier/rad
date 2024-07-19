#pragma once

//!  Derived class to configure HepMC root files for electroproduction

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for HepMC files with fixed particle order
*/
#include "ElectroIonReaction.h"


namespace rad{
  namespace config{

    //! Class definition

    class HepMCElectro : public ElectroIonReaction {


    public:

    HepMCElectro(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={}) : ElectroIonReaction{treeName,fileNameGlob,columns} {
	
      }
    HepMCElectro(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ElectroIonReaction{treeName,filenames,columns} {

      }

     void AliasMomentumComponents(){
      AddType("mc");
      setBranchAlias("particles.momentum.m_v1","mc_px");
      setBranchAlias("particles.momentum.m_v2","mc_py");
      setBranchAlias("particles.momentum.m_v3","mc_pz");
      setBranchAlias("particles.mass","mc_m");

      DefineForAllTypes("phi", Form("rad::ThreeVectorPhi(mc_px,mc_py,mc_pz)"));
      DefineForAllTypes("theta", Form("rad::ThreeVectorTheta(mc_px,mc_py,mc_pz)"));
      DefineForAllTypes("eta", Form("rad::ThreeVectorEta(mc_px,mc_py,mc_pz)"));
      DefineForAllTypes("pmag", Form("rad::ThreeVectorMag(mc_px,mc_py,mc_pz)"));
    }

    };
  }
}
