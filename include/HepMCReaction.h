#pragma once

//!  Derived class to configure HepMC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for HepMC files with fixed particle order
*/
#include "ConfigReaction.h"


namespace rad{
  namespace config{

  //! Class definition

    class HepMCReaction : public ConfigReaction {


  public:

  HepMCReaction(const std::string_view treeName, const std::string_view fileNameGlob) : ConfigReaction{treeName,fileNameGlob} {

    }

    void AliasMomentumComponents(){
      AddType("mc")
      setBranchAlias("particles.momentum.m_v1","mc_px");
      setBranchAlias("particles.momentum.m_v2","mc_py");
      setBranchAlias("particles.momentum.m_v3","mc_pz");
      setBranchAlias("particles.mass","mc_m");

      DefineForAllTypes("pmag", Form("compute::VecMag(components_p3)"));
    }


    };//HepMCReaction
    
  }//config
}//rad
