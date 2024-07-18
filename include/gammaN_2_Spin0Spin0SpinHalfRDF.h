#pragma once
#include "gammaN_2_Spin0Spin0SpinHalf.h"


namespace rad{
  namespace rdf{
    
    void PhotoCMAngles(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::gn2s0s0s12::PhotoCMDecay(%s,components_p4)",names::ReactionMap().data()));
    }

    void HelicityAngles(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::gn2s0s0s12::PhotoHelicityDecay(%s,components_p4)",names::ReactionMap().data()));
    }

    void GJAngles(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::gn2s0s0s12::PhotoGJDecay(%s,components_p4)",names::ReactionMap().data()));
    }

    
  }
}
