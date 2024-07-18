#pragma once
#include "ReactionKinematics.h"

void ReactionKinematicsRDF(){}

namespace rad{
  namespace rdf{

    void MissMass(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::MissMass(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    void MissMass2(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::MissMass2(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void TBot(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TBot(%s,components_p4)",names::ReactionMap().data()));
    }
    void TTop(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TTop(%s,components_p4)",names::ReactionMap().data()));
    }
    void TPrimeTop(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TPrimeTop(%s,components_p4)",names::ReactionMap().data()));
    }
    void TPrimeBot(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TPrimeBot(%s,components_p4)",names::ReactionMap().data()));
    }



  }//rdf
}//rad
