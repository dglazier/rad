#pragma once
#include "ReactionKinematics.h"

void ReactionKinematicsRDF(){}

namespace rad{
  namespace rdf{


  void MissMass(config::ConfigReaction& cr,const string& name, const string_view neg){
    cr.DefineForAllTypes(name, Form("rad::MissMass(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
   }
  void TBot(config::ConfigReaction& cr,const string& name){
    cr.DefineForAllTypes(name, Form("rad::TBot(%s,components_p4)",names::ReactionMap().data()));
  }
  void TTop(config::ConfigReaction& cr,const string& name){
    cr.DefineForAllTypes(name, Form("rad::Top(%s,components_p4)",names::ReactionMap().data()));
  }
  void TPrimeTop(config::ConfigReaction& cr,const string& name){
    cr.DefineForAllTypes(name, Form("rad::TPrimeTop(%s,components_p4)",names::ReactionMap().data()));
  }
  void TPrime(config::ConfigReaction& cr,const string& name){
    cr.DefineForAllTypes(name, Form("rad::TPrimeBot(%s,components_p4)",names::ReactionMap().data()));
  }



  }//rdf
}//rad
