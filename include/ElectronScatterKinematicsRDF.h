#pragma once
#include "ElectronScatterKinematics.h"

void ElectronScatterKinematicsRDF(){}

namespace rad{
  namespace rdf{
  
    void Q2(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::Q2(%s,components_p4)",names::ReactionMap().data()));
    }
    void CosThetaCM(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::CosThetaCM(%s,components_p4)",names::ReactionMap().data()));
    }
    void PhiCM(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::PhiCM(%s,components_p4)",names::ReactionMap().data()));
    }
    void CMAngles(config::ConfigReaction& cr,const string& name){
      CosThetaCM(cr,name+"_CosTheta");
      PhiCM(cr,name+"_Phi");
    }
     void CosThetaProtonRest(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::CosThetaProtonRest(%s,components_p4)",names::ReactionMap().data()));
    }
    void PhiProtonRest(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::PhiProtonRest(%s,components_p4)",names::ReactionMap().data()));
    }


  }//rdf
}//rad
