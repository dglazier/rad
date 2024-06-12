#pragma once
#include "ElectronScatterKinematics.h"

void ElectronScatterKinematicsRDF(){}

namespace rad{
  namespace rdf{

    // void CosThetaCM(config::ConfigReaction& cr , const string& name){
    //   cr.DefineForAllTypes(name+"_CosTheta", []( const rad::ElCMDecay_t& cm){ return cm.CosTheta;} , {name});
    // }
    // void PhiCM(config::ConfigReaction& cr,const string& name){
    //   cr.DefineForAllTypes(name+"_Phi", []( const rad::ElCMDecay_t& cm){ return cm.Phi;} , {name});
    // }
    
    void CosThetaCM(config::ConfigReaction& cr,const string_view& name){
      cr.DefineForAllTypes(name, Form("rad::CosThetaCM(%s,components_p4)",names::ReactionMap().data()));
    }
    void PhiCM(config::ConfigReaction& cr,const string_view& name){
      cr.DefineForAllTypes(name, Form("rad::PhiCM(%s,components_p4)",names::ReactionMap().data()));
    }
    void CMAngles(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::ElectroCMDecay(%s,components_p4)",names::ReactionMap().data()));
      CosThetaCM(cr,name+"_CosTheta");
      PhiCM(cr,name+"_Phi");
      //CosThetaCM(cr,name);
      // PhiCM(cr,name);
    }

  }//rdf
}//rad
