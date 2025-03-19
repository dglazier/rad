#pragma once
#include "BasicKinematics.h"
#include "ConfigReaction.h"


////Link DataFrame to Computations
namespace rad{
  namespace rdf{
 
    void Mass(config::ConfigReaction& cr,const string& name, const string& pos, const string& neg="{}"){
      //Note originally the rad:: function was called Mass
      //however this broke when a branch in the input tree
      // had a LEAF (not just branch) called Mass!!!
      // Function names need to be slightly unusual
      cr.DefineForAllTypes(name, (Form("rad::FourVectorMassCalc(%s,%s,components_p4)",pos.data(),neg.data() )) );
    }

  void PrintParticles(config::ConfigReaction& cr,const string& type="tru"){
      //cr.DefineForAllTypes(name, (Form("rad::PrintParticles(rdfentry_,components_p4)")) );
      auto cf = cr.CurrFrame();
      std::vector<std::string> cols ={"rdfentry_",type+"pid",type+"px",type+"py",type+"pz",type+"m"};
      cf.Foreach([](ULong64_t entry,const ROOT::RVecI &pid,const ROOT::RVecF &px, const ROOT::RVecF &py, const ROOT::RVecF &pz, const ROOT::RVecD &m){
	rad::PrintParticles(entry,pid,px,py,pz,m);},
	cols);
      
      cr.setCurrFrame(cf);
    }
 
  }//namespace rdf
}//namespace rad
