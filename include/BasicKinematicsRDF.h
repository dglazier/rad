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
    void DeltaPhi(config::ConfigReaction& cr,const string& name, const string& idxs){
      cr.DefineForAllTypes(name, (Form("rad::DeltaPhi(%s,components_p3)",idxs.data() )) );
    }
    void DeltaTheta(config::ConfigReaction& cr,const string& name, const string& idxs){
      cr.DefineForAllTypes(name, (Form("rad::DeltaTheta(%s,components_p3)",idxs.data() )) );
    }
    void DeltaP(config::ConfigReaction& cr,const string& name, const string& idxs){
      cr.DefineForAllTypes(name, (Form("rad::DeltaP(%s,components_p3)",idxs.data() )) );
    }

    
    void PrintParticles(config::ConfigReaction& cr,const string& type=rad::names::data_type::Truth()){
    
      //Note can't sue templated lambdas so use JITing, but can't use Foreach for this
      //So just use Filter returning true.
      cr.Filter(Form("rad::PrintParticles(\"%s\",rdfentry_,%spid,%spx,%spy,%spz,%sm);return true;",type.data(),type.data(),type.data(),type.data(),type.data(),type.data()),type+"print");
      //cr.setCurrFrame(cf);
    }
 
  }//namespace rdf
}//namespace rad
