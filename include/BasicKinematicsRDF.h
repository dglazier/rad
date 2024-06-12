#pragma once
#include "BasicKinematics.h"
#include "ConfigReaction.h"


////Link DataFrame to Computations
namespace rad{
  namespace rdf{
 
    void Mass(config::ConfigReaction& cr,const string_view name, const string_view pos, const string_view neg="{}"){
      cr.DefineForAllTypes(name, (Form("rad::Mass(%s,%s,components_p4)",pos.data(),neg.data() )) );
    }
 
  }//namespace rdf
}//namespace rad
