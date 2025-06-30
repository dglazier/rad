#pragma once
#include "ParticleGenerator.h"
#include "StringUtilities.h"

//void ParticleGeneratorRDF(){}

namespace rad{
  namespace generator{
    
    void GenerateTwoBody(config::ConfigReaction& cr, const std::vector<std::string> &names, const string &parent){
      cr.DefineForAllTypes( "twobody_parent",Form("rad::generator::TwoBody(%s,%s,components_p4)",rad::utils::combineVectorToString(names).data() ,parent.data() ));
    }
    
  }
}
