#pragma once
#include "ConfigReaction.h"

namespace rad{
  namespace reaction{

    void CountParticles(rad::config::ConfigReaction* rad, const std::string& type){
      
	rad->Define(type+"Ngamma",Form("rad::helpers::Count(%spid,22)",type.data()) );
	rad->Define(type+"Npip",Form("rad::helpers::Count(%spid,211)",type.data()) );
	rad->Define(type+"Npim",Form("rad::helpers::Count(%spid,-211)",type.data()) );
	rad->Define(type+"NKp",Form("rad::helpers::Count(%spid,321)",type.data()) );
	rad->Define(type+"NKm",Form("rad::helpers::Count(%spid,-321)",type.data()) );
	rad->Define(type+"Nele",Form("rad::helpers::Count(%spid,11)",type.data()) );
	rad->Define(type+"Npos",Form("rad::helpers::Count(%spid,-11)",type.data()) );
	rad->Define(type+"Npro",Form("rad::helpers::Count(%spid,2212)",type.data()) );
	rad->Define(type+"Nneutron",Form("rad::helpers::Count(%spid,2112)",type.data()) );
    }

  }//reaction
}//rad
