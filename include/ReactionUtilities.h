#pragma once

//#include "ConfigReaction.h"
#include <TString.h>
#include <ROOT/RDFHelpers.hxx>

namespace rad{
  namespace reaction{
    namespace util{

    
      // void CountParticles(rad::config::ConfigReaction* rad, const std::string& type){
      template<typename T> //use template so can #include this in ConfigReaction
      void CountParticles(T* rad, const std::string& type){
      
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
    
      //////////////////////////////////////////////////////////////////
      std::string ColumnsToString(const ROOT::RDF::ColumnNames_t &cols) {
	if(cols.empty()==true) return "{}";
      
	string toString ="{";
	for(const auto& p:cols){
	  toString=(toString + p + ",");
	}
	toString.pop_back(); //remove last ,
	toString+='}';
	return toString;
      }
    }
  }//reaction
}//rad
