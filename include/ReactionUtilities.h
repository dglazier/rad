#pragma once

//#include "ConfigReaction.h"
#include "DefineNames.h"

#include <TString.h>
#include <ROOT/RDFHelpers.hxx>

namespace rad{
  namespace reaction{
    namespace util{
      
      
    
      // void CountParticles(rad::config::ConfigReaction* rad, const std::string& type){
      template<typename T> //use template so can #include this in ConfigReaction
      void CountParticles(T* rdf, const std::string& type){
      
	rdf->Define(type+"Ngamma",Form("rad::helpers::Count(%spid,22)",type.data()) );
	rdf->Define(type+"Npip",Form("rad::helpers::Count(%spid,211)",type.data()) );
	rdf->Define(type+"Npim",Form("rad::helpers::Count(%spid,-211)",type.data()) );
	rdf->Define(type+"NKp",Form("rad::helpers::Count(%spid,321)",type.data()) );
	rdf->Define(type+"NKm",Form("rad::helpers::Count(%spid,-321)",type.data()) );
	rdf->Define(type+"Nele",Form("rad::helpers::Count(%spid,11)",type.data()) );
	rdf->Define(type+"Npos",Form("rad::helpers::Count(%spid,-11)",type.data()) );
	rdf->Define(type+"Npro",Form("rad::helpers::Count(%spid,2212)",type.data()) );
	rdf->Define(type+"Nneutron",Form("rad::helpers::Count(%spid,2112)",type.data()) );
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
    

      using rad::names::data_type::Rec;
      using rad::names::data_type::Truth;

      /**
       * calculate the difference in reconsutructed and truth variables
       * Case Reconstructed and truth synched via AliasColumnsAndMatchWithMC()
       */
      template<typename T> //use template so can #include this in ConfigReaction
      void Resolution(T* const  rdf,const string& var){
	// Define(string("res_")+var,[](const ROOT::RVec<T> &rec,const ROOT::RVec<T> &tru){
	//   return (rec - tru);
	// },{string(Rec())+var,string(Truth())+var});
	rdf->Define(string("res_")+var,Form("%s-%s",(Truth()+var).data(),(Rec()+var).data() ));
      }
      template<typename T> //use template so can #include this in ConfigReaction
      void ResolutionFraction(T* const rdf,const string& var){
	// Define(string("res_")+var,[](const ROOT::RVec<T> &rec,const ROOT::RVec<T> &tru){
	//    return (rec - tru)/tru;
	// },{Rec()+var,Truth()+var});
	rdf->Define(string("res_")+var,Form("(%s-%s)/%s",(Truth()+var).data(),(Rec()+var).data(),(Truth()+var).data() ));
      }
    }
  }//reaction
}//rad
