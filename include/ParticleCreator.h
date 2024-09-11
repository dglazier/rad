#pragma once

#include "ConfigReaction.h"

namespace rad{
  namespace config{

    class ParticleCreator {

    public:
      
      ParticleCreator() = default;
    ParticleCreator(ConfigReaction& cr):_reaction{&cr}{};
      
      
      //////////////////////////////////////////////////////////////////
      void Sum(const string& name,const std::vector<std::string>& parts){
	//adding ( at end of string allows extra arguments to be added, like for Miss
	DefineParticle(name,parts,"rad::reactkine::ParticleCreateBySum(");
      }
       //////////////////////////////////////////////////////////////////
      void Miss(const string& name,const std::vector<std::string>& parts){
	//THIS DOES NOT WORK YET!!!!!
	DefineParticle(name,parts,Form("rad::reactkine::ParticleCreateByMiss(%s,",BeamIndices().data()));
      }
      //////////////////////////////////////////////////////////////////
      void DefineParticle(const string& name,const std::vector<std::string> parts,const string& funcExpr){
	
	//store all defined particle names
	std::vector<std::string> names;

	//loop over ConfigReaction types and define this particle for each
	auto types = _reaction->GetTypes();
	for(auto &atype:types){
	  //Make sure any created particle uses type_idx
	  //This ensures its Create function is called prior to this one
	  string sum ="{";
	  for(string p:parts){
	    if( std::find(_created.begin(),_created.end(),p)!=_created.end() ){
	      std::string type_p = atype.first + p;
	      p=type_p;
	    }
	    sum=(sum+p+",");
	  }
	  sum.pop_back(); //remove last ,
	  sum+='}';

	  //Note we give all created particles as argument to ensure creation order
	  auto type_created=_created;
	  for(auto& col: type_created) col=atype.first+col;
	  auto after_cols  = VectorToString(type_created);
	  
	  //format args "func(idxs,name,components,after_idxs")
 	  TString type_expr = Form("%s%s,%s,%s)",funcExpr.data(),sum.data(),atype.second["components_p4"].data(),after_cols.data());
	  std::cout<<"DefineParticle "<<type_expr<<std::endl;
	  names.push_back(atype.first + name.data());
	  _reaction->Define(atype.first + name.data(),type_expr.Data());
	}
	_created.push_back(name);
	
	auto snames = VectorToString(names);
	//define name as the first type entry in names
	//this function ensures all type create particles are called at same time
	// std::cout<<"Sum names "<<snames<<" "<<Form("ROOT::RVecU%s[0]",snames.data())<<std::endl;
	_reaction->Define(name.data(),Form("ROOT::RVecI%s[0]",snames.data()));
      }
    
      //////////////////////////////////////////////////////////////////
      std::string VectorToString(const std::vector<std::string>& parts){
	if(parts.empty()==true) return "{}";
      
	string toString ="{";
	for(const auto& p:parts){
	  toString=(toString+p+",");
	}
	toString.pop_back(); //remove last ,
	toString+='}';
	return toString;
      }
     
    private:
      
      ConfigReaction* _reaction=nullptr;
      std::vector<string> _created;
      
    };
    
  }
}
