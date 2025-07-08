#pragma once

#include "ParticleCreator.h"
#include "ConfigReaction.h"
#include "StringUtilities.h"

namespace rad{
  namespace generator{
    using ROOT::RVecD;
    
    template<typename Tp, typename Tm>
      RVec<double> CalculateTwoBody( const int &pidx, const RVec<Tm> &masses, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
           
      //get boost
      auto parent = FourVector(pidx,px,py,pz,m);
      auto decBoost = parent.BoostToCM();
      
      //pmag of 2body decay
      auto Mpar = parent.M();
      auto term1 = pow(Mpar,2) - pow(masses[0]+masses[1],2);
      auto term2 = pow(Mpar,2) - pow(masses[0]-masses[1],2);
      auto p = sqrt(term1*term2)/(2*Mpar);
      
      //random thetaphi in parent CM frame
      auto costheta = gRandom->Uniform(-1,1);
      auto sintheta = sqrt(1 - costheta*costheta);
      auto phi = gRandom->Uniform( 0, 2.0*TMath::Pi() );
      
      //define momentum components
      auto dpx = p * sintheta * cos(phi);
      auto dpy = p * sintheta * sin(phi);
      auto dpz = p * costheta;
      
      //define first decay product
      PxPyPzMVector cmpart1 = {dpx,dpy,dpz,masses[0]};
      //boost back to lab frame
      auto part1 = boost(cmpart1,-decBoost);
      
      RVec<double> result={part1.X(),part1.Y(),part1.Z(),part1.M()};
      return result;
    }
    
    template<typename Tp, typename Tm>
      int ParticleCreateTwoBody(const int &pidx, const RVec<Tm> masses, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m, const RVecI &iafter){
      
      auto result = CalculateTwoBody(pidx,masses,px,py,pz,m);
      
      auto idx = px.size();
	
      px.push_back(result[0]);
      py.push_back(result[1]);
      pz.push_back(result[2]);
      m.push_back(result[3]);
      
      return idx;
    }
    

    template<typename T>
      std::string fVecToString(const RVec<T> vec){
      std::string result = "{";
      
      //go through vec and make string of each element
      for (auto iter:vec){
	auto p = std::to_string(iter);
	result.append(p);
	result.append(",");
      }

      //remove last "," easily
      if(!result.empty()){
	result.pop_back();
      }
      //close the curly brackets {}
      result.append("}");
      return result;
    }
    class ParticleGenerator : public rad::config::ParticleCreator{
      
    public:
      
      ParticleGenerator() = default;
    ParticleGenerator(rad::config::ConfigReaction &cr):rad::config::ParticleCreator{cr}{};
      
      template<typename Tm>
	void GenerateTwoBody(const std::vector<std::string> &names, const RVec<Tm> masses, const string &parent){
        
	std::string smasses=fVecToString(masses);
        
	TString expr0 = Form( "rad::generator::ParticleCreateTwoBody(%s,%s",parent.data(), smasses.data());
	DefineParticle(names[0],std::vector<string>(),expr0.Data());
	Diff(names[1],std::vector<string>{parent.data()},std::vector<string>{names[0]});
	
      }
      
    };// end Class ParticleGenerator

  }//end namespace generator
}//end namespace rad
