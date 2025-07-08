#pragma once

#include "ParticleCreator.h"
#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "Random.h"

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

      //Calculate 2-body breakup rest-frame momentum
      auto sum_masses = masses[0]+masses[1];
      auto dif_masses = masses[0]+masses[1];
      auto term1 = Mpar*Mpar - sum_masses*sum_masses;
      auto term2 = Mpar*Mpar - dif_masses*dif_masses;
      auto p = TMath::Sqrt(term1*term2)/(2*Mpar);
      
      //random thetaphi in helicity frame
      double costheta = rad::random::Generator().Uniform(-1.,1.);
      auto sintheta = TMath::Sqrt(1 - costheta*costheta);
      auto phi =  rad::random::Generator().Uniform(0, 2.0*TMath::Pi());
      
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
      
      int idx = px.size();
      
      px.push_back(result[0]);
      py.push_back(result[1]);
      pz.push_back(result[2]);
      m.push_back(result[3]);
      
      return idx;
    }
    
    
  

    
    class ParticleGenerator : public rad::config::ParticleCreator{
      
    public:
      
      ParticleGenerator() = default;
      ParticleGenerator(rad::config::ConfigReaction &cr,size_t seed=1):
	rad::config::ParticleCreator{cr}
      {
	
      Reaction()->InitRandom(seed);

      };
      

      void GenerateTwoBody(const std::vector<std::string> &names,
			   const ROOT::RVecD& masses, const string &parent){

	//combine masses vector to string
	std::string smasses=rad::utils::combineAnyVectorToString(masses);
	
	//create Define function string
	auto expr = Form( "rad::generator::ParticleCreateTwoBody(%s,%s",parent.data(), smasses.data());
	
	//create first decay particle
	DefineParticle(names[0],std::vector<string>(),expr);
	
	//create second decay particle
	Diff(names[1],{parent},{names[0]});
      }
      
    };// end Class ParticleGenerator

  }//end namespace generator
}//end namespace rad
