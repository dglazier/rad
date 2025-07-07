#pragma once

#include "ParticleCreator.h"
#include "ConfigReaction.h"
#include "StringUtilities.h"

namespace rad{
  namespace generator{
    using ROOT::RVecD;
    
    template<typename Tp, typename Tm>
      RVec<double> CalculateTwoBody( const int &pidx,  const int &bidx, const int &gidx, const RVec<Tm> &masses, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      
      //cout << "idx meson, baryon, gamma: " << pidx << " " << bidx << " " << gidx << endl;
      
      //parent meson, recoil baryon, virt photon 4vecs
      auto meson = FourVector(pidx,px,py,pz,m);
      auto baryon = FourVector(bidx,px,py,pz,m);
      auto gamma = FourVector(gidx,px,py,pz,m);
      
      //get boost
      auto decBoost = meson.BoostToCM();
      //vectors in rest/decay frame of meson
      auto decBar=boost(baryon,decBoost);
      auto decGamma=boost(gamma,decBoost);
      
      XYZVector  zV=-decBar.Vect().Unit();
      XYZVector  yV=(decBar.Vect().Cross(decGamma.Vect())).Unit();
      XYZVector  xV=yV.Cross(zV).Unit();
      
      //pmag of 2body decay
      auto Mpar = meson.M();
      auto term1 = pow(Mpar,2) - pow(masses[0]+masses[1],2);
      auto term2 = pow(Mpar,2) - pow(masses[0]-masses[1],2);
      auto p = sqrt(term1*term2)/(2*Mpar);
      
      //random thetaphi in helicity frame
      auto costheta = gRandom->Uniform(-1,1);
      auto sintheta = sqrt(1 - costheta*costheta);
      auto phi = gRandom->Uniform( 0, 2.0*TMath::Pi() );
      
      //define momentum in helicity frame using axes
      auto dpx = p * sintheta * cos(phi);
      auto dpy = p * sintheta * sin(phi);
      auto dpz = p * costheta;
      
      //
      XYZVector V1 = dpx*xV + dpy*yV + dpz*zV;
      XYZVector V2 = (-dpx)*xV + (-dpy)*yV + (-dpz)*zV;
      
      PxPyPzMVector cmpart1 = {V1.X(),V1.Y(),V1.Z(),masses[0]};
      PxPyPzMVector cmpart2 = {V2.X(),V2.Y(),V2.Z(),masses[1]};
      
      auto part1 = boost(cmpart1,-decBoost);
      RVec<double> result={part1.X(),part1.Y(),part1.Z(),part1.M()};
      return result;
    }
    
    template<typename Tp, typename Tm>
      int ParticleCreateTwoBody(const int &pidx, const int &bidx, const int &gidx, const RVec<Tm> masses, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m, const RVecI &iafter){
      
      auto result = CalculateTwoBody(pidx,bidx,gidx,masses,px,py,pz,m);
      
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
	void GenerateTwoBody(const std::vector<std::string> &names, const RVec<Tm> masses, const string &parent, const string &baryon){
        
	std::string smasses=fVecToString(masses);
        
	Diff("gen_virt_gam",{rad::names::BeamEle()},{rad::names::ScatEle()});
	TString expr0 = Form( "rad::generator::ParticleCreateTwoBody(%s,%s,%s,%s",parent.data(), baryon.data(), "gen_virt_gam", smasses.data());
	DefineParticle(names[0],std::vector<string>(),expr0.Data());
	Diff(names[1],std::vector<string>{parent.data()},std::vector<string>{names[0]});
	
      }
      
    };// end Class ParticleGenerator

  }//end namespace generator
}//end namespace rad
