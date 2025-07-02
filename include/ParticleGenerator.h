#pragma once

#include "ParticleCreator.h"
#include "ConfigReaction.h"
#include "StringUtilities.h"

namespace rad{
  namespace generator{
    using ROOT::RVecD;
    
    template<typename Tp, typename Tm>
      RVec<RVecD> CalculateTwoBody( const int &pidx, const int &bidx, const int &gidx, const RVec<Tm> &masses, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      
      
      //parent meson, recoil baryon, virt photon 4vecs
      auto meson = FourVector(pidx,px,py,pz,m);
      auto baryon = FourVector(bidx,px,py,pz,m);
      auto gamma = FourVector(gidx,px,py,pz,m);
      /* cout << "New event" << endl; */
      /* auto n=px.size(); */
      /* for (int iter=0; iter<n; iter++){ */
      /* 	auto thisvec = FourVector(iter,px,py,pz,m); */
      /* 	cout << iter << ": " << thisvec << endl; */
      /* } */
      /* cout << pidx << " " << meson << endl; */
      /* cout << bidx << " " << baryon << endl; */
      /* cout << gidx << " "<< gamma << endl; */
      /* cout << endl; */
      //get boost
      auto decBoost = meson.BoostToCM();
      //vectors in rest/decay frame of meson
      auto decBar=boost(baryon,decBoost);
      auto decGamma=boost(gamma,decBoost);
      
      XYZVector  zV=decGamma.Vect().Unit();
      XYZVector  yV=(decBar.Vect().Cross(decGamma.Vect())).Unit();
      XYZVector  xV=yV.Cross(zV).Unit();
      
      //pmag of 2body decay
      auto Mpar = meson.M();
      auto term1 = pow(Mpar,2) - pow(masses[0]+masses[1],2);
      auto term2 = pow(Mpar,2) - pow(masses[0]-masses[1],2);
      auto p = sqrt(term1*term2)/(2*Mpar);
      
      /* std::cout << "decBar direction: " << decBar.Vect().Unit() << std::endl; */
      /* std::cout << "decGamma direction: " << decGamma.Vect().Unit() << std::endl; */
      /* std::cout << "zV: " << zV << "  |zV| = " << zV.R() << std::endl; */
      /* std::cout << "yV: " << yV << "  |yV| = " << yV.R() << std::endl; */
      /* std::cout << endl; */
      
      //random thetaphi in helicity frame
      auto costheta = 2.0 * gRandom->Rndm() - 1.0;
      auto sintheta = sqrt(1 - costheta*costheta);
      auto phi = 2.0 * TMath::Pi() * gRandom->Rndm();
      
      //define momentum in helicity frame using axes
      auto dpx = p * sintheta * cos(phi);
      auto dpy = p * sintheta * sin(phi);
      auto dpz = p * costheta;
      
      //
      XYZVector V1 = dpx*xV + dpy*yV + dpz*zV;
      XYZVector V2 = (-dpx)*xV + (-dpy)*yV + (-dpz)*zV;
      
      PxPyPzMVector  cmpart1 = {V1.X(),V1.Y(),V1.Z(),masses[0]};
      PxPyPzMVector cmpart2 = {V2.X(),V2.Y(),V2.Z(),masses[1]};
      
      //PxPyPzMVector  cmpart1 = {dpx,dpy,dpz,masses[0]};
      //PxPyPzMVector cmpart2 = {-dpx,-dpy,-dpz,masses[1]};
      
      auto part1 = boost(cmpart1,-decBoost);
      auto part2 = boost(cmpart2,-decBoost);
      
      RVec<double> result1={part1.X(),part1.Y(),part1.Z(),masses[0]};
      RVec<double> result2={part2.X(),part2.Y(),part2.Z(),masses[1]};
      
      RVec<RVecD> results;
      results.push_back(result1);
      results.push_back(result2);
      
      return results;
    }
    
    template<typename Tp, typename Tm>
      int ParticleCreateTwoBody(const int &id, const int &pidx, const int &bidx, const int &gidx, const RVec<Tm> masses, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m, const RVecI &iafter){
      
      RVec<RVecD> result = CalculateTwoBody(pidx,bidx,gidx,masses,px,py,pz,m);
      
      auto idx = px.size();
	
      px.push_back(result[id][0]);
      py.push_back(result[id][1]);
      pz.push_back(result[id][2]);
      m.push_back(result[id][3]);
      
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
        TString expr0 = Form( "rad::generator::ParticleCreateTwoBody(0,%s,%i,%i,%s",parent.data(), names::BaryonsIdx(), names::VirtGammaIdx(), smasses.data());
	TString expr1 = Form( "rad::generator::ParticleCreateTwoBody(1,%s,%i,%i,%s",parent.data(), names::BaryonsIdx(), names::VirtGammaIdx(), smasses.data() );
	
	DefineParticle(names[0],std::vector<string>(),expr0.Data());
	DefineParticle(names[1],std::vector<string>(),expr1.Data());
      }
      
    };// end Class ParticleGenerator

  }//end namespace generator
}//end namespace rad
