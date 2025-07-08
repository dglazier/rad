#pragma once

#include "ParticleCreator.h"
#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "Random.h"
//#include <TRandom.h>

namespace rad{
  namespace generator{
    using ROOT::RVecD;
    
    template<typename Tp, typename Tm>
      RVec<RVecD> CalculateTwoBody( const int &pidx,  const int &bidx, const int &gidx, const RVec<Tm> &masses, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      
      
      //parent meson, recoil baryon, virt photon 4vecs
      auto nbidx=5;
      auto ngidx=2;
      auto meson = FourVector(pidx,px,py,pz,m);
      auto baryon = FourVector(nbidx,px,py,pz,m);
      auto gamma = FourVector(ngidx,px,py,pz,m);
      auto CM = meson+baryon;
    
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
     
      // double costheta = gRandom->Uniform(-1.0,1);
      // double costheta = getUniformRandom(-1,1);
      double costheta = rad::random::Generator().Uniform(-1,1);
      //auto costheta = TMath::Cos(CMParent.Theta());
      //costheta = -0.001;
      auto sintheta = sqrt(1 - costheta*costheta);
      // auto phi = gRandom->Uniform( 0, 2.0*TMath::Pi() );
      //auto phi =  getUniformRandom(0, 2.0*TMath::Pi());
      auto phi =  rad::random::Generator().Uniform(0, 2.0*TMath::Pi());
      
      //define momentum in helicity frame using axes
      auto dpx = p * sintheta * cos(phi);
      auto dpy = p * sintheta * sin(phi);
      auto dpz = p * costheta;
      
      //
      XYZVector V1 = dpx*xV + dpy*yV + dpz*zV;
      //XYZVector V2 = (-dpx)*xV + (-dpy)*yV + (-dpz)*zV;
      
       PxPyPzMVector cmpart1 = {V1.X(),V1.Y(),V1.Z(),masses[0]};
      //PxPyPzMVector cmpart2 = {V2.X(),V2.Y(),V2.Z(),masses[1]};
      
       // PxPyPzMVector  cmpart1 = {dpx,dpy,dpz,masses[0]};
      //PxPyPzMVector cmpart2 = {-dpx,-dpy,-dpz,masses[1]};
      
      //cmpart1 = rot1*cmpart1;
      //cmpart2 = meson-cmpart1;
       //cout << "cm "<<cmpart1 << " " <<TMath::Cos(cmpart1.Theta())<<  " rand gen "<<costheta<<" "<<sintheta<<endl;
      if(TMath::Abs(TMath::Cos(cmpart1.Theta()))>1||TMath::Abs(costheta)>1||TMath::Abs(sintheta)>1) exit(0);
      
      auto part1 = boost(cmpart1,-decBoost);
      //auto part2 = boost(cmpart2,-decBoost);
      auto part2 = meson-part1;
      auto part3=(part1+part2);
      // cout << meson << " " << part3 << " "<<part1<<part2<<endl;
      
      RVec<double> result1={part1.X(),part1.Y(),part1.Z(),masses[0]};
      RVec<double> result2={part2.X(),part2.Y(),part2.Z(),masses[1]};
      
      RVec<RVecD> results;
      results.push_back(result1);
      results.push_back(result2);
      
      return results;
    }
    
    template<typename Tp, typename Tm>
      RVecI ParticleCreateTwoBody(const int &id, const int &pidx, const int &bidx, const int &gidx, const RVec<Tm> masses, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m, const RVecI &iafter){
      
      RVec<RVecD> result = CalculateTwoBody(pidx,bidx,gidx,masses,px,py,pz,m);
      
      int idx = px.size();
	
      px.push_back(result[0][0]);
      py.push_back(result[0][1]);
      pz.push_back(result[0][2]);
      m.push_back(result[0][3]);
      
      px.push_back(result[1][0]);
      py.push_back(result[1][1]);
      pz.push_back(result[1][2]);
      m.push_back(result[1][3]);
      
      return ROOT::RVecI{idx,idx+1};
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
      ParticleGenerator(rad::config::ConfigReaction &cr,size_t seed=1):rad::config::ParticleCreator{cr}{
      Reaction()->InitRandom(seed);
    };
      
      template<typename Tm>
	void GenerateTwoBody(const std::vector<std::string> &names, const RVec<Tm> masses, const string &parent){
        
	std::string smasses=fVecToString(masses);
        TString expr0 = Form( "rad::generator::ParticleCreateTwoBody(0,%s,%i,%i,%s,mc_px,mc_py,mc_pz,mc_m,{})",parent.data(),0,1 , smasses.data());
	auto diparticle = names[0]+"_"+names[1];

	Reaction()->Define(diparticle,expr0.Data());
	Reaction()->Define(names[0],Form("%s[0]",diparticle.data() ));
	Reaction()->Define(names[1],Form("%s[1]",diparticle.data() ));
		       
      }
      
    };// end Class ParticleGenerator

  }//end namespace generator
}//end namespace rad
