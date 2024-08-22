#pragma once
#include "Beams.h"
#include "BasicKinematics.h"
#include "ConfigReaction.h"

void ReactionKinematics(){}

namespace rad{
  namespace reactkine{
  
    ///\brief create a new particle and add it to the momentum vectors
    ///return the index to be used to access the components
    ///this also allows it to be used with Define which requires a return
    template<typename Tp, typename Tm>
    int ParticleCreateBySum(const RVecI& isum, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m,const RVecI& iafter){
      //sum the 4-vectors
      // std::cout<<"ParticleCreateBySum "<<isum<<m<<" "<<m.size()<<std::endl;
      PxPyPzMVector p4;
      SumFourVector(p4,isum,px,py,pz,m);
      //make particle id = last entry
      auto idx = px.size();
      //add new components
      px.push_back(p4.X());
      py.push_back(p4.Y());
      pz.push_back(p4.Z());
      m.push_back(p4.M());
      return idx;
    }
    ///\brief create a new particle and add it to the momentum vectors
    ///return the index to be used to access the components
    ///this also allows it to be used with Define which requires a return
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    int  ParticleCreateByMiss(const config::RVecIndexMap react,const RVecI& ineg, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m,const RVecI& iafter){
      //sum the 4-vectors
      auto p4 = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      p4+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
        
      SubtractFourVector(p4,ineg,px,py,pz,m);

      //make particle id = last entry
      auto idx = px.size();

      //add new components
      px.push_back(p4.X());
      py.push_back(p4.Y());
      pz.push_back(p4.Z());
      m.push_back(p4.M());
      return idx;
    }
 
    
   ///\brief missing mass fo reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    Tp FourVectorMissMassCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.M();
    }

   ///\brief missing mass sqaured of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    Tp FourVectorMissMass2Calc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.M2();
    }

   ///\brief functions to compute standard reaction kinematics
    template<typename Tp, typename Tm>
    PxPyPzMVector CMVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      RVecI icm=react[names::BaryonsIdx()];
      icm.insert(icm.end(),react[names::MesonsIdx()].begin(),react[names::MesonsIdx()].end());
    
      return FourVector(icm,px,py,pz,m);
    }

    /**
     * reaction baryon 4-vector
     */
    template<typename Tp, typename Tm>
    PxPyPzMVector BaryonFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      return FourVector(react[names::BaryonsIdx()],px,py,pz,m);
    }
    /**
     * reaction meson 4-vector
     */
    template<typename Tp, typename Tm>
    PxPyPzMVector MesonFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      return FourVector(react[names::MesonsIdx()],px,py,pz,m);
    }

   
    template<typename Tp, typename Tm>
    Tp T0(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      //   auto tar = beams::BeamIonFourVector(react[names::BeamIonIdx()][0],px,py,pz,m);
      auto tar = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      auto bar = FourVector(react[names::BaryonsIdx()],px,py,pz,m);
      //generate CM from sum of final state meson and baryon particles
      auto cm = CMVector(react,px,py,pz,m);
      auto cmBoost = cm.BoostToCM();
      PxPyPzMVector  CMTar=boost(tar,cmBoost);
      PxPyPzMVector  CMBar=boost(bar,cmBoost);
    
      //return  M1*M1 + M3*M3  - 2 * ( E1*E3 -p1*p3*costh );
      Double_t t0 = CMBar.M2() + CMTar.M2() - 2*(CMBar.E()*CMTar.E() - CMBar.P()*CMTar.P());
    
      return t0;
    }
  
   
    ///\brief return 4 momentum transfer squared of "in particles" - "out particles" on bottom vertex
    template<typename Tp, typename Tm>
    Tp TBot(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //std::cout<<" TBot "<<react[names::InitialBotIdx()][0]<<" "<<react[names::BaryonsIdx()]<<" "<<pz<<m<<std::endl;
      auto psum = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      //std::cout<<psum<<std::endl;
      SubtractFourVector(psum,react[names::BaryonsIdx()],px,py,pz,m);
      //std::cout<<psum<<" "<<- (psum.M2())<<std::endl;
      return - (psum.M2());
    }
    ///\brief return 4 momentum transfer squared of "in particles" - "out particles" on top vertex
    template<typename Tp, typename Tm>
    Tp TTop(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
    
      //std::cout<<" TTop "<<react[names::InitialTopIdx()][0]<<" "<<react[names::MesonsIdx()]<<" "<<pz<<m<<std::endl;
      auto psum = PhotoFourVector(react,px,py,pz,m);//beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      //std::cout<<psum<<" "<< FourVector(react[names::MesonsIdx()],px,py,pz,m)<<std::endl;
      
      SubtractFourVector(psum,react[names::MesonsIdx()],px,py,pz,m);
      //      std::cout<<psum<<" "<<- (psum.M2())<<std::endl;
      return - (psum.M2());
    }

  
    ///\brief return 4 momentum transfer squared, t, minus t0 (or tmin) of "in particles" - "out particles"
    template<typename Tp, typename Tm>
    Tp TPrimeBot(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
  
      return TBot(react,px,py,pz,m) - T0(react,px,py,pz,m);
    }
    ///\brief return 4 momentum transfer squared, t, minus t0 (or tmin) of "in particles" - "out particles"
    template<typename Tp, typename Tm>
    Tp TPrimeTop(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
  
      return TTop(react,px,py,pz,m) - T0(react,px,py,pz,m);
    }
  


  }

}
