#pragma once
#include "Beams.h"
#include "BasicKinematics.h"
#include "ConfigReaction.h"

void ReactionKinematics(){}

namespace rad{
  namespace reactkine{
  
   
   ///\brief missing mass fo reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    Tp MissMass(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      //auto psum = beams::InitialFourVector(react[names::ElectroIonIdx()][0],px,py,pz,m);
      //psum+=beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.M();
    }

   ///\brief missing mass sqaured of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    Tp MissMass2(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      //auto psum = beams::InitialFourVector(react[names::ElectroIonIdx()][0],px,py,pz,m);
      //psum+=beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      
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
    
      auto psum = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);;   
      SubtractFourVector(psum,react[names::BaryonsIdx()],px,py,pz,m);
      return - (psum.M2());
    }
    ///\brief return 4 momentum transfer squared of "in particles" - "out particles" on top vertex
    template<typename Tp, typename Tm>
    Tp TTop(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
    
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,react[names::MesonsIdx()],px,py,pz,m);
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
