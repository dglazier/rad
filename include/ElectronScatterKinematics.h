#pragma once
#include "ReactionKinematics.h"


namespace rad{
  namespace electro{
  

    template<typename Tp, typename Tm>
    PxPyPzMVector PhotonVector(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto phot =  beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      SubtractFourVector(phot,react[names::ScatEleIdx()],px,py,pz,m);
      return phot;
    }

    template<typename Tp, typename Tm>
    Tp Q2(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto phot = PhotoFourVector(react,px,py,pz,m);
      return -phot.M2();
    
    }

    //would like to use a struct here, but define cannot take structs or classes
    //must be basic types or RVecs of basic types
    //HAve added method to deal with structs by defining seperate columns for
    //each member. However this seems slower than doing the calc twice!

    struct ElCMDecay_t{
      double CosTheta=0.;
      double Phi=0.;
    };

 
    template<typename Tp, typename Tm>
    XYZVector ElectroCMDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //  ElCMDecay_t ElectroCMDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //CM frame defined by e-scattering
      auto cm = reactkine::CMVector(react,px,py,pz,m);
      auto cmBoost = cm.BoostToCM();
      auto beam = beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      auto mes = FourVector(react[names::MesonsIdx()],px,py,pz,m);
      auto photon = PhotoFourVector(react,px,py,pz,m);
  
      PxPyPzMVector CMBeam=boost(beam,cmBoost);
      PxPyPzMVector CMMes=boost(mes,cmBoost);
      PxPyPzMVector CMGamma=boost(photon,cmBoost);
  
      XYZVector zV=CMGamma.Vect().Unit();
      XYZVector yV=CMGamma.Vect().Cross(CMBeam.Vect()).Unit();
      XYZVector xV=yV.Cross(zV).Unit();
  
      XYZVector angles(CMMes.Vect().Dot(xV),CMMes.Vect().Dot(yV),CMMes.Vect().Dot(zV));
      // ElCMDecay_t result;
      // result.CosTheta=(TMath::Cos(angles.Theta()));
      // result.Phi=angles.Phi();
      // return result;
      return angles;
    }
  
    template<typename Tp, typename Tm>
    Tp CosThetaCM(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = ElectroCMDecay(react,px,py,pz,m);
      return TMath::Cos(angles.Theta());
    }
  
    template<typename Tp, typename Tm>
    Tp PhiCM(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = ElectroCMDecay(react,px,py,pz,m);
      return angles.Phi();
    }


  

  }
}
