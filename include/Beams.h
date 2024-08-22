#pragma once
#include <ROOT/RVec.hxx>
#include <array>

namespace rad{
  namespace beams{
    
    using ROOT::Math::PxPyPzMVector ;
    using ROOT::RVec;

    constexpr std::array<double,4>  InitTopComponents();
    constexpr std::array<double,4>  InitBotComponents();
    
    constexpr int InitTopFix(){return 10000;}
    constexpr int InitBotFix(){return 20000;}
    
   
    ///\brief return 4-vector of beam ion
    template<typename Tp, typename Tm>
    PxPyPzMVector InitialFourVector(const uint idx,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //Check if we have been provided with fixed beam parameters
      //or of we have a valid index idx
      switch(idx){
      case InitTopFix() :
	return 	PxPyPzMVector(InitTopComponents()[0],InitTopComponents()[1],
			      InitTopComponents()[2],InitTopComponents()[3]);
      case InitBotFix() :
	return PxPyPzMVector(InitBotComponents()[0],InitBotComponents()[1],
			     InitBotComponents()[2],InitBotComponents()[3]);
	
      default:	
	return PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]);
      }
      
    }
    /*
    ///\brief return 4-vector of beam ion
    template<typename Tp, typename Tm>
    PxPyPzMVector BeamIonFourVector(const uint idx,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //If beam ion 4-vector fixed use fixed values
      //if not use idx given
      auto check = (idx!= BeamIonFix());
      //if(check==1&&idx>px.size()) //should probably assert here...
      return check ?
	PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]) :
	PxPyPzMVector(BeamIonComponents()[0],BeamIonComponents()[1],
		      BeamIonComponents()[2],BeamIonComponents()[3]);
    }
   
   
    ///\brief return 4-vector of beam electron
    template<typename Tp, typename Tm>
    PxPyPzMVector BeamEleFourVector(const uint idx,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //If beam ion 4-vector fixed use fixed values
      //if not use idx given
      auto check = (idx!= BeamEleFix());
      return check?
	PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]) :
	PxPyPzMVector(BeamEleComponents()[0],BeamEleComponents()[1],
		      BeamEleComponents()[2],BeamEleComponents()[3]);
    }
    ///\brief return 4-vector of beam ion
    template<typename Tp, typename Tm>
    PxPyPzMVector TargetIonFourVector(const uint idx,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //If beam ion 4-vector fixed use fixed values
      //if not use idx given
      auto check = (idx!= TargetIonFix());
      //if(check==1&&idx>px.size()) //should probably assert here...
      return check ?
	PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]) :
	PxPyPzMVector(TargetIonComponents()[0],TargetIonComponents()[1],
		      TargetIonComponents()[2],TargetIonComponents()[3]);
    }
   
   
    ///\brief return 4-vector of beam electron
    template<typename Tp, typename Tm>
    PxPyPzMVector BeamGamFourVector(const uint idx,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //If beam ion 4-vector fixed use fixed values
      //if not use idx given
      auto check = (idx!= BeamGamFix());
      return check?
	PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]) :
	PxPyPzMVector(BeamGammaComponents()[0],BeamGammaComponents()[1],
		      BeamGammaComponents()[2],BeamGammaComponents()[3]);
    }
    */
  }//beams
}//rad
