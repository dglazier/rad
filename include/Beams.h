#pragma once
#include <ROOT/RVec.hxx>
#include <array>

namespace rad{
  namespace beams{
    
   using ROOT::Math::PxPyPzMVector ;
   using ROOT::RVec;

   constexpr std::array<double,4>  BeamIonComponents();// {return {0.,0.,100.,0.93827210};}
   constexpr std::array<double,4>  BeamEleComponents();// {return {0.,0.,-10.,0.00051099900};}
   constexpr int BeamIonFix(){return 10000;}
   constexpr int BeamEleFix(){return 20000;}
   
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
   
  }//beams
}//rad
