#pragma once

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

void BasicKinematics(){}

//namespace rad{
//RDF interpreter can only take 1 level of namespace...
  namespace rad{
    ///\brief Helper functions and functors for RDF processing
    ///       combining momentum components into 4-vectors
    using ROOT::Math::PxPyPzMVector ;
    using ROOT::Math::XYZVector ;
    using ROOT::RVecF;
    using ROOT::RVecI;
    using ROOT::RVec;
    using ROOT::Math::VectorUtil::boost;

    ///\brief return 4-vector of particle idx
    template<typename Tp, typename Tm>
    PxPyPzMVector FourVector(const uint idx,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      return PxPyPzMVector(px[idx], py[idx], pz[idx], m[idx]);
    }
  
    ///\brief add 4-vectors of particles ip to p4
    template<typename Tp, typename Tm>
    void SumFourVector(PxPyPzMVector& p4, const RVecI &ip,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
    
      auto np = ip.size();
      for (size_t i = 0; i < np; ++i) {
	p4 += PxPyPzMVector(px[ip[i]], py[ip[i]], pz[ip[i]], m[ip[i]]);
      }
    }
  
    ///\brief subtract 4-vectors of particles ip from p4
    template<typename Tp, typename Tm>
    void SubtractFourVector(PxPyPzMVector& p4, const RVecI &ip,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto np = ip.size();
      for (size_t i = 0; i <np ; ++i) {
	p4 -= PxPyPzMVector(px[ip[i]], py[ip[i]], pz[ip[i]], m[ip[i]]);
      }
  
    }
  
    ///\brief return 4-vector of summed particles ipart
    template<typename Tp, typename Tm>
    PxPyPzMVector FourVector(const RVecI &ipart,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
  
      PxPyPzMVector psum(0,0,0,0);
      SumFourVector(psum,ipart,px,py,pz,m);
      return psum;
    
    }
    ///\brief return mass of combined 4-vectors, adding particles with indices ipos and subtracting ineg
    template<typename Tp, typename Tm>
    Tp FourVectorMassCalc(const RVecI &ipos, const RVecI &ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      PxPyPzMVector psum(0,0,0,0);
      SumFourVector(psum,ipos,px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      // std::cout<<"FourVectorMassCalc "<<ipos<<psum<<std::endl
      //       <<std::endl<<pz<<m<<std::endl;
      return psum.M();
    }
    ///\brief return magnitude of momentum
    template<typename T>
    RVec<T> ThreeVectorMag(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      return sqrt(x * x + y * y + z * z);
    }
    ///\brief return eta of momentum
    template<typename T>
    RVec<T> ThreeVectorTheta(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      auto mag = ThreeVectorMag(x,y,z);
      auto costh = z/mag;
      auto test = acos(costh);
      // if(test.size()>4)
      // 	if(test[4]<2.1){
      // 	  cout<<"**********************ThreeVectorTheta"<<test<<x<<y<<z<<endl;
      // 	  exit(0);
      // 	}
      return acos(costh);
    }
    ///\brief return eta of momentum
    template<typename T>
    RVec<T> ThreeVectorPhi(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      return atan2(y,x); //will use vectorised version
    }
   ///\brief return eta of momentum
    template<typename T>
    RVec<T> ThreeVectorEta(const RVec<T> &x, const RVec<T> &y, const RVec<T> &z){
      auto theta = ThreeVectorTheta(x,y,z);
      return -log(tan(0.5 * theta));//will use vectorised version
    }
    ///\brief return x-component
    template<typename T>
    RVec<T> ThreeVectorX(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      return p*sin(theta)*cos(phi);
    }
    ///\brief return x-component
    template<typename T>
    RVec<T> ThreeVectorY(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      return p*sin(theta)*sin(phi);
    }
    ///\brief return x-component
    template<typename T>
    RVec<T> ThreeVectorZ(const RVec<T> &p, const RVec<T> &theta, const RVec<T> &phi){
      return p*cos(theta);
    }

    /* NOTE : we might want to change to using edm4hep functions. Then VecMag would change to 
       template <typename T>
       auto VecMag = [](ROOT::VecOps::RVec<T> momenta) {
       return ROOT::VecOps::Map(momenta, [](const T& p) { return edm4hep::utils::magnitude(p.momentum); });
       };
       //and we would call like Define("pmag", VecMag<edm4hep::MCParticleData> ,{"MCParticles"})
       */


  }//compute
//}//rad
