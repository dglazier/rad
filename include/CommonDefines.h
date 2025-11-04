#pragma once

namespace rad{

  //kinematics functions etc should return common type
  using ResultType_t = Double_t;
  using RVecResultType = ROOT::RVec<ResultType_t>;
    
  using Indices_t = ROOT::RVecI;
  using RVecCombis = ROOT::RVec<Indices_t>; //indices : RVecCombis[particle][combo]
  using RVecIndices = ROOT::RVec<Indices_t>; //indices : RVecIndices[set][particle]
  using RVecIndexMap = ROOT::RVec<Indices_t>;

}
