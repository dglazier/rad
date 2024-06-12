#pragma once

#include <ROOT/RVec.hxx>
#include <utility>

namespace util{
  template<typename... Args> 
    ROOT::RVecI createRVecI(Args&&... args) {
    ROOT::RVecI vec;
    (vec.emplace_back(std::forward<Args>(args)), ...); // fold-expression
    return vec;
  }
}
