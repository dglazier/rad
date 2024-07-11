#pragma once
#include <ROOT/RVec.hxx>
#include <algorithm>

namespace rad{
  namespace helpers{
    
     //! Code simplifications
  
    using ROOT::RVecI;
    
    /**
     * Find the next occurance of value in vec after vecstart
     */
    inline size_t findIndex(const ROOT::RVecI& vec,const int value,const int* vecstart=nullptr) {
      if(vecstart==nullptr) vecstart = vec.begin();
      return distance(vec.begin(), find(vecstart, vec.end(), value ));
    }

     /**
     * Find the nth occurance of value in vec. Note n_occurance should be 1 for first occurance
     */
    int findNthIndex(const ROOT::RVecI& vec, int n_occurance,const int value){
      if(n_occurance == 0 ) return -1;
      //find first occurance
      auto currIdx = findIndex(vec,value);
      auto vsize= vec.size();
      //check if we reached the end of the vector
      if(currIdx==vsize) return -1;

      //Look for the n_occurance case
      n_occurance--;
      while( (n_occurance--) >0){
	currIdx=findIndex(vec,value,&vec[currIdx+1]);
	//check if we reached the end of the vector
	if(currIdx==vsize) return -1;
      }
      return currIdx;
    }
    /**
     * Rearrange vec in order of elements in imatch
     */
    template<typename T>
    ROOT::VecOps::RVec<T> Rearrange(const ROOT::VecOps::RVec<T>& vec,const ROOT::RVecU& imatch){
       return ROOT::VecOps::Take(vec,imatch);
    }
    /**
     *reorder vec to elements of matchid, with missing elements set to 0
     */
    template<typename T>
    ROOT::VecOps::RVec<T> Reorder(const ROOT::VecOps::RVec<T>& vec0,const ROOT::RVecU& iorder0,const ROOT::RVecU& iorder1,const size_t n){
      //create new vector size  n
      ROOT::VecOps::RVec<T> vec1(n); //create new vector same sized as old
      //need to loop over order0
      size_t target_size = iorder1.size();
      for(size_t i=0;i<target_size;++i){
	//add value of vec0 at iorder1
	vec1[iorder1[i]]=vec0[iorder0[i]]; //give vec1[iorder1] value of vec0[iorder0]
      }
      return vec1;
    }

    template <typename T>
     ROOT::RVec<T> Enumerate(size_t size)
    //  ROOT::RVec<typename  ROOT::RVec<T>::size_type> Enumerate(size_t size)
    {
      ROOT::RVec<T> ret;
      ret.reserve(size);
      for (auto i = 0UL; i < size; ++i) {
	ret.emplace_back(i);
      }
      return ret;
    }
   /**
     * Truncate vec at n
     */
    template<typename T>
    ROOT::VecOps::RVec<T> Truncate(const ROOT::VecOps::RVec<T>& vec,const size_t size){
      ROOT::RVec<T> ret;
      ret.reserve(size);
      for (auto i = 0UL; i < size; ++i) {
	ret.emplace_back(vec[i]);
      }
      return ret;
    }

    template<typename T>
   size_t Count(const ROOT::VecOps::RVec<T>& vec,const T& val){
     return std::count(vec.begin(), vec.end(), val);
   }
    
  }//helpers
  
}//rad
