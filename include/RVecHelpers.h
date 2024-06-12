#pragma once
#include <ROOT/RVec.hxx>

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
    ROOT::VecOps::RVec<T> Reorder(const ROOT::VecOps::RVec<T>& vec0,const ROOT::RVecU& iorder0,const ROOT::RVecU& iorder1,const RVecI& size_vec){
      //create new vector same sized as size_vec branch e.g. tru_pid 
      size_t target_size = size_vec.size();
      ROOT::VecOps::RVec<T> vec1(target_size); //create new vector same sized as old
      //Fill new vector by looping over matched indexes and checking if
      //particle exists in recID
       size_t Nentries = iorder1.size(); //just loop over order1 indexes
      for(size_t i=0;i<Nentries;++i){
	vec1[iorder1[i]]=vec0[iorder0[i]]; //give vec1[iorder1] value of vec0[iorder0] 
      }
      return vec1;
    }
    
  }//helpers
  
}//rad
