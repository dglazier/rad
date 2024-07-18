#pragma once
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <iostream>

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
     * Reorder vec0 moving entries at iorder0 to iorder1, 
     * with missing elements set to 0
     */
    template<typename T>
    ROOT::VecOps::RVec<T> Reorder(const ROOT::VecOps::RVec<T>& vec0,const ROOT::RVecU& iorder0,const ROOT::RVecU& iorder1,const size_t n){
      //create new vector size  n
      // cout<<"reorder "<< vec0<<iorder0<<iorder1<<" "<<n<<endl;
      ROOT::VecOps::RVec<T> vec1(n); //create new vector size n=iorder1.size
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

    /**
     * count the instances of val in vec
     */
    template<typename T>
    size_t Count(const ROOT::VecOps::RVec<T>& vec,const T& val){
      return std::count(vec.begin(), vec.end(), val);
    }
    /**
     * increment the values in vec
     */
    template<typename T>
    void Increment(ROOT::VecOps::RVec<T>& vec, long off){
      if(vec.size()==0) return ;
      std::for_each(vec.begin(), vec.end(), [&off](T &val) { val+=off; });
     
    }
    /**
     * increment the values in vec and return a copy
     */
    template<typename T>
    ROOT::VecOps::RVec<T> IncrementCopy(ROOT::VecOps::RVec<T> vec, long off){
      if(vec.size()==0) return vec;
      std::for_each(vec.begin(), vec.end(), [&off](T &val) { val+=off; });
      return vec;
    }
    /**
     * Concatenate all given vectors in args pack
     */
    /*
#include <type_traits>
    
    template <typename T>
    void AppendToVec(T t,T& vec) {
      std::copy(t.begin(), t.end(), std::back_inserter(vec));
      //return t;
    }
    

    
    // The base case: we just have a single number.
    template <typename T>
    T Append(T t) {
      return t;
    }
    
    // The recursive case: we take a number, alongside
    // some other numbers, and produce their sum.
    template <typename T, typename... Rest>
    T Append(T t, Rest... rest) {
      //     std::copy(t.begin(), t.end(), std::back_inserter(vec));
      T next = Append(rest...);
      std::copy(next.begin(), next.end(), std::back_inserter(t));
      return t;
    }
    
    template<typename... Args>
    ROOT::RVecD Concatenate(Args &... vs){
      constexpr const auto nArgs = sizeof...(Args);
      const std::size_t sizes[] = {vs.size()...};
      size_t capacity=0;
      for (auto i = 0UL; i < nArgs; i++) {
	capacity+=sizes[i];
      }
      using CT = std::common_type_t<Args...>;
      ROOT::VecOps::RVec<CT> res;
      res.reserve(capacity);
      Append(res,vs...);
      return res;
    }
    */
   template <typename T>
   ROOT::RVec<T> Concatenate(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1)
   {

     ROOT::RVec<T> res;
     res.reserve(v0.size() + v1.size());
     std::copy(v0.begin(), v0.end(), std::back_inserter(res));
     std::copy(v1.begin(), v1.end(), std::back_inserter(res));
      return res;
   }
   template <typename T>
   ROOT::RVec<T> Concatenate(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1,const ROOT::RVec<T> &v2)
   {

     ROOT::RVec<T> res;
     res.reserve(v0.size() + v1.size() + v2.size() );
     std::copy(v0.begin(), v0.end(), std::back_inserter(res));
     std::copy(v1.begin(), v1.end(), std::back_inserter(res));
     std::copy(v2.begin(), v2.end(), std::back_inserter(res));
     return res;
   }
    template <typename T>
   ROOT::RVec<T> Concatenate(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1,const ROOT::RVec<T> &v2, const ROOT::RVec<T> &v3)
   {

     ROOT::RVec<T> res;
     res.reserve(v0.size() + v1.size() + v2.size() + v3.size());
     std::copy(v0.begin(), v0.end(), std::back_inserter(res));
     std::copy(v1.begin(), v1.end(), std::back_inserter(res));
     std::copy(v2.begin(), v2.end(), std::back_inserter(res));
     std::copy(v3.begin(), v3.end(), std::back_inserter(res));
     return res;
   }
    template <typename T>
    ROOT::VecOps::RVec<ROOT::RVec<T>> Concatenate(const ROOT::VecOps::RVec<ROOT::RVec<T>> vec_all)
    {
     const size_t nelements = vec_all.size();
     size_t length=1;
     for(auto i =0UL ; i<nelements;++i){
       length*=vec_all[i].size();
     }
 
     ROOT::RVec<T> res(length);

     size_t entry=0;
     for(auto i =0UL ; i<nelements;++i){
       size_t curr_size=vec_all[i].size();
       for(auto j =0UL ; j<curr_size;++j){
	 res[entry++]=vec_all[i][j];
       }
     }

     return res;
    }
    /** combine v0 and v1 into 1 vector with elements
     * which are vectors of length 2.
     */
    template <typename T>
  ROOT::RVec<T> CombineElements(const ROOT::RVec<T> &v0, const ROOT::RVec<T> &v1)
   {
     const size_t length = v0.size();
     ROOT::VecOps::RVec<ROOT::RVec<T>> res(length, ROOT::RVec<T>(2));
     for(auto i =0UL ; i<length;++i){
       res[i][0]=v0[i];
       res[i][1]=v1[i];
     }
     return res;
   }
    /** combine the n vectors in vec_all into 1 vector with elements
     * which are vectors of length n.
     */
    template <typename T>
    ROOT::VecOps::RVec<ROOT::RVec<T>> CombineElements(const ROOT::VecOps::RVec<ROOT::RVec<T>> vec_all)
   {
     const size_t nelements = vec_all.size();
     const size_t length = vec_all[0].size();
     ROOT::VecOps::RVec<ROOT::RVec<T>> res(length, ROOT::RVec<T>(nelements));
     for(auto i =0UL ; i<length;++i){
       for(auto j =0UL ; j<nelements;++j){
	 res[i][j]=vec_all[j][i];
       }
     }
     return res;
   }
  
 
    /**
     * Return the positions where c is true
     */
    ROOT::RVecI PositionsWhere(const ROOT::RVecB& c)
    {
      const size_t length = c.size();
      ROOT::RVecI r(Count(c,true));//new vector sized to non-zero elements
      uint entry=0;
      for (auto i=0UL; i<length; ++i) {
	if(c[i]==true)
	  r[entry++]=i;
      }
      return r;
    }
  }//helpers
  
}//rad
