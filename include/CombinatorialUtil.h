#pragma once

#include <ROOT/RVec.hxx>
#include <type_traits> // For std::invoke_result_t
#include <utility>     // For std::forward
//#include "ConfigReaction.h" // For config::RVecIndexMap definition
#include "CommonDefines.h" 

namespace rad {
  namespace util {

    //---------------------------------------------------------
    // Generic Combinatorial Helper Pattern
    //---------------------------------------------------------
 
    /**
     * @brief Vectorizes a single-combination kinematic function across all unique combinations in an event.
     * * This template acts as the central execution engine for combinatorial analysis. It reconstructs 
     * the exact particle index arguments for a given combination from the structured input data, 
     * and calls the user-provided kinematic function once for each combination, returning a vector
     * of results.
     * * @tparam F The type of the single-combination kinematic function (must be explicitly instantiated 
     * with C++ types, e.g., FourVectorMassCalc<double, double>).
     * @tparam Args The types of any additional, non-indexed kinematic arguments (e.g., RVec<double> px, py, etc.).
     * * @param singleComboFunc The function to apply (must accept const RVecIndices& as its first argument).
     * @param combo_indices The structured input holding all candidate indices: 
     * RVec<RVecCombis>[set][particle][combo].
     * @param args Any additional columns/data needed by the kinematic function (e.g., momentum components).
     * * @return ROOT::RVec<T> A vector where T is the return type of singleComboFunc, containing one result per combination.
     * * @note This function ensures that the number of combinations is consistent across all index sets.
     * @see RVecCombis, RVecIndices
     */ 
    template <typename F, typename... Args>
    auto ApplyCombinations(
			   F&& singleComboFunc,
			   const ROOT::RVec<RVecCombis>& combo_indices,
			   Args&&... args) {

      // Determine the return type of the single-combo function
      using T = std::invoke_result_t<F, const RVecIndices&, Args...>;       
 
      //check how many sets of indices this functions requires
      std::vector<UInt_t> n_particles;
      UInt_t n_combos = 0;
      for(const auto& indices:combo_indices){
	n_particles.push_back(indices.size()); //indices.size = Number of particle indices
	if(n_particles.back()==0 ) continue;
       
	auto n = indices[0].size();
	if(n_combos>0 && n!=n_combos){
	  throw std::runtime_error("ApplyCombinations : we have indices with different numbers of combos ! :" + std::to_string(n)+ ","+std::to_string(n_combos));
	}
       
	n_combos = n;  //indices[0].size = Number of combinations,
	//must be same for every element of combo_indices
      }
     
      // Initialize the result vector with the size of the combinations vector
      // and deduced type T
      ROOT::RVec<T> results(n_combos);

      //should check combos1[0].size()==combos2[0].size()
      const auto n_idx_sets=n_particles.size();
     
      // Loop over all combinations 
      for (size_t i_combo = 0; i_combo < n_combos; ++i_combo) {
	//Loop over sets of indices and apply the function

	// Build the consolidated RVec<RVecIndices> argument for F
	RVecIndices indices_for_F(n_idx_sets);
    
	for (size_t i_idx = 0; i_idx < n_idx_sets; ++i_idx) {
	  Indices_t indices(n_particles[i_idx]);
	 
	  for(size_t i_particle=0;i_particle<n_particles[i_idx];++i_particle){
	    indices[i_particle]=combo_indices[i_idx][i_particle][i_combo];
	  }//particle

	  indices_for_F[i_idx]=indices;
	}//indices
	//now have all indices sorted for this combi
	results[i_combo] = std::forward<F>(singleComboFunc)(
							    indices_for_F, 
							    std::forward<Args>(args)...
							    );
	
      }//combi
      
      return results;
    }
    
  } // namespace util
} // namespace rad
