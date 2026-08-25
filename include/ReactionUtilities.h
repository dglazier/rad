#pragma once

#include "DefineNames.h"

#include <TString.h>
#include <ROOT/RDFHelpers.hxx>

namespace rad{
    namespace util{

      using rad::consts::data_type::Rec;
      using rad::consts::data_type::Truth;

      /**
       * @brief Defines columns counting the occurrences of standard PIDs (gamma, pi, K, p, n, e).
       * * Creates columns like `rec_Ngamma`, `rec_Npip`, etc.
       * @tparam T The RDataFrame type.
       * @param rdf Pointer to the RDataFrame interface.
       * @param type Data prefix (e.g. "rec_").
       */
      template<typename T>
      void CountParticles(T* rdf, const std::string& type);

      /**
       * @brief Calculates simple resolution: (Truth - Rec).
       * * Requires inputs to be aliased/matched first.
       * Creates column `res_{var}`.
       */
      template<typename T> 
      void Resolution(T* const rdf, const std::string& var);
      
      /**
       * @brief Calculates fractional resolution: (Truth - Rec) / Truth.
       * Creates column `res_{var}`.
       */
      template<typename T> 
      void ResolutionFraction(T* const rdf, const std::string& var);

      /**
     * @brief Returns mass in GeV for detector-stable particles.
     * @details Evaluated entirely at compile-time. Uses absolute value 
     * to automatically handle anti-particles perfectly.
     */
      constexpr double PdgToMass(int pdg);

      /**
       * @brief returns synched vector of masses given PDG codes
       */      
      ROOT::RVecD AssignMasses( const ROOT::RVecI &pid);
      
    } // namespace util
} // namespace rad

// =================================================================================
// IMPLEMENTATION
// =================================================================================

namespace rad {
  namespace util {
    
    /**
     * @brief Debugging tool to count how many events satisfy a specific PID topology.
     * * Useful for diagnosing "0 events" issues by checking if the raw input arrays 
     * actually contain the required particles before combinatorics runs.
     * * * **Usage:**
     * @code
     * // Check for 3 Electrons (11), 1 Positron (-11), 1 Pion (211)
     * rad::util::CountTopo(&rad_df, Rec(), {{3, 11}, {1, -11}, {1, 211}});
     * @endcode
     * * @tparam T The Reaction class type (e.g. ePICReaction).
     * @param rdf Pointer to the reaction interface.
     * @param type Data prefix (e.g. "rec_" or "tru_").
     * @param topo Vector of pairs: `{Required_Count, PID}`.
     */
    template<typename T>
    inline void CountTopo(T* rdf, const std::string& type, const ROOT::RVec<std::pair<unsigned int, int>>& topo) {
      
      std::string col_name = type + "has_topo_debug";

      // Define boolean pass/fail based on PID counts
      rdf->Define(col_name, [topo](const ROOT::RVecI& pid) {
        bool pass = true;
        for(const auto& p : topo) {
          unsigned int required_count = p.first;
          int target_pid = p.second;
          
          // Count instances of this PID
          int actual_count = ROOT::VecOps::Sum(pid == target_pid);
          
          // Logic: We need AT LEAST the required amount. 
          // Strict equality (==) is too restrictive (e.g. 3 electrons is valid for a 2-electron topology).
          if(actual_count < required_count) {
            pass = false;
            break;
          }
        }
        return pass;
      }, {type + "pid"});
      
      // Execute Count immediately
      auto count_result = rdf->CurrFrame().Filter(col_name).Count();
      
      std::cout << "------------------------------------------------------------\n"
                << " [CountTopo] Checking Topology " << type << ":\n";
      for(const auto& p : topo) {
	std::cout << "   - PID " << p.second << ": Require >= " << p.first << "\n";
      }
      std::cout << " => Events Passing: " << *count_result << " / " << *rdf->CurrFrame().Count() << "\n"
                << "------------------------------------------------------------" << std::endl;
    }

    /**
     * @brief Prints yield statistics for each particle independently.
     * * Useful for debugging: "Do I have ANY positrons?"
     * * **Lazy Execution:** Queues definitions and runs the event loop only ONCE at the end.
     */
    template<typename T>
    inline void CheckCandidateStats(T* reaction, const std::string& prefix) {
      
      auto names = reaction->ParticleNames();
      if(names.empty()) return;

      // Local struct to hold the "Promises" (Smart Pointers) cleanly
      struct StatHolder {
          std::string name;
          ROOT::RDF::RResultPtr<double> mean;
          ROOT::RDF::RResultPtr<double> gt0;
          ROOT::RDF::RResultPtr<double> gt1;
      };
      
      ROOT::RVec<StatHolder> stats;
      stats.reserve(names.size());

      auto rdf = reaction->CurrFrame();

      // 1. Lazy Booking Loop
      for(const auto& name : names) {
         std::string col_name = prefix + name;
         if(!reaction->ColumnExists(col_name)) continue;

         std::string n_col = "n_" + prefix + name; 
         
         // Define count column (cast to double for Sum/Mean compatibility)
         rdf = rdf.Define(n_col, [](const ROOT::RVecI& inds){ return (double)inds.size(); }, {col_name});

         // Queue the actions (Does not run loop yet)
         stats.push_back({
             name,
             rdf.Mean(n_col),                                         // Mean count per event
             rdf.Define(n_col+"_gt0", n_col+" > 0").Sum(n_col+"_gt0"), // Events with >0 candidates
             rdf.Define(n_col+"_gt1", n_col+" > 1").Sum(n_col+"_gt1")  // Events with >1 candidates
         });
      }
      
      // Queue Total Event Count
      auto total_ptr = rdf.Count(); 

      // 2. Trigger Execution & Print
      // Accessing *total_ptr triggers the RDataFrame loop once for all calculations.
      auto total_events = *total_ptr;

      std::cout << "\n==========================================================\n";
      std::cout << " [RAD] Candidate Yield Report (" << prefix << ")\n";
      std::cout << "==========================================================\n";
      std::cout << std::left << std::setw(15) << "Particle" 
                << std::setw(15) << "Events (>0)" 
                << std::setw(15) << "Events (>1)" 
                << "Mean/Event" << std::endl;
      std::cout << "----------------------------------------------------------\n";

      for(auto& s : stats) {
          // Dereference pointers safely to get values
	double mean_val = *(s.mean);
	double gt0_val  = *(s.gt0);
	double gt1_val  = *(s.gt1);

          std::cout << std::left << std::setw(15) << s.name 
                    << std::setw(15) << (long long)gt0_val 
                    << std::setw(15) << (long long)gt1_val 
                    << std::setprecision(2) << mean_val << std::endl;
      }

      std::cout << "----------------------------------------------------------\n";
      std::cout << " Total Events Processed: " << total_events << "\n";
      std::cout << "==========================================================\n" << std::endl;
    }
      template<typename T>
      inline void CountParticles(T* rdf, const std::string& type){
        rdf->Define(type+"Ngamma",   Form("rad::util::Count(%spid,22)",   type.data()) );
        rdf->Define(type+"Npip",     Form("rad::util::Count(%spid,211)",  type.data()) );
        rdf->Define(type+"Npim",     Form("rad::util::Count(%spid,-211)", type.data()) );
        rdf->Define(type+"NKp",      Form("rad::util::Count(%spid,321)",  type.data()) );
        rdf->Define(type+"NKm",      Form("rad::util::Count(%spid,-321)", type.data()) );
        rdf->Define(type+"Nele",     Form("rad::util::Count(%spid,11)",   type.data()) );
        rdf->Define(type+"Npos",     Form("rad::util::Count(%spid,-11)",  type.data()) );
        rdf->Define(type+"Npro",     Form("rad::util::Count(%spid,2212)", type.data()) );
        rdf->Define(type+"Nneutron", Form("rad::util::Count(%spid,2112)", type.data()) );
      }

      template<typename T> 
      inline void Resolution(T* const rdf, const std::string& var){
        rdf->Define(string("res_")+var, Form("%s-%s", (Truth()+var).data(), (Rec()+var).data() ));
      }

      template<typename T> 
      inline void ResolutionFraction(T* const rdf, const std::string& var){
        rdf->Define(string("res_")+var, Form("(%s-%s)/%s", (Truth()+var).data(), (Rec()+var).data(), (Truth()+var).data() ));
      }

   

    inline constexpr ResultType_t PdgToMass(int pdg) {
        // Safe constexpr absolute value
        int abs_pdg = (pdg < 0) ? -pdg : pdg;

        switch (abs_pdg) {
            case 22:   return 0.0;              // Photon
            case 11:   return 0.00051099895;    // Electron / Positron
            case 13:   return 0.10565837;       // Muon / Anti-muon
            case 211:  return 0.13957039;       // Charged Pion
            case 321:  return 0.493677;         // Charged Kaon
            case 2212: return 0.93827208;       // Proton / Anti-proton
            case 2112: return 0.93956541;       // Neutron / Anti-neutron
            
            // Light Nuclei
            case 45:         return 1.875612;   // Deuteron (CLAS12 old convention)
            case 1000010020: return 1.875612;   // Deuteron (Strict PDG convention)
            case 1000020030: return 2.80839;    // Helium-3
            case 1000020040: return 3.72738;    // Alpha
            
            default: return 0.0; // Fallback for unmatched/undefined PID
        }
    }
    

    ///////////////////////////////////////////////////////
    inline RVecResultType AssignMasses( const ROOT::RVecI &pid){
      auto n = pid.size();
      RVecResultType masses(n);
      for(size_t i=0;i<n;++i){
	masses[i] = PdgToMass(pid[i]);
      }
      return masses;
    }

  }
}
