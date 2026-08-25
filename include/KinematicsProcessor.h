/**
 * @file KinematicsProcessor.h
 * @brief The Computational Engine for combinatorial kinematic analysis.
 * @details
 * This class orchestrates the combinatorial logic. It:
 * 1. Consumes the combinatorial indices from the `ParticleCreator`.
 * 2. Runs the "Event Loop" (via `operator()`) to generate 4-vectors for every valid combination.
 * 3. Provides an API to register calculations (Mass, Pt, etc.) on these combinations.
 * 4. Prepares flattened data columns for `SnapshotCombi`.
 */

#pragma once

#include "ConfigReaction.h"
#include "StringUtilities.h"
#include "RVecHelpers.h"
#include "ParticleCreator.h"
#include "KinematicsDispatch.h" // Contains DefineKinematicsProcessor template
#include "KineCalculation.h" 
#include "BasicKinematics.h"
#include "ParticleModifier.h" 

#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <iomanip>

namespace rad {

  using ROOT::RVec;
  using RVecRVecD = ROOT::RVec<ROOT::RVecD>;
  using RVecRVecI = ROOT::RVec<ROOT::RVecI>;

  /**
   * @class KinematicsProcessor
   * @brief Manages the kinematic calculations and RDF graph construction.
   * @details 
   * This class serves as the bridge between the physics topology (defined in `ParticleCreator`)
   * and the RDataFrame execution graph. It manages a registry of calculations 
   * so that the `AnalysisManager` knows exactly which columns to write to disk.
   */
  class KinematicsProcessor {

  public:
    /// Type alias for the complex nested vector structure [Combination][Component][Particle]
    using CombiOutputVec_t = RVec<RVec<RVecResultType>>;

    // =================================================================================
    // Lifecycle & Initialization
    // =================================================================================

    /**
     * @brief Constructor.
     * @param cr Pointer to the Reaction configuration.
     * @param prefix The data prefix (e.g., "rec_", "tru_").
     * @param suffix Optional suffix for systematics (e.g., "_sysUp").
     */
    KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix = "");
    
    virtual ~KinematicsProcessor() = default;

    /** * @brief Initializes the processor and defines columns in the RDF graph. 
     * @details This function is idempotent (can be called multiple times safely).
     * It triggers the `ParticleCreator` setup, defines the main execution kernel via
     * `DefineKinematicsProcessor`, and runs definitions for all registered calculations.
     */
    void Init();
    
    // =================================================================================
    // Core Execution (Functor)
    // =================================================================================

    /**
     * @brief The main computational kernel.
     * @details Generates Px, Py, Pz, M for every combination defined in `indices`.
     * Applies Pre/Post modifiers (Momentum corrections, Energy loss, etc.).
     * * @tparam Tp Type of Momentum input (RVecF or RVecD).
     * @tparam Tm Type of Mass input (RVecF or RVecD).
     * @param indices The combinatorial indices [Particle][Combination].
     * @param px Input Px vector.
     * @param py Input Py vector.
     * @param pz Input Pz vector.
     * @param m Input Mass vector.
     * @param aux_pre_d Auxiliary Doubles (Pre-creation).
     * @param aux_pre_i Auxiliary Ints (Pre-creation).
     * @param aux_post_d Auxiliary Doubles (Post-creation).
     * @param aux_post_i Auxiliary Ints (Post-creation).
     * @return The nested structure containing kinematics for all combinations.
     */
    template<typename Tp, typename Te> 
    CombiOutputVec_t operator()(const RVecIndices& indices, 
                                const Tp& px, const Tp& py, const Tp& pz, const Te& e,
                                const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
                                const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const;

    /** * @brief Defines flattened columns for Px, Py, Pz, E for every particle. 
     * @details 
     * Essential for `SnapshotCombi`. This creates scalar RVecs (one entry per combination)
     * effectively flattening the nested `CombiOutputVec_t` structure.
     * It registers these variables (e.g., "ele_px") so `AnalysisManager` knows to save them.
     */
    void DefineNewComponentVecs();

    // =================================================================================
    // Configuration Helpers
    // =================================================================================

    /** @brief Defines the Auxiliary data pack columns in the RDF graph. */
    void DefineAux();

    /** @brief Define a generic particle group override. */
    void SetGroup(const std::string& groupName, const ROOT::RVec<std::string>& particles);

    /** @brief Override the standard "Mesons" group. */
    void SetMesonParticles(const ROOT::RVec<std::string>& particles);

    /** @brief Override the standard "Baryons" group. */
    void SetBaryonParticles(const ROOT::RVec<std::string>& particles);

    /** @brief Re-target this processor to a new Reaction (used internally). */
    void SetReaction(ConfigReaction* cr) { 
      _reaction = cr; 
      _creator.SetReaction(cr); 
    }
    
    // =================================================================================
    // Accessors
    // =================================================================================

    ParticleCreator& Creator();
    const ParticleCreator& Creator() const;
    
    ParticleModifier& PreModifier();
    ParticleModifier& PostModifier();

    ConfigReaction* Reaction() const;
    std::string GetSuffix() const;
    std::string GetPrefix() const;
    
    /** @brief Construct a full column name: prefix + base + suffix. */
    std::string FullName(const std::string& baseName) const;
    std::string CheckedFullName(const std::string& baseName) const;

    /** * @brief Get list of all variables defined via RegisterCalc/Mass/Pt etc.
     * @details This is used by `AnalysisManager::Snapshot` to auto-detect columns.
     */
    const ROOT::RVec<std::string>& GetDefinedNames() const { 
      return _registered_vars; 
    }

    // =================================================================================
    // Definition API & Calculation Registration
    // =================================================================================
    
    /** * @brief Registers a calculation based on the RVecIndexMap (e.g. Mass, Missing Mass).
     * @param name The output variable name (base name only).
     * @param func The kernel function to execute.
     */
    void RegisterCalc(const std::string& name, KineCalculation::MapKernel func);

    /** * @brief Registers a calculation based on explicit particle indices. 
     * @param name The output variable name.
     * @param func The kernel function.
     * @param particles List of particle names needed by the kernel.
     */
    void RegisterCalc(const std::string& name, KineCalculation::IndexKernel func, ROOT::RVec<ParticleNames_t> particles);

    // --- Legacy String-Based Definition ---
    
    void Define(const std::string& name, const std::string& func);
    void Define(const std::string& name, const std::string& func, const ROOT::RVec<ParticleNames_t>& particles);

    /** @brief Helper to define an arbitrary lambda kernel in RDF. */
    template <typename Lambda>
    void DefineKernel(const std::string& name, Lambda&& func);

    void DefineTruthFlag();
    
    /**
     * @brief Passes raw detector arrays through the combinatorial engine to the flat output tree.
     * @param pName The target particle (e.g., "ele").
     * @param rawArray The raw data array (e.g., "rec_charge").
     * @param compSuffix The name to append to the particle for the output (e.g., "_charge").
     */
    void PassThrough(const std::string& pName, const std::string& rawArray, const std::string& compSuffix);
    // =================================================================================
    // Physics Shortcuts
    // =================================================================================

    void Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    void Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    void Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});
    void Energy(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg={});

    void ParticleTheta(const ParticleNames_t& particles);
    void ParticlePhi(const ParticleNames_t& particles);
    void ParticleP(const ParticleNames_t& particles);
    void ParticleEta(const ParticleNames_t& particles);
    
    void PrintReactionMap() const;


    //////////////////////
    ///Diagnostics
    ///\brief Print registered calculations.
    void PrintCalculations() const;

    ///\brief Print registered output variables.
    void PrintRegisteredVariables() const;

    ///\brief Print group overrides.
    void PrintGroupOverrides() const;

    ///\brief Print comprehensive processor diagnostics.
    void PrintProcessorDiagnostics() const;

  private:
    ConfigReaction* _reaction = nullptr;
    std::string _prefix; 
    std::string _suffix; 
    
    bool _isInitialized = false; 
    
    ParticleCreator _creator;
    ParticleModifier _preModifier;
    ParticleModifier _postModifier;

    ROOT::RVec<KineCalculation> _calculations;
    ROOT::RVec<std::string> _registered_vars; // Registry of variable names for Snapshot
    
    struct GroupOverride {
        std::string name;
        ROOT::RVec<std::string> particles;
    };
    ROOT::RVec<GroupOverride> _groupOverrides;

    void ApplyGroupOverrides();
    
    struct PassThroughDef {
        std::string pName;
        std::string rawArray;
        std::string compSuffix;
    };
    ROOT::RVec<PassThroughDef> _passThroughs; // Stores deferred passthrough requests

    // =========================================================================
    // Thread-Local Pre-allocated Buffers & Index Caches
    // =========================================================================
    
    /** 
     * @brief Pre-allocated memory buffers for kinematics.
     * @details Marked mutable to allow in-place updates within the const operator().
     */
    mutable RVecResultType _bufferPx;
    mutable RVecResultType _bufferPy;
    mutable RVecResultType _bufferPz;
    mutable RVecResultType _bufferM;
    mutable RVecResultType _bufferE;
    
    /**
     * @brief Pre-computed SIMD energy buffer for the flat detector arrays.
     */
    mutable RVecResultType _baseE;

    /** 
     * @brief Pre-allocated memory buffers for auxiliary data caches.
     * @details Marked mutable to allow in-place updates within the const operator().
     */
    mutable AuxCacheD _cachePreD;
    mutable AuxCacheI _cachePreI;
    mutable AuxCacheD _cachePostD;
    mutable AuxCacheI _cachePostI;

    /** 
     * @brief Fast caches for topology indexing to prevent redundant lookups 
     *        and improve inner-loop cache locality.
     */
    mutable std::vector<size_t> _ipartiCache;
    mutable std::vector<int> _ogIndicesCache;
    
  };

  // =================================================================================
  // IMPLEMENTATION: KinematicsProcessor
  // =================================================================================

  inline KinematicsProcessor::KinematicsProcessor(ConfigReaction* cr, const std::string& prefix, const std::string& suffix) 
    : _reaction{cr}, _prefix{prefix}, _suffix{suffix}, _creator{cr, prefix, suffix} 
  {}

  inline void KinematicsProcessor::Init() {
    if (_isInitialized) return; 
    _isInitialized = true;

    _reaction->ValidateType(_prefix);
    // Clear the registry to ensure no stale/duplicate entries from config phase
    _registered_vars.clear();

    ApplyGroupOverrides();
    Creator().InitMap(); 
    _preModifier.Init(Creator());
    _postModifier.Init(Creator());
    DefineAux();

    // Calls KinematicsDispatch::DefineKinematicsProcessor
    // IMPORTANT: This internally calls this->DefineNewComponentVecs() !
    DefineKinematicsProcessor(*_reaction, *this, _prefix);

    
    for(auto& calc : _calculations) {
      calc.Define(this); 
    }

    // Execute deferred PassThroughs AFTER the clear!
    for(const auto& pt : _passThroughs) {
        std::string outNameBase = pt.pName + pt.compSuffix; 
        std::string colName = FullName(outNameBase);  
        std::string candCol = CheckedFullName(pt.pName); 
        
        // Use FastTake to avoid the ROOT condition vector crash!
        _reaction->Define(colName, "ROOT::VecOps::Take(" + pt.rawArray + ", " + candCol + ")");
        
        _registered_vars.push_back(outNameBase);
    }
  }

  inline void KinematicsProcessor::ApplyGroupOverrides() {
      for(const auto& group : _groupOverrides) {
          Creator().OverrideGroup(group.name, group.particles); 
      }
  }

  inline void KinematicsProcessor::SetGroup(const std::string& groupName, const ROOT::RVec<std::string>& particles) {
      _groupOverrides.push_back({groupName, particles});
  }

  inline void KinematicsProcessor::SetMesonParticles(const ROOT::RVec<std::string>& particles) {
      SetGroup(rad::consts::Mesons(), particles);
  }

  inline void KinematicsProcessor::SetBaryonParticles(const ROOT::RVec<std::string>& particles) {
      SetGroup(rad::consts::Baryons(), particles);
  }

  inline void KinematicsProcessor::RegisterCalc(const std::string& name, KineCalculation::MapKernel func) {
      _calculations.emplace_back(name, func);
  }

  inline void KinematicsProcessor::RegisterCalc(const std::string& name, KineCalculation::IndexKernel func, ROOT::RVec<ParticleNames_t> particles) {
    
      _calculations.emplace_back(name, func, particles);
  }

  inline void KinematicsProcessor::DefineAux() {
      auto define_pack = [&](const std::string& name, const ROOT::RVec<std::string>& cols, bool is_int) {
          if(cols.empty()) {
              if(is_int) _reaction->Define(name, [](){ return RVecRVecI{}; }, {});
              else       _reaction->Define(name, [](){ return RVecRVecD{}; }, {});
              return;
          }
          _reaction->Define(name, rad::util::createPackVectorString(cols));
      };
      define_pack(_prefix + "aux_pre_d" + _suffix + DoNotWriteTag(), _preModifier.GetAuxDoubleCols(), false);
      define_pack(_prefix + "aux_pre_i" + _suffix + DoNotWriteTag(), _preModifier.GetAuxIntCols(), true);
      define_pack(_prefix + "aux_post_d" + _suffix + DoNotWriteTag(), _postModifier.GetAuxDoubleCols(), false);
      define_pack(_prefix + "aux_post_i" + _suffix + DoNotWriteTag(), _postModifier.GetAuxIntCols(), true);
  }
  
  inline ParticleCreator& KinematicsProcessor::Creator() { return _creator; }
  inline const ParticleCreator& KinematicsProcessor::Creator() const { return _creator; }
  
  inline ParticleModifier& KinematicsProcessor::PreModifier() { return _preModifier; }
  inline ParticleModifier& KinematicsProcessor::PostModifier() { return _postModifier; }

  inline ConfigReaction* KinematicsProcessor::Reaction() const { return _reaction; }
  inline std::string KinematicsProcessor::GetSuffix() const { return _suffix; }
  inline std::string KinematicsProcessor::GetPrefix() const { return _prefix; }
  
  inline std::string KinematicsProcessor::FullName(const std::string& baseName) const { 
    return _prefix + baseName + _suffix; 
  }
  inline std::string KinematicsProcessor::CheckedFullName(const std::string& baseName) const {
    auto fullName = _prefix + baseName + _suffix;
    if(_reaction->ColumnExists(fullName)==false){
      fullName = _prefix + baseName;
         if(_reaction->ColumnExists(fullName)==false){
	   fullName = baseName;
	 }
	 if(_reaction->ColumnExists(fullName)==false){
	   throw std::runtime_error("KinematicsProcessor::FullName, Column '" + fullName + "' does not exist with any known perfic or suffix.");
	 }
    }
    return fullName; 
  }

  // --- Core Operator ---
  template<typename Tp, typename Tm> 
  inline KinematicsProcessor::CombiOutputVec_t KinematicsProcessor::operator()(
        const RVecIndices& indices, const Tp& px, const Tp& py, const Tp& pz, const Tm& m,
        const RVecRVecD& aux_pre_d, const RVecRVecI& aux_pre_i,
        const RVecRVecD& aux_post_d, const RVecRVecI& aux_post_i) const 
  {
    const auto Ncomponents = 4; // x, y, z, m
    const auto Nparticles0 = indices.size(); // Number of input particles
    const auto Nparticles = Nparticles0 + _creator.GetNCreated(); 
    
    if (Nparticles == 0) return CombiOutputVec_t(Ncomponents); 
          
    const auto Ncombis = indices[0].size(); 
    CombiOutputVec_t result(Ncombis, RVec<RVecResultType>(Ncomponents, RVecResultType(Nparticles)));

    // Resize thread-local buffers once per initialization if particle count changes
    if (_bufferPx.size() != Nparticles) {
        _bufferPx.resize(Nparticles);
        _bufferPy.resize(Nparticles);
        _bufferPz.resize(Nparticles);
        _bufferM.resize(Nparticles);
        _bufferE.resize(Nparticles);
    }

    // Resize fast index caches and pre-compute destination indices ONCE per event
    if (_ipartiCache.size() != Nparticles0) {
        _ipartiCache.resize(Nparticles0);
        _ogIndicesCache.resize(Nparticles0);
        for (size_t ip = 0; ip < Nparticles0; ++ip) {
            _ipartiCache[ip] = _creator.GetReactionIndex(ip);
        }
    }

    // Resize auxiliary caches if their outer size changed or inner size does not match Nparticles
    if (_cachePreD.size() != aux_pre_d.size() || (!_cachePreD.empty() && _cachePreD[0].size() != Nparticles)) {
        _cachePreD = AuxCacheD(aux_pre_d.size(), ROOT::RVecD(Nparticles));
    }
    if (_cachePreI.size() != aux_pre_i.size() || (!_cachePreI.empty() && _cachePreI[0].size() != Nparticles)) {
        _cachePreI = AuxCacheI(aux_pre_i.size(), ROOT::RVecI(Nparticles));
    }
    if (_cachePostD.size() != aux_post_d.size() || (!_cachePostD.empty() && _cachePostD[0].size() != Nparticles)) {
        _cachePostD = AuxCacheD(aux_post_d.size(), ROOT::RVecD(Nparticles));
    }
    if (_cachePostI.size() != aux_post_i.size() || (!_cachePostI.empty() && _cachePostI[0].size() != Nparticles)) {
        _cachePostI = AuxCacheI(aux_post_i.size(), ROOT::RVecI(Nparticles));
    }

    const ResultType_t invalid_val = consts::InvalidEntry<ResultType_t>();

    // ========================================================================
    // SIMD PRE-COMPUTATION
    // ========================================================================
    // Calculate energy for all raw input tracks exactly ONCE per event.
    // The compiler will auto-vectorize this flat loop across the SoA arrays.
    if (_baseE.size() != px.size()) {
        _baseE.resize(px.size());
    }
    for(size_t i = 0; i < px.size(); ++i) {
        _baseE[i] = std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] + m[i]*m[i]);
    }

    //cout<< "KinematicsProcessor::operator() "<<Nparticles0<<" "<<Nparticles<<" "<<Ncombis<<endl;
    // ========================================================================
    // OPTIMIZED COMBINATORIAL LOOP
    // ========================================================================
    for (size_t icombi = 0; icombi < Ncombis; ++icombi) {
      
      // Initialize standard kinematics to InvalidEntry once per combination 
      std::fill(_bufferPx.begin(), _bufferPx.end(), invalid_val);
      std::fill(_bufferPy.begin(), _bufferPy.end(), invalid_val);
      std::fill(_bufferPz.begin(), _bufferPz.end(), invalid_val);
      std::fill(_bufferM.begin(),  _bufferM.end(),  invalid_val);
      std::fill(_bufferE.begin(),  _bufferE.end(),  invalid_val);

      // Replicate original AuxCache zero-initialization on creation
      for (auto& vec : _cachePreD) std::fill(vec.begin(), vec.end(), 0.0);
      for (auto& vec : _cachePreI) std::fill(vec.begin(), vec.end(), 0);
      for (auto& vec : _cachePostD) std::fill(vec.begin(), vec.end(), 0.0);
      for (auto& vec : _cachePostI) std::fill(vec.begin(), vec.end(), 0);

      // 1. Core Kinematics & Source Index Extraction
      for (size_t ip = 0; ip < Nparticles0; ++ip) {
        
        // Use fast O(1) cache instead of redundant map lookups
        size_t iparti = _ipartiCache[ip];                
        const int og_idx = indices[ip][icombi];     
        
        // Save source index for cache-friendly aux processing
        _ogIndicesCache[ip] = og_idx;

        _bufferPx[iparti] = px[og_idx];
        _bufferPy[iparti] = py[og_idx];
        _bufferPz[iparti] = pz[og_idx];
        _bufferM[iparti]  = m[og_idx];

        // Important: Only for REAL input particles
        // i.e. M2>0 and relatavistic condition:
        // E2 = p2 + m2 HOLDS!
        // O(1) lookup: Fetch pre-calculated energy instead of doing scalar sqrt()
        _bufferE[iparti]  = _baseE[og_idx];
      }

      // 2. Cache-Friendly Auxiliary Processing
      // Iterating arrays (v) first, then tracks (ip), preserves memory contiguity
      for(size_t v=0; v<aux_pre_d.size(); ++v) {
          for(size_t ip=0; ip<Nparticles0; ++ip) {
              _cachePreD[v][_ipartiCache[ip]] = aux_pre_d[v][_ogIndicesCache[ip]];
          }
      }
      for(size_t v=0; v<aux_pre_i.size(); ++v) {
          for(size_t ip=0; ip<Nparticles0; ++ip) {
              _cachePreI[v][_ipartiCache[ip]] = aux_pre_i[v][_ogIndicesCache[ip]];
          }
      }

      _preModifier.Apply(_bufferPx, _bufferPy, _bufferPz, _bufferE, _cachePreD, _cachePreI);
      
      _creator.ApplyCreation(_bufferPx, _bufferPy, _bufferPz, _bufferE);
      
      _postModifier.Apply(_bufferPx, _bufferPy, _bufferPz, _bufferE, _cachePostD, _cachePostI);
 
      result[icombi][OrderX()] = _bufferPx;
      result[icombi][OrderY()] = _bufferPy;
      result[icombi][OrderZ()] = _bufferPz;
      result[icombi][OrderE()] = _bufferE;
    }
    // auto particle_names = _creator.GetParticleNames();
    // for(const auto& pname:particle_names){
    //   cout<<" "<<pname<<" "<<_creator.GetReactionIndex(pname);
    // }
    // cout<<endl;
    // cout<<"KinematicPRoessor "<<endl<<result<<endl;
    return result;
  }

 // --- Snapshot Support ---
  inline void KinematicsProcessor::DefineNewComponentVecs() {

       auto particle_names = _creator.GetParticleNames();
       
       const ROOT::RVec<std::pair<std::string, int>> components = {
         {"_px", 0}, {"_py", 1}, {"_pz", 2}, {"_e", 3}
       };
       
       std::string resultColName = _prefix + consts::KineComponents() + _suffix;

       for (const auto& pName : particle_names) {
           size_t idx = _creator.GetReactionIndex(pName);

           for (const auto& comp : components) {
               std::string compSuffix = comp.first; // e.g. "_px"
               auto compIdx = comp.second;
        
               // 1. Register base name: "Y_px"
               _registered_vars.push_back(pName + compSuffix); 

               // 2. Construct Column Name: "rec_Y_px_0"
               // Logic: Prefix + Particle + CompSuffix + StreamSuffix
               std::string colName = _prefix + pName + compSuffix + _suffix;

               _reaction->Define(colName, 
                 [idx,compIdx](const CombiOutputVec_t& res) {
                   ROOT::RVec<double> out(res.size());
                   for(size_t i=0; i<res.size(); ++i) {
                     out[i] = res[i][compIdx][idx];
                   }
                   return out;
                 }, 
                 {resultColName}
                 );
           }
       }
  }
 
  // --- Definitions ---
  
  /**
     * @brief Generates a stream-specific Truth Match vector.
     * @details 
     * If truth matching is set up, this creates a boolean vector (e.g. "rec_isTruth_loose")
     * corresponding exactly to the combinations in this stream ("rec_Indices_loose").
     */
    // void DefineTruthFlag() {
    //     // 1. Only generate for Reconstruction streams (skip "tru_" streams)
    //     // Heuristic: Check if prefix starts with "rec"
    //   if(_prefix.find(Rec()) == std::string::npos) return;

    //     // 2. Construct names
    //     std::string baseName = consts::TruthMatchedCombi(); // "isTruth"
    //     std::string outputCol = FullName(baseName);         // "rec_isTruth_loose"
    //     std::string indicesCol = Creator().GetMapName();    // "rec_Indices_loose"

    //     // 3. Define the column via the Reaction's Truth Registry
    //     // This function (in ConfigReaction) handles looking up the Truth IDs
    //     // and matching them against the provided indicesCol.
    //     // It returns true if successful (i.e., truth matching is configured).
        
    //     bool success = _reaction->DefineTruthMatch(outputCol, indicesCol);
        
    //     if(success) {
    //         // 4. Register for Snapshotting
    //         // AnalysisManager will see "isTruth", convert it to "rec_isTruth_loose",
    //         // and save it to the tree.
    //         _registered_vars.push_back(baseName);
    //     }
    // }

  inline void KinematicsProcessor::PassThrough(const std::string& pName, const std::string& rawArray, const std::string& compSuffix) {
    _passThroughs.push_back({pName, rawArray, compSuffix});
  }
 
  inline void KinematicsProcessor::Define(const std::string& name, const std::string& func) {
      std::string colName = FullName(name); 
      _reaction->Define(colName, util::createFunctionCallStringFromVec("rad::util::ApplyKinematics", 
                                       {func, Creator().GetMapName(), (_prefix + consts::KineComponents() + _suffix)}));
      _registered_vars.push_back(name);
  }
  
  template <typename Lambda>
  inline void KinematicsProcessor::DefineKernel(const std::string& name, Lambda&& func) {
   
      ROOT::RDF::ColumnNames_t cols = { Creator().GetMapName(), _prefix + consts::KineComponents() + _suffix };
      auto apply_func = [func](const RVecIndexMap& map, const ROOT::RVec<ROOT::RVec<RVecResultType>>& comps){
        return rad::util::ApplyKinematics(func, map, comps);
      };
      _reaction->Define(FullName(name), apply_func, cols);
      
      // calculated variables should be registered!
      _registered_vars.push_back(name);
  }

 
  // --- Shortcuts ---
  inline void KinematicsProcessor::Mass(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    RegisterCalc(name, rad::FourVectorMassCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Mass2(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    RegisterCalc(name, rad::FourVectorMass2Calc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Pt(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    RegisterCalc(name, rad::FourVectorPtCalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::Energy(const std::string& name, const ParticleNames_t& particles_pos, const ParticleNames_t particles_neg) {
    RegisterCalc(name, rad::FourVectorECalc<rad::RVecResultType, rad::RVecResultType>, {particles_pos, particles_neg});
  }
  inline void KinematicsProcessor::ParticleTheta(const ParticleNames_t& particles) {
    //here we actually perform a loop over combies
    //so we need to call a dunction which returns
    //a single entry each time.
    //rad::ThreeVectorTheta has an overwite which only uses the
    //first entry of the first entry in the given RVecIndices list
    for(const auto& p: particles){
      RegisterCalc(p+"_theta", rad::ThreeVectorTheta, {{p}});
    }
  }
 inline void KinematicsProcessor::ParticlePhi(const ParticleNames_t& particles) {
    for(const auto& p: particles){
      RegisterCalc(p+"_phi", rad::ThreeVectorPhi, {{p}});
    }
  }
 inline void KinematicsProcessor::ParticleP(const ParticleNames_t& particles) {
    for(const auto& p: particles){
      RegisterCalc(p+"_pmag", rad::ThreeVectorMag, {{p}});
    }
  }
 inline void KinematicsProcessor::ParticleEta(const ParticleNames_t& particles) {
    for(const auto& p: particles){
      RegisterCalc(p+"_eta", rad::ThreeVectorEta, {{p}});
    }
 }
 
  inline void KinematicsProcessor::PrintReactionMap() const {
    std::cout << "\n=== KinematicsProcessor [" << _prefix << "] " << _suffix << " Reaction Map ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Particle Name" << "Index" << std::endl;
    std::cout << std::string(30, '-') << std::endl;
    for(const auto& [name, idx] : _creator.GetIndexMap()) {
      std::cout << std::left << std::setw(20) << name << idx << std::endl;
    }
    std::cout << "========================================\n" << std::endl;
  }

  // =================================================================================
  // IMPLEMENTATION: KineCalculation::Define
  // =================================================================================
  
  inline void KineCalculation::Define(KinematicsProcessor* processor) {
    if (_kern_type == KernelType::Map) {
          processor->DefineKernel(_name, _mapFunc);
      } 
      else if (_kern_type == KernelType::Index) {
          RVecIndices resolved_indices;
          for(const auto& group : _particles) {
              Indices_t idxs;
              for(const auto& pname : group) {
                  idxs.push_back(processor->Creator().GetReactionIndex(pname));
              }
              resolved_indices.push_back(idxs);
          }

          auto func_ptr = _indexFunc; 
          auto adapter = [resolved_indices, func_ptr](const RVecIndexMap&, 
                                                      const RVecResultType& px, const RVecResultType& py, 
                                                      const RVecResultType& pz, const RVecResultType& e) 
          {
              return func_ptr(resolved_indices, px, py, pz, e);
          };

          processor->DefineKernel(_name, adapter);
      }
  }

} // namespace rad
#include "Diagnostics_KinematicProcessor.hxx"
