/**
 * @file AnalysisManager.h
 * @brief High-level orchestration of the RAD analysis chain.
 * @details
 * The AnalysisManager is the central entry point for the user. It manages:
 * 1. Data Streams (Rec, Truth, Systematics).
 * 2. Configuration Recipes (Kinematics, Cuts, Histograms).
 * 3. Execution phases (Initialization, Snapshotting, Running).
 * 
 * It uses a "Three-Pass Initialization" strategy to handle dependencies:
 * - Pass 1: Define Kinematics (Variables).
 * - Pass 2: Compile Selections (Cuts depending on Variables).
 * - Pass 3: Book Histograms (Plots depending on Variables and Cuts).
 */

#pragma once

#include "ConfigReaction.h"
#include "KinematicsProcessor.h"
#include "PhysicsSelection.h"
#include "Histogrammer.h"
#include "CommonDefines.h"
#include "Indicing.h"
#include "Diagnostics.h"

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include <functional>

namespace rad {

// =====================================================================
// Internal Detail Namespace (Hides data structs from public API)
// =====================================================================
namespace detail {

    /**
     * @struct CrossStreamDef
     * @brief Configuration storage for deferred cross-stream calculations.
     */
    struct CrossStreamDef {
        std::string stream1;
        std::string stream2;
        std::string varBaseName;
        std::string suffix;
        std::string outName;
        bool isSum; // true for Addition (+), false for Subtraction (-)
    };

    /**
     * @struct CrossStreamAliasDef
     * @brief Aliasing cross stream variables into streams
     */
    struct CrossStreamAliasDef {
        std::string sourceStream;
        std::string targetStream;
        std::string varBaseName;
        std::string suffix;
    };

    /**
     * @struct AnalysisStream
     * @brief Holds the components for a single analysis pipeline.
     */
    template <typename ReactionClass, typename ProcessorClass>
    struct AnalysisStream {
        std::string fullName;   // Unique ID (e.g. "rec_loose")
        std::string source;     // Data Source (e.g. "rec")
        std::string label;      // Variation Label (e.g. "loose")
        
        std::unique_ptr<ProcessorClass> kine;
        std::unique_ptr<PhysicsSelection> sel;
        std::unique_ptr<histo::Histogrammer> hist;
        
        bool hasHistograms = false; // Tracks if a histogram recipe was applied

        // Dumb Constructor: All string logic is handled securely by AnalysisManager
        AnalysisStream(ReactionClass* reaction, const std::string& src, const std::string& lbl,
                       const std::string& fName, const std::string& inPrefix, const std::string& outSuffix) 
            : fullName(fName), source(src), label(lbl),
              kine(std::make_unique<ProcessorClass>(reaction, inPrefix, outSuffix)),
              sel(std::make_unique<PhysicsSelection>(*kine)),
              hist(std::make_unique<histo::Histogrammer>(*kine, sel.get())) 
        {}
    };

    /**
     * @brief Helper to check if a stream matches a pattern.
     * @details Exact matching on Source, Label, or FullName.
     */
    template <typename ReactionClass, typename ProcessorClass>
    inline bool StreamMatches(const AnalysisStream<ReactionClass, ProcessorClass>& stream, const std::string& pattern) {
        if(pattern == stream.fullName) return true; 
        if(pattern == stream.source)   return true; 
        if(pattern == stream.label)    return true; 
        return false;
    }

} // namespace detail


  /**
   * @class AnalysisManager
   * @brief The main driver class for RAD analysis.
   * @tparam ReactionClass The concrete reaction class (e.g. ePICReaction).
   * @tparam ProcessorClass The concrete processor class (e.g. KinematicsProcessor).
   */
  template <typename ReactionClass, typename ProcessorClass>
  class AnalysisManager {

  public:
    // Recipe Signatures
    using KineRecipe  = std::function<void(ProcessorClass&)>;
    using SelRecipe   = std::function<void(PhysicsSelection&)>;
    using HistoRecipe = std::function<void(histo::Histogrammer&)>;

    // =====================================================================
    // Setup
    // =====================================================================
    
    /**
     * @brief Constructor.
     * @param name Name of the analysis (used for output filenames).
     * @param treeName Input TTree name.
     * @param fileGlob Input file path/pattern.
     */
    AnalysisManager(const std::string& name, const std::string& treeName, const std::string& fileGlob);
    
    /**
     * @brief Constructor (multiple input files).
     */
    AnalysisManager(const std::string& name, const std::string& treeName, const ROOT::RVec<std::string>& files);
    
    /**
     * @brief Constructor (from existing reaction).
     */
    AnalysisManager(const std::string& name, const ReactionClass& baseReaction);
    
    /**
     * @brief Constructor clone method.
     */
    AnalysisManager Clone(const std::string& newName, bool copyStreams = false) const;
    
    /** @brief Sets and creates the output directory. */
    void SetOutputDir(const std::string& dir);

    // =====================================================================
    // Stream Management
    // =====================================================================

    /**
     * @brief Adds a specific analysis stream.
     * @param dataSource The data prefix to read (e.g. "rec", "tru").
     * @param label Analysis label (e.g. "loose", "tight", "sysUp"). Leave empty for default.
     */
    void AddStream(const std::string& dataSource, const std::string& label = "");

    /** @return Reference to the underlying Reaction object. */
    ReactionClass& Reaction();

    // =====================================================================
    // Configuration (Pattern Matching Enabled)
    // =====================================================================

    /** @brief Apply a Kinematics recipe to ALL active streams. */
    void ConfigureKinematics(KineRecipe recipe);
    /** @brief Apply a Kinematics recipe to streams matching the pattern. */
    void ConfigureKinematics(const std::string& pattern, KineRecipe recipe);

    /** @brief Apply a Selection recipe to ALL active streams. */
    void ConfigureSelection(SelRecipe recipe);
    /** @brief Apply a Selection recipe to streams matching the pattern. */
    void ConfigureSelection(const std::string& pattern, SelRecipe recipe);

    /** @brief Apply a Histogram recipe to ALL active streams. */
    void ConfigureHistograms(HistoRecipe recipe);
    /** @brief Apply a Histogram recipe to streams matching the pattern. */
    void ConfigureHistograms(const std::string& pattern, HistoRecipe recipe);

    // =====================================================================
    // Execution
    // =====================================================================

    /**
     * @brief Initializes the entire analysis chain using a 3-Pass system.
     */
    void Init();

    /**
     * @brief Snapshots data to separate flat TTrees (One per Stream).
     * @param addCols Additional columns to save (e.g. Truth Scalars) in ALL trees.
     * @param filenameBase Optional override for base filename.
     */
    void Snapshot(const ROOT::RDF::ColumnNames_t& addCols = {}, const std::string& filenameBase = "");
  
    /**
     * @brief Runs the analysis and writes histograms.
     * @param suffix Suffix for the histogram file (default "_Hist.root").
     */
    void Run(const std::string& suffix = "_Hist.root");

    /**
     * @brief Print comprehensive diagnostics for the entire analysis setup.
     */
    void PrintDiagnostics(int verbosity = 1);

    // =====================================================================
    // Cross-Stream Operations
    // =====================================================================

    /**
     * @brief Queues the difference between variables in two streams (Stream1 - Stream2).
     */
    void CrossStreamDifference(const std::string& stream1, const std::string& stream2,
                               const std::string& varBaseName, const std::string& suffix = "",
                               const std::string& outName = "");

    /**
     * @brief Queues the sum of variables in two streams (Stream1 + Stream2).
     */
    void CrossStreamSum(const std::string& stream1, const std::string& stream2,
                        const std::string& varBaseName, const std::string& suffix = "",
                        const std::string& outName = "");

    /**
     * @brief Batch queues differences between variables in two streams (Stream1 - Stream2).
     * @param tracks List of particle names (e.g., {"proton", "pip"}).
     * @param vars List of variable suffixes (e.g., {"pmag", "theta"}).
     */
    void CrossStreamDifferences(const std::string& stream1, const std::string& stream2,
                                const std::vector<std::string>& tracks, 
                                const std::vector<std::string>& vars,
                                const std::string& suffix = "");

    /**
     * @brief Batch queues sums of variables in two streams (Stream1 + Stream2).
     */
    void CrossStreamSums(const std::string& stream1, const std::string& stream2,
                         const std::vector<std::string>& tracks, 
                         const std::vector<std::string>& vars,
                         const std::string& suffix = "");
    
    /**
     * @brief Aliases a variable from one stream into another (e.g., tru_P -> rec_tru_P).
     */
    void CrossStreamAlias(const std::string& sourceStream, const std::string& targetStream,
                          const std::string& varBaseName, const std::string& suffix = "");

  private:
    ReactionClass _reaction;
    std::string _name;
    std::string _outputDir;
    bool _initialized = false;
    std::string _primaryStream;
    std::map<std::string, detail::AnalysisStream<ReactionClass, ProcessorClass>> _streams;
    
    std::vector<detail::CrossStreamDef> _crossStreams;
    std::vector<detail::CrossStreamAliasDef> _crossAliases;

    // --- Private Helpers ---
    std::string MakePath(const std::string& filename);
    ROOT::RDF::ColumnNames_t CollectStreamColumns(const ProcessorClass& kine);
    void ProcessCrossStreams();

    // Snapshot Helpers
    ROOT::RDF::ColumnNames_t CollectGlobalTruth();
    std::string ConstructSnapshotFilename(const std::string& filenameBase, const std::string& streamFullName);
    std::string GenerateFallbackMask(detail::AnalysisStream<ReactionClass, ProcessorClass>& stream);
  };

  // ===========================================================================
  // IMPLEMENTATION: AnalysisManager
  // ===========================================================================

  template <typename R, typename P>
  inline AnalysisManager<R,P>::AnalysisManager(const std::string& name,
                                               const std::string& treeName,
                                               const std::string& fileGlob) 
    : _reaction(treeName, fileGlob), _name(name) {}
  
  template <typename R, typename P>
  inline AnalysisManager<R,P>::AnalysisManager(const std::string& name,
                                               const std::string& treeName,
                                               const ROOT::RVec<std::string>& files)
    : _reaction(treeName, files), _name(name) {}

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::SetOutputDir(const std::string& dir) {
      _outputDir = dir;
      if (!_outputDir.empty() && !std::filesystem::exists(_outputDir)) {
          std::filesystem::create_directories(_outputDir);
      }
  }
  
  template <typename R, typename P>
  inline AnalysisManager<R,P>::AnalysisManager(const std::string& name, const R& baseReaction)
    : _reaction(baseReaction), _name(name) {}
  
  template <typename R, typename P>
  inline AnalysisManager<R,P> AnalysisManager<R,P>::Clone(const std::string& newName, bool copyStreams) const {
      AnalysisManager<R,P> out(newName, _reaction);
      if (copyStreams) {
          for (auto& [key, s] : _streams) out.AddStream(s.source, s.label);
      }
      out.SetOutputDir(_outputDir);
      return out;
  }
  
  // --- Stream Management ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::AddStream(const std::string& dataSource, const std::string& label) {
      std::string key = dataSource;
      if(!label.empty()) key += "_" + label;

      if(_streams.find(key) != _streams.end()) {
          std::cerr << "AnalysisManager Warning: Stream '" << key << "' already exists." << std::endl;
          return;
      }

      // Isolate string formatting from the constructor logic
      std::string fullName = dataSource;
      std::string outputSuffix = "";
      if(!label.empty()) {
          fullName += label;
          outputSuffix = "_" + label;
      } else {
          if(!fullName.empty()) fullName.pop_back(); // Remove trailing underscore for default names
      }
      
      std::string inputPrefix = dataSource;
      if(inputPrefix.back() != '_') inputPrefix += "_";

      _streams.try_emplace(key, &_reaction, dataSource, label, fullName, inputPrefix, outputSuffix);

      if(_primaryStream.empty()) _primaryStream = key;
  }

  template <typename R, typename P>
  inline R& AnalysisManager<R,P>::Reaction() { return _reaction; }

  // --- Configuration ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureKinematics(KineRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.kine) recipe(*stream.kine); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureKinematics(const std::string& pattern, KineRecipe recipe) {
      bool found = false;
      for(auto& [key, stream] : _streams) {
          if(detail::StreamMatches(stream, pattern) && stream.kine) {
              recipe(*stream.kine);
              found = true;
          }
      }
      if(!found) std::cerr << "AnalysisManager Warning: No streams matched pattern '" << pattern << "'" << std::endl;
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureSelection(SelRecipe recipe) {
      for(auto& [key, stream] : _streams) { if(stream.sel) recipe(*stream.sel); }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureSelection(const std::string& pattern, SelRecipe recipe) {
      bool found = false;
      for(auto& [key, stream] : _streams) {
          if(detail::StreamMatches(stream, pattern) && stream.sel) {
              recipe(*stream.sel);
              found = true;
          }
      }
      if(!found) std::cerr << "AnalysisManager Warning: No streams matched pattern '" << pattern << "'" << std::endl;
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureHistograms(HistoRecipe recipe) {
      for(auto& [key, stream] : _streams) { 
          if(stream.hist) {
              recipe(*stream.hist); 
              stream.hasHistograms = true;
          }
      }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::ConfigureHistograms(const std::string& pattern, HistoRecipe recipe) {
      bool found = false;
      for(auto& [key, stream] : _streams) {
          if(detail::StreamMatches(stream, pattern) && stream.hist) {
              recipe(*stream.hist);
              stream.hasHistograms = true;
              found = true;
          }
      }
      if(!found) std::cerr << "AnalysisManager Warning: No streams matched pattern '" << pattern << "'" << std::endl;
  }

  // --- Execution ---

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Init() {
      if(_initialized) return;
      if(_primaryStream.empty()) throw std::runtime_error("[AnalysisManager] No streams defined! Call AddStream().");
       
      // PASS 1: Initialize Kinematics (Create Variables)
      for(auto& [key, stream] : _streams) stream.kine->Init();
      
      // PROCESS CROSS-STREAMS (After Kinematics, Before Selections/Histograms)
      ProcessCrossStreams();
        
      // PASS 2: Compile Selections (Create Masks)
      for(auto& [key, stream] : _streams) { if(stream.sel) stream.sel->Init(); }

      // PASS 3: Initialize Histograms (Book Actions)
      for(auto& [key, stream] : _streams) { if(stream.hist) stream.hist->Init(); }

      _initialized = true;
  }

  // --- Snapshot Helpers ---

  template <typename R, typename P>
  inline ROOT::RDF::ColumnNames_t AnalysisManager<R,P>::CollectGlobalTruth() {
      using namespace rad::consts::data_type;
      ROOT::RDF::ColumnNames_t globalTruthCols;
      for(auto& [key, stream] : _streams) {
          if(stream.source.find(Truth()) == 0) {
              auto tCols = CollectStreamColumns(*stream.kine);
              globalTruthCols.insert(globalTruthCols.end(), tCols.begin(), tCols.end());
          }
      }
      return globalTruthCols;
  }

  template <typename R, typename P>
  inline std::string AnalysisManager<R,P>::ConstructSnapshotFilename(const std::string& filenameBase, const std::string& streamFullName) {
      if(filenameBase.empty()) return _name + "_" + streamFullName + "_Tree.root";
      
      auto dotPos = filenameBase.find_last_of(".");
      if(dotPos != std::string::npos) {
          std::string base = filenameBase.substr(0, dotPos);
          std::string ext  = filenameBase.substr(dotPos);
          return base + "_" + streamFullName + ext;
      }
      return filenameBase + "_" + streamFullName + ".root";
  }

  template <typename R, typename P>
  inline std::string AnalysisManager<R,P>::GenerateFallbackMask(detail::AnalysisStream<R, P>& stream) {
      std::string mask = stream.sel->GetMaskColumn();
      if(mask.empty()) {
          auto pNames = stream.kine->Creator().GetParticleNames();
          if(!pNames.empty()) {
              mask = stream.kine->GetPrefix() + "Analysis_AllIndices" + stream.kine->GetSuffix();
              std::string ref = stream.kine->FullName(pNames[0] + "_" + consts::NamePx());
              if(!_reaction.ColumnExists(mask)) {
                  _reaction.Define(mask, [](const ROOT::RVecD& r){ 
                      return rad::util::EnumerateIndicesFrom(0, r.size()); 
                  }, {ref});
              }
          }
      }
      return mask;
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Snapshot(const ROOT::RDF::ColumnNames_t& addCols, const std::string& filenameBase) {
      using namespace rad::consts::data_type;
      Init();
      
      ROOT::RDF::ColumnNames_t globalTruthCols = CollectGlobalTruth();

      for(auto& [key, stream] : _streams) {
          std::string specificFile = ConstructSnapshotFilename(filenameBase, stream.fullName);
          std::string finalPath = MakePath(specificFile);

          auto cols = addCols;
          cols.push_back("rdfentry_");
          
          if((stream.source.find(Rec()) == 0) && (stream.source.find(Truth()) == 0)) {
              cols.push_back(consts::TruthMatchedCombi());
              cols.push_back(Rec() + consts::TruthMatchedCombi());
          }
          
          auto streamCols = CollectStreamColumns(*stream.kine);
          cols.insert(cols.end(), streamCols.begin(), streamCols.end());

          // Append Truth columns to Reconstruction trees
          if(stream.source.find(Rec()) == 0) {
              cols.insert(cols.end(), globalTruthCols.begin(), globalTruthCols.end());
          }

          std::string mask = GenerateFallbackMask(stream);

          std::cout << "[AnalysisManager] Snapshotting stream '" << stream.fullName 
                    << "' to " << finalPath << " (with Mask: " << mask << ")" << std::endl;
          _reaction.BookSnapshotCombi(finalPath, "tree", cols, mask);
      }
  }
  
  template <typename R, typename P>
  inline void AnalysisManager<R,P>::Run(const std::string& suffix) {
      Init();
      _reaction.TriggerSnapshots();
      
      if(!_streams.empty()) {
          for(auto& [key, stream] : _streams) {
              if(!stream.hasHistograms) {
                  std::cout << "[AnalysisManager] Warning: No histograms configured for stream '" 
                            << stream.fullName << "'. Skipping file creation." << std::endl;
                  continue;
              }

              std::string fname = _name + "_" + stream.fullName + suffix;
              std::string finalPath = MakePath(fname);
              
              stream.hist->File(finalPath, "RECREATE");
              std::cout << "[AnalysisManager] Wrote " << stream.fullName << " histograms to " << finalPath << std::endl;
          }
      }
  }

  template <typename R, typename P>
  inline void AnalysisManager<R,P>::PrintDiagnostics(int verbosity) {
      if (verbosity <= 0) return; 
      Init(); 
      
      diag::DiagnosticsPrinter::PrintSectionHeader("ANALYSIS MANAGER DIAGNOSTICS", '=', 90);
      diag::DiagnosticsPrinter::PrintKeyValue("Analysis Name", _name);
      diag::DiagnosticsPrinter::PrintKeyValue("Output Directory", _outputDir.empty() ? "[none]" : _outputDir);
      diag::DiagnosticsPrinter::PrintKeyValue("Primary Stream", _primaryStream);
      diag::DiagnosticsPrinter::PrintBlank();
      
      _reaction.PrintReactionDiagnostics();
      diag::DiagnosticsPrinter::PrintBlank();
      
      if (!_crossStreams.empty() || !_crossAliases.empty()) {
          diag::DiagnosticsPrinter::PrintSectionHeader("CROSS-STREAM CONFIGURATIONS", '=', 90);
          
          if (!_crossStreams.empty()) {
              diag::DiagnosticsPrinter::PrintSubsection("Cross-Stream Math (Resolutions / Sums)", '-', 80);
              for (const auto& def : _crossStreams) {
                  std::string opStr = def.isSum ? "Sum (+)" : "Difference (-)";
                  std::string defaultPrefix = def.isSum ? "sum_" : "res_";
                  std::string out = def.outName.empty() ? def.stream1 + defaultPrefix + def.varBaseName + def.suffix : def.outName;
                  std::string equation = def.stream1 + def.varBaseName + def.suffix + " " + (def.isSum ? "+" : "-") + " " + def.stream2 + def.varBaseName + def.suffix;
                  
                  diag::DiagnosticsPrinter::PrintKeyValue("Output Column", out);
                  diag::DiagnosticsPrinter::PrintKeyValue("Equation", equation);
                  diag::DiagnosticsPrinter::PrintBlank();
              }
          }
          
          if (!_crossAliases.empty()) {
              diag::DiagnosticsPrinter::PrintSubsection("Cross-Stream Aliases", '-', 80);
              for (const auto& alias : _crossAliases) {
                  std::string sourceCol = alias.sourceStream + alias.varBaseName + alias.suffix;
                  std::string targetCol = alias.targetStream + alias.sourceStream + alias.varBaseName + alias.suffix;
                  
                  diag::DiagnosticsPrinter::PrintKeyValue("Target Alias", targetCol);
                  diag::DiagnosticsPrinter::PrintKeyValue("Source Column", sourceCol);
                  diag::DiagnosticsPrinter::PrintBlank();
              }
          }
      }

      diag::DiagnosticsPrinter::PrintSectionHeader("STREAM CONFIGURATIONS", '=', 90);
      std::cout << "Registered Streams: " << _streams.size() << std::endl;
      
      for (const auto& [key, stream] : _streams) {
          diag::DiagnosticsPrinter::PrintSubsection("Stream: " + stream.fullName, '-', 80);
          diag::DiagnosticsPrinter::PrintKeyValue("Source", stream.source);
          diag::DiagnosticsPrinter::PrintKeyValue("Label", stream.label.empty() ? "[default]" : stream.label);
          diag::DiagnosticsPrinter::PrintBlank();
          
          if(stream.kine && verbosity >= 2) stream.kine->PrintProcessorDiagnostics(); 
      }
  }

  // --- Cross-Stream Math & Aliasing ---
  
  template<typename R, typename P>
  inline void AnalysisManager<R, P>::CrossStreamDifference(
      const std::string& stream1, const std::string& stream2,
      const std::string& varBaseName, const std::string& suffix, const std::string& outName) 
  {
      _crossStreams.push_back({stream1, stream2, varBaseName, suffix, outName, false});
  }

  template<typename R, typename P>
  inline void AnalysisManager<R, P>::CrossStreamSum(
      const std::string& stream1, const std::string& stream2,
      const std::string& varBaseName, const std::string& suffix, const std::string& outName) 
  {
      _crossStreams.push_back({stream1, stream2, varBaseName, suffix, outName, true});
  }

  template<typename R, typename P>
  inline void AnalysisManager<R, P>::CrossStreamDifferences(
      const std::string& stream1, const std::string& stream2,
      const std::vector<std::string>& tracks, const std::vector<std::string>& vars,
      const std::string& suffix) 
  {
      for (const auto& track : tracks) {
          for (const auto& var : vars) {
              CrossStreamDifference(stream1, stream2, track + "_" + var, suffix, "");
          }
      }
  }

  template<typename R, typename P>
  inline void AnalysisManager<R, P>::CrossStreamSums(
      const std::string& stream1, const std::string& stream2,
      const std::vector<std::string>& tracks, const std::vector<std::string>& vars,
      const std::string& suffix) 
  {
      for (const auto& track : tracks) {
          for (const auto& var : vars) {
              CrossStreamSum(stream1, stream2, track + "_" + var, suffix, "");
          }
      }
  }
  
  template<typename R, typename P>
  inline void AnalysisManager<R, P>::CrossStreamAlias(
      const std::string& sourceStream, const std::string& targetStream,
      const std::string& varBaseName, const std::string& suffix) 
  {
      _crossAliases.push_back({sourceStream, targetStream, varBaseName, suffix});
  }
  
  template<typename R, typename P>
  inline void AnalysisManager<R, P>::ProcessCrossStreams() 
  {
      // 1. Process Math (Resolutions)
      for (const auto& def : _crossStreams) {
          std::string col1 = def.stream1 + def.varBaseName + def.suffix;
          std::string col2 = def.stream2 + def.varBaseName + def.suffix;
          
          std::string opStr = def.isSum ? "Sum" : "Difference";
          std::string defaultPrefix = def.isSum ? "sum_" : "res_";
          std::string out = def.outName.empty() ? def.stream1 + defaultPrefix + def.varBaseName + def.suffix : def.outName;
          std::string mathOp = def.isSum ? " + " : " - ";

          if (!_reaction.ColumnExists(col1)) throw std::runtime_error("\n[RAD CROSS-STREAM ERROR] Missing Stream 1 Column: " + col1 + "\n");
          if (!_reaction.ColumnExists(col2)) throw std::runtime_error("\n[RAD CROSS-STREAM ERROR] Missing Stream 2 Column: " + col2 + "\n");

          // ROOT::RVecD cast prevents temporary Expression Templates from causing reallocation segfaults during event loop
          std::string formula = col2 + ".size() > 0 ? ROOT::RVecD(" + col1 + mathOp + col2 + "[0]) : ROOT::RVecD(" + col1 + " * 0.0 - 9999.0)";
          _reaction.Define(out, formula);
      }

      // 2. Process Aliases
      for (const auto& alias : _crossAliases) {
          std::string sourceCol = alias.sourceStream + alias.varBaseName + alias.suffix;
          std::string targetCol = alias.targetStream + alias.sourceStream + alias.varBaseName + alias.suffix;
          std::string refCol = alias.targetStream + alias.varBaseName + alias.suffix;

          if (!_reaction.ColumnExists(sourceCol)) throw std::runtime_error("\n[RAD ALIAS ERROR] Source column missing: " + sourceCol + "\n");
          if (!_reaction.ColumnExists(refCol)) throw std::runtime_error("\n[RAD ALIAS ERROR] Target reference column missing: " + refCol + "\n");

          // ROOT::RVecD cast prevents temporary Expression Templates from causing reallocation segfaults during event loop
          std::string formula = sourceCol + ".size() > 0 ? ROOT::RVecD(" + refCol + " * 0.0 + " + sourceCol + "[0]) : ROOT::RVecD(" + refCol + " * 0.0 - 9999.0)";
          _reaction.Define(targetCol, formula);
      }
  }

  // --- Private Helpers ---

  template <typename R, typename P>
  inline std::string AnalysisManager<R,P>::MakePath(const std::string& filename) {
      if(!_outputDir.empty()) {
          return (std::filesystem::path(_outputDir) / filename).string();
      }
      return filename;
  }

  template <typename R, typename P>
  inline ROOT::RDF::ColumnNames_t AnalysisManager<R,P>::CollectStreamColumns(const P& kine) {
      ROOT::RDF::ColumnNames_t cols;
      for(const auto& var : kine.GetDefinedNames()) {
          cols.push_back(kine.FullName(var));
      }
      std::string sigCol = kine.GetPrefix() + rad::consts::TruthMatchedCombi();
      if(_reaction.ColumnExists(sigCol)) {
          cols.push_back(sigCol);
      }
      return cols;
  }

} // namespace rad
