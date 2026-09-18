# RAD: Reaction Analysis & Design

**A Declarative, Vectorized RDataFrame Framework for High-Performance Physics Analysis**

[![Documentation](https://img.shields.io/badge/docs-rad-blue.svg)](https://dglazier.github.io/rad/)

RAD is a header-only C++ framework designed for high-intensity physics experiments. It abstracts complex multi-particle combinatorics, kinematics, and Monte Carlo truth matching into a high-level declarative API, leveraging ROOT's **RDataFrame** for parallel execution. 

For comprehensive guides and API references, please visit our [Official Documentation](https://dglazier.github.io/rad/).

---

## ⚡ Key Features

* **Vectorized Core (SoA):** Replaces heavy objects (like `TLorentzVector`) with Structure-of-Arrays (SoA) `ROOT::RVec<double>` data layouts, enabling compiler auto-vectorization and high cache locality.
* **Recipe-Based Architecture:** Compartmentalizes user code into safe, multi-threading-friendly C++ lambdas using the `AnalysisManager`.
* **Zero-Copy Combinatorics:** Evaluates multi-particle topologies via integer lookup arrays (`ReactionMap`) without deep-copying memory.
* **Lock-Free Parallel I/O:** Safely writes multi-threaded TTrees using `SnapshotCombi`, utilizing thread-local buffers to prevent segmentation faults during parallel execution.
* **Lazy Masking:** Uses non-destructive boolean mask columns for physics cuts instead of standard filtering, preserving array alignment for complex sideband analyses.

---

## 🚀 Quick Start (Recipe-Based API)

RAD utilizes an `AnalysisManager` to orchestrate execution, ensuring stability across multiple CPU threads. Here is how you define a complete topology using declarative recipes and a generic data source (e.g., HepMC3).

```cpp
#include "AnalysisManager.h"
#include "HepMCElectro.h"

void ProcessHepMCZCombi() {
    // 1. Enable Parallel Processing
    ROOT::EnableImplicitMT();

    // 2. Initialize the Manager & Data Source
    rad::AnalysisManager<rad::HepMCElectro, rad::KinematicsProcElectro> mgr{
        "Analysis", "hepmc3_tree", "data_file.root"
    };
    
    auto& reaction = mgr.Reaction();
    reaction.SetupMC();
    reaction.SetBeamElectronIndex(0); 
    
    // 3. Define Candidates & Map Combinations
    reaction.SetParticleCandidates("ele", 2 /*Role*/, rad::index::FilterIndices(11), {"pid"});
    reaction.SetParticleCandidates("pos", 3 /*Role*/, rad::index::FilterIndices(-11), {"pid"});
    reaction.MakeCombinations();
    
    mgr.AddStream(rad::data_type::Rec());

    // 4. The Topology Recipe (Injected into all CPU threads)
    auto topology_recipe = [](rad::KinematicsProcElectro& p) {
        // Define Composite Particles (J/psi -> e+ e-)
        p.Creator().Sum("Jpsi", {{"ele", "pos"}}); 
        
        // Calculate Physics      
        p.Mass("M_Jpsi", {"Jpsi"});
        p.Q2();                 
    };

    // 5. Inject, Snapshot, and Run!
    mgr.ConfigureKinematics(topology_recipe);
    mgr.Snapshot("output.root"); 
    mgr.Run();      
}
```

---

## 🏗 Modular Architecture

RAD relies on a strict, hierarchical pipeline to ensure high performance and lazy evaluation:

* **The Orchestrator (`AnalysisManager`):** Manages data streams, lambda recipes, and lifecycle execution.
* **The Interface (`ConfigReaction`):** Wraps `ROOT::RDataFrame` and handles combinatorics and data extraction.
* **The Engine (`KinematicsProcessor`):** Operates on raw SoA columns (Px, Py, Pz) using SIMD instructions. It features modular injection points:
    * **Injector:** For data normalization.
    * **Creator:** For topology and grouping.
    * **Modifier:** For calibrations and auxiliary data.

---

## 📁 Parallel Output & Masking

RAD provides robust tools for saving your analysis-level data without crashing the event loop:

* **SnapshotCombi:** Pre-allocates `thread_local` buffers and guards `gDirectory` to safely write TTrees from multiple threads.
* **Lazy Flattening:** Automatically flattens nested combinatorial structures (Event $\to$ Combination) into an analysis-friendly N-Tuple (one row per candidate). Event-level scalars are automatically broadcasted.
* **Lazy Masking:** Instead of using standard `Filter()` nodes which drop events entirely, RAD creates boolean mask columns to preserve array alignment for complex sideband analyses.