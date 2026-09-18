## RAD Framework FAQ

**What is a stream?**
A stream dictates the primary data source for your analysis. It defines the input prefix that the framework uses to read raw arrays from the TTree or HIPO file. The most common examples are `Rec()` for reconstructed detector data and `Truth()` for Monte Carlo generator data.

**What is a label for?**
A label (or suffix) represents an alternate hypothesis, systemic variation, or calibration applied to an existing data stream. By applying a label (like `_loose` or `_corrected`), the framework appends it to output variables so they do not overwrite the baseline calculations. This allows the new hypothesis to efficiently share the exact same underlying combinatorial track indices as the base stream without wasting CPU cycles.

**How do I analyze multiple topologies?**
The approach depends entirely on whether your topologies share the same base detector tracks:
*   **Different Base Tracks (Different Combinatorial Sizes):** If one topology requires 3 detected tracks ($e' \pi^+ \pi^-$) and another requires 4 tracks ($e' p \pi^+ \pi^-$), the combinatorial arrays will be different sizes. Because a manager is locked to a single combinatorial matrix, you **must** instantiate a separate `AnalysisManager` for each topology. You can efficiently copy your setup by cloning the base manager (`mgr.Clone("New_Topology")`).
*   **Same Base Tracks (Different Intermediate States):** If both topologies use the exact same base tracks but reconstruct different intermediate states (e.g., calculating a $\rho$ meson vs a $\Delta^{++}$ baryon from the same event), you can use the *same* `AnalysisManager`. You simply add a stream with a new label and use pattern matching to apply different kinematic recipes to each label.

**How do I analyze multiple different reactions?**
Similar to analyzing topologies with different base tracks, changing the physical reaction (e.g., Deeply Virtual Compton Scattering vs. Meson Production) fundamentally changes the required particle candidates and the Cartesian expansion. Therefore, each distinct physical reaction requires its own dedicated cloned `AnalysisManager` instance (`mgr.Clone("New_Reaction")`). By cloning the Analysis Manager we use same underlying dataframe thereby ensuring optimal data reading.

## The Strict Order of Operations

**What is the correct order to set up my analysis script?**
Because reaction aware dataframes utilize lazy evaluation and dynamically generate data columns, you must build your analysis in a very specific sequence. If you define elements out of order, you may trigger a "Reaction Class Error" (e.g., trying to define a particle without a registered stream) or find that your detector columns are missing.

Here is the strict initialization lifecycle required for your scripts:

*   **Step 1: Register Streams.** Immediately after initializing the `AnalysisManager`, you must explicitly register your data streams (e.g., using `SetupReconstructed()` or `mgr.AddStream()`) so the framework knows what data types and prefixes to expect.
*   **Step 2: Define Candidates.** Once the stream is registered, define your physical particle candidates using `SetParticleCandidates`. 
*   **Step 3: Make Combinations.** You must invoke `MakeCombinations()` to expand your track lists into the Cartesian matrix of all valid permutations.
*   **Step 4: Build Detectors (`clas12-rad` specific).** If you are using the `CLAS12DetectorBuilder`, it *must* be invoked immediately after `MakeCombinations()`. The builder relies on projection operations (`ROOT::VecOps::Take`) to map fundamental detector hit arrays directly onto the newly created combinatorial tracks, which requires the combinations to already exist.
*   **Step 5: Inject Recipes.** Finally, inject your topology, selection, and histogram recipes into the manager and execute the run.

## Truth Matching FAQ

**What is truth matching?**
Truth matching is the process of linking a reconstructed detector track back to the original Monte Carlo generated particle that created it. In the framework, this is done by assigning integer "Roles" to your particle candidates (e.g., `Role_ScatEle = 5`) which are their position or index in the mc particle record. The framework automatically traverses the underlying association arrays to map the reconstructed track to the exact MC particle that fulfills that role.

**How does it relate to streams?**
Reconstructed detector data and Monte Carlo generator data exist in completely separate streams (typically `rec_` and `tru_`). Truth matching acts as the bridge between them. When matching is enabled, the framework can safely pull variables from the truth stream and alias them into the reconstructed stream. This allows you to evaluate both sets of kinematics simultaneously within the same combinatorial matrix.

**How can I calculate resolutions?**
Resolutions (such as $\Delta P = P_{Rec} - P_{Tru}$) are calculated by subtracting the truth stream variables from the reconstructed stream variables. You can generate these automatically using the `AnalysisManager`'s cross-stream batch wrappers, such as `mgr.CrossStreamDifferences(Rec(), Truth(), {"proton"}, {"pmag", "theta"})`. The framework handles the heavy lifting of broadcasting the single truth value across your entire array of reconstructed combinations to prevent data-size crashes.

**How do I know which is the correct combi?**
Because the combinatorial engine expands a single event into multiple possible track permutations, you will often end up with a large array of candidate states. When truth matching is active, the framework automatically generates a boolean flag for each combination (typically accessed via `TruthMatchedCombi()`). By applying a selection cut requiring `isTruth == 1`, you instantly mask out the combinatorial background, leaving only the exact track combination that perfectly mirrors the Monte Carlo generated state.

## Plotting Data: 1D and 2D Histograms

**How do I make simple 1D and 2D histograms from my calculated variables?**

In RAD, you generate plots by defining a **Histogram Recipe** and injecting it into your `AnalysisManager`. This recipe uses the `rad::histo::Histogrammer` to safely book histograms for the combinatorial output.

### The "Drop the Prefix" Rule
When you ask the `Histogrammer` to plot a variable, **you must drop the stream prefix** (e.g., do not write `"rec_RhoMass"` or `"tru_RhoMass"`). Instead, you only provide the base name of the variable (e.g., `"RhoMass"`). 

**Why?** Because you configure the histogram recipe for a specific stream inside the `AnalysisManager` (e.g., `mgr.ConfigureHistograms(Rec(), histogram_recipe)`). The framework inherently knows which stream it is operating on and automatically resolves your base names to the correct, fully-qualified column names in the background. 

This is a powerful feature: it means you can write **one** histogram recipe and apply it to both your Reconstructed and Truth streams independently without rewriting any code!

### Example

Here is how to book basic 1D and 2D histograms using the base names generated by your topology recipe:

```cpp
// 1. Define your Histogram Recipe
auto histogram_recipe = [](rad::histo::Histogrammer& h) {
    
    // 1D Histogram
    // Syntax: Create(Name, Title; X-axis, Bins, Min, Max, TargetColumn)
    h.Create("hMissMass", "Missing Mass; MM [GeV/c^{2}]", 200, -0.5, 1.5, "MissMass");
             
    // 2D Histogram
    // Syntax: Create2D(Name, Title; X-axis; Y-axis, X-Bins, X-Min, X-Max, Y-Bins, Y-Min, Y-Max, X-Col, Y-Col)
    h.Create2D("hQ2_vs_W", "Q^{2} vs W; W [GeV]; Q^{2} [GeV^{2}]", 
               200, 0, 5, 200, 0, 10, "W", "Q2");
               
    // Example using a specific particle's property
    h.Create("hProtonP", "Proton Momentum; P [GeV/c]", 100, 0, 10, "proton_pmag");
};

// 2. Inject it into the Analysis Manager for your chosen stream(s)
// The manager will automatically prepend "rec_" to all the variables above!
mgr.ConfigureHistograms(rad::consts::data_type::Rec(), histogram_recipe);
```

*(Note: If you applied a Selection Recipe to this stream, the histograms will automatically evaluate the lazy mask and only plot the events that passed your cuts!)*


**How do I create a series of split histograms based on a category?**

Often in an analysis, you want to plot a distribution, but you need a separate histogram for every distinct bin or category of another variable (e.g., creating one momentum plot for each Detector Sector, or plotting missing mass in slices of $Q^2$).

To do this, you must first define a **Split Axis** using the `AddSplit` method on the `Histogrammer`. You must provide the exact binning limits (number of bins, minimum, and maximum) for your split variable so the framework knows how to segment the data.

*(Remember the Golden Rule: Just like your target variables, you must **drop the stream prefix** for your split variable. The manager handles the prefixes dynamically!)*

### Example

Here is how you define a split axis and book your histograms:

```cpp
auto histogram_recipe = [](rad::histo::Histogrammer& h) {
    
    // 1. Define the Split Axis
    // Syntax: AddSplit(Name, TargetColumn, Bins, Min, Max)
    // Example: Splitting by 6 CLAS12 sectors (discrete values 1 through 6)
    h.AddSplit("Sector", "proton_Sector", 6, 0.5, 6.5);
    
    // 2. Create your histograms
    // The framework will automatically generate an array of these histograms 
    // (e.g., hProtonP_Sector_1, hProtonP_Sector_2) mapped to the split axis above!
    h.Create("hProtonP", "Proton Momentum; P [GeV/c]", 100, 0, 10, "proton_pmag");
    
    h.Create2D("hQ2_vs_W", "Q^{2} vs W; W [GeV]; Q^{2} [GeV^{2}]", 
               200, 0, 5, 200, 0, 10, "W", "Q2");
};

// Inject the recipe into the manager as usual
mgr.ConfigureHistograms(rad::consts::data_type::Rec(), histogram_recipe);
```

## Applying Cuts and Filters to Combinatorial Events

**How do I apply a cut or filter to my combinatorial events?**

In standard `ROOT::RDataFrame` workflows, using a basic `Filter()` drops the entire event from the processing chain, which is bad for sideband analysis or when evaluating multiple combinatorial candidates within the same event. 

To solve this safely, RAD uses **Masking**. 
* Instead of deleting the row, `PhysicsSelection` creates a boolean column (a mask).
* Any connected Histograms and Trees will dynamically check this mask and only fill entries where `mask == true`.
* This allows you to run multiple, completely different selections (e.g., Signal vs Sideband) in a single pass without destroying the underlying data by using different labeled streams.

You can apply these masks by defining a selection recipe and injecting it into your `AnalysisManager`. Here is a general example:

```cpp
// 1. Define your Selection Recipe
auto selection_recipe = [](rad::PhysicsSelection& s) {
    // Add a minimum threshold cut (e.g., Track momentum > 0.3 GeV)
    s.AddCutMin("p_cut", "proton_pmag", 0.3); 
    
    // Add a range cut (e.g., Exclusivity cut on Missing Mass Squared around 0)
    s.AddCutRange("mm2_cut", "MissMass2", -0.05, 0.05);
};

// 2. Inject it into the Analysis Manager for your specific stream
mgr.ConfigureSelection(rad::consts::data_type::Rec(), selection_recipe);
```
## Understanding the Diagnostic Output

Because reaction aware dataframes utilize lazy evaluation, the computational graph is essentially a black box until the event loop executes. Running `mgr.PrintDiagnostics()` provides a transparent blueprint of your setup, but the output can be overwhelmingly dense. 

Here is how to decipher the most critical sections of that dump to troubleshoot your analysis:

### 1. Controlling the Clutter (Verbosity Levels)
Printing hundreds of lines per run is a massive help for debugging but is overkill for production. You can control this output by passing an integer to `PrintDiagnostics()`:

*   **`mgr.PrintDiagnostics(0);` (Silent):** Perfect for batch jobs and keeping log files completely clean.
*   **`mgr.PrintDiagnostics(1);` (Summary):** Prints a clean, high-level summary of your registered streams and input types.
*   **`mgr.PrintDiagnostics(2);` (Deep Dive):** Dumps the complete array maps, input aliases, and every registered output variable.

### 2. Combinatoric Structure
```text
Type: rec_
  Combo column: rec_reaction_combos__dnwtag
  Combinatoric generation: GenerateAllCombinations
--------------------------------------------------------------------------------
Candidates for type: rec_
--------------------------------------------------------------------------------
Particle                  Definition                                         
-----------------------------------------------------------------------------
pim                       ROOT::VecOps::Nonzero(rec_pid == -211 && rec_re... 
beam_ele                  [lambda: ]                                         
pip                       [lambda: rec_pid, rec_region]                      
proton                    [lambda: rec_pid]                                  
```
> **How to use it:** This section verifies exactly how the framework identifies candidate tracks for each particle. It confirms whether a simple `[lambda: rec_pid]` was used or if complex regional string cuts were successfully registered before the combinatorics engine runs.

### 3. Deciphering the Reaction Map 
```text
================================================================================
  Reaction Map (ParticleCreator)
================================================================================
Particle Name             Index      Type            
-----------------------------------------------------
beam_ele                  1          Input            
proton                    2          Input            
pip                       3          Input            
rho                       7          Created          
Miss                      9          Created          
```
> **How to use it:** To avoid slow, deep copies of particle objects, the framework separates the topology (indices) from the data (flat momentum arrays). The `ReactionMap` is simply a matrix of integer indices pointing to those flat arrays in memory. `Input` particles link directly to raw TTree data, while `Created` particles are dynamically calculated during the event loop.

### 4. Column Aliases 
```text
================================================================================
  Column Aliases (Input Particles)
================================================================================
Alias Name (In Analysis)            Source Column (In Tree)             
------------------------------------------------------------------------
rec_proton_base                     rec_proton                          
rec_pip_base                        rec_pip                             
rec_scat_ele_base                   rec_scat_ele                        
```
> **How to use it:** This section explicitly shows how your internal analysis names map to the raw source columns in the underlying TTree. This is invaluable for tracking down missing branch errors, proving exactly what the framework is searching for in your input file. 

### 5. Registered Calculations
```text
================================================================================
  Registered Calculations
================================================================================
Calculation Name                    Type            
----------------------------------------------------
rec_RhoMass_base                    IndexKernel     
rec_Whad_base                       IndexKernel     
rec_Q2_base                         MapKernel       
```
> **How to use it:** Lists all formulas and kernels queued for execution. Checking this ensures that your custom topology recipes and physics functions were correctly ingested by the KinematicsProcessor.

## Understanding Kernel Types: IndexKernel vs. MapKernel

When reviewing the diagnostics for the "reaction aware dataframes" framework, you will notice that registered calculations are categorized as either an `IndexKernel` or a `MapKernel`. This refers to how the underlying C++ function accesses the combinatorial particle tracks inside the event loop.

**What is an `IndexKernel`?**
An `IndexKernel` is heavily optimized for speed. When you register this type of calculation, you explicitly pass the list of particle names it requires (e.g., `{"pip", "pim"}`). The framework pre-resolves these string names into fixed integer indices *before* the event loop ever starts. During execution, the kernel simply receives a fast, pre-compiled array of integers (`RVecIndices`) pointing directly to the required tracks. Standard kinematic shortcuts like `Mass`, `Pt`, and `ParticleTheta` use this under the hood.

**What is a `MapKernel`?**
A `MapKernel` is optimized for flexibility. Instead of pre-resolving specific tracks, this kernel receives the entire `RVecIndexMap` (the complete topology map of the combination) during the hot loop. The kernel itself performs lookups on the fly using constants or names (e.g., `map[rad::consts::OrderScatEle()]`). This is ideal for global event variables—like $Q^2$, $W$, or defining custom Helicity frames—where you might need to flexibly query the whole topology at once without passing long lists of specific particles.

**Are there other types?**
Currently, no. Within the internal `KineCalculation` registry of the "reaction aware dataframes" framework, `KernelType::Map` and `KernelType::Index` are the only two defined execution paths for custom physics variables.

### 6. Variables Registered for Output
```text
================================================================================
  Variables Registered for Output
================================================================================
ID    Base Name                      Full Column Name                         
------------------------------------------------------------------------------
0     W_sys_px                       rec_W_sys_px_base                        
44    RhoMass                        rec_RhoMass_base                         
47    MissMass2                      rec_MissMass2_base                       
```
> **How to use it:** The framework lists every single variable it plans to write out alongside its exact final column name (e.g., translating `RhoMass` into `rec_RhoMass_base`). This saves you from having to guess your branch names when opening a `TBrowser`.

## Writing Custom Physics Kernels

**How do I write and use a custom physics calculation?**
To add new physics variables (like specific scattering angles, custom invariant masses, or helicity frames), you write a standalone C++ function and register it inside your Topology Recipe. 

Because the framework uses lazy evaluation, you never call this function directly to get a return value. Instead, you register it with the processor, which seamlessly injects it into the execution graph as a new named column.

**Important:** The C++ implementation of your function **must** be placed in a separate header file (or a separate cell, if using a Jupyter Notebook) from your main analysis script. 

### Step 1: The Implementation (Separate File / Cell)
Your custom function must accept the combinatorial index map (`RVecIndexMap`) and the flat data arrays (like `px`, `py`) as arguments. 

```cpp
// ---------------------------------------------------------
// File: MyPhysicsKernels.h (or executed in a prior Notebook cell)
// ---------------------------------------------------------
#include "CommonDefines.h"

namespace rad {
  namespace physics {
      
    // The function signature MUST return an RVecResultType (an array of doubles)
    inline RVecResultType CalcMyAngle(const RVecIndexMap& map, 
                                      const RVecResultType& px, 
                                      const RVecResultType& py) 
    {
      // 1. Look up the specific particle indices for the current topology
      auto idx_ele = map[rad::consts::OrderScatEle()];
      auto idx_pro = map[rad::consts::OrderBaryon()];

      // 2. Perform the math across the vectorized arrays
      auto result = px[idx_ele] * py[idx_pro]; 
      
      return result;
    }
    
  } // namespace physics
} // namespace rad
```

### Step 2: The User Script (Main Analysis Cell)
Once the compiler knows about your function, you can register it in your Topology Recipe and immediately use the resulting column name in your cuts or histograms.

```cpp
// ---------------------------------------------------------
// File: MainAnalysis.C
// ---------------------------------------------------------
#include "MyPhysicsKernels.h"

// 1. Inject it via the Topology Recipe Lambda
auto topology_recipe = [](Processor& p) {
    // Registers the calculation as a new column named "MyAngleVal"
    p.RegisterCalc("MyAngleVal", rad::physics::CalcMyAngle);
};

// 2. Call the registered calculation by its string name in downstream recipes
auto histogram_recipe = [](rad::histo::Histogrammer& h) {
    // Use the "MyAngleVal" column to populate a 1D histogram
    h.Create("hMyAngle", "Custom Angle; #theta [rad]", 100, 0, 3.14, "MyAngleVal");
};

auto selection_recipe = [](rad::PhysicsSelection& s) {
    // Or use it to generate a lazy mask (e.g., keep events where angle > 0.5)
    s.AddCutMin("AngleCut", "MyAngleVal", 0.5);
};
```