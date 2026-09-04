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
When you ask the `Histogrammer` to plot a variable, **you must drop the stream prefix** (e.g., do not write `"rec_RhoMass"` or `"tru_RhoMass"`). Instead, you only provide the base name of the variable (e.g., `"RhoMass"`)[cite: 5]. 

**Why?** Because you configure the histogram recipe for a specific stream inside the `AnalysisManager` (e.g., `mgr.ConfigureHistograms(Rec(), histogram_recipe)`). The framework inherently knows which stream it is operating on and automatically resolves your base names to the correct, fully-qualified column names in the background[cite: 5]. 

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

## Writing Custom Physics Kernels

**How do I write and use a custom physics calculation?**
To add new physics variables (like specific scattering angles, custom invariant masses, or helicity frames), you write a standalone C++ function and register it inside your Topology Recipe[cite: 9]. 

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
Once the compiler knows about your function, you can register it in your Topology Recipe and immediately use the resulting column name in your cuts or histograms[cite: 9].

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