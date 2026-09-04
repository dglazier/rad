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