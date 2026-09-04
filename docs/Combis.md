# RAD Framework: Combinatorics and Topology FAQ

## Introduction: The Combinatorial Problem
In high-energy physics, analyzing experimental data is rarely as simple as matching one detected track to one requested particle. The reality of detector data is messy. 

The need for "combinatorial events" arises when your event topology does not have a clear, unambiguous selection of particles. This can occur when:
1. There are more particles of a particular type in the detector than your reaction requires.
2. There are multiple identical particles in the final state (e.g., several photons or pions), and you don't know which ones came from which parent decay.
3. There is no clear Particle Identification (PID) for a track, so you want to test multiple mass hypotheses based purely on the track's charge.

In these cases, you could just reject the event entirely (losing valuable data)—or you can try **all possible combinations** of particles that fit your hypothesis. The RAD framework does the latter automatically, removing the need for you to write complex, deeply nested `for` loops that are often the source of hard-to-track bugs.

---

## What is a Combinatorial?
A "combinatorial" is simply one valid mathematical arrangement of detected tracks that satisfies your requested final state. 

If your topology requires exactly 1 electron and 1 proton, but the detector recorded **3 electrons and 2 protons** in a single collision event, the framework cannot know which tracks are the "correct" ones. Instead of guessing, it creates a Cartesian matrix of all possible pairings. In this example, it creates **$3 \times 2 = 6$ distinct combinations** for that single event.

Later in your analysis, you will calculate physics variables (like Missing Mass or Invariant Mass) for all 6 combinations and use selections (cuts) to filter out the false background.

---

## The Three Kinds of Combinatorics

### 1. The Over-Abundance Combinatorial
This occurs when the detector simply records more tracks of a specific species than your final state requires.
*   **Example A:** Your final state requires 1 electron, but the event has 3. -> *Generates 3 combinations.*
*   **Example B:** Your final state requires 1 electron and 1 proton, but the event has 3 electrons and 2 protons. -> *Generates 6 combinations.*

### 2. The Identical Species (Permutation) Combinatorial
This occurs when you have multiple particles of the same species in your final state, and you must evaluate permutations based on their parent decays. 
*   **Example A:** Your final state requires 2 photons from a single $\pi^0$ decay. -> *Generates 1 combination (order doesn't matter).*
*   **Example B:** Your final state requires 3 photons resulting from an $\omega \to \pi^0 + \gamma$ decay. The detector doesn't know which two photons belong to the $\pi^0$ and which is the standalone $\gamma$. -> *Generates 3 combinatorial permutations.*
*   **Example C:** Your final state requires $2\pi^+$ and $2\pi^-$ tracks from a complex decay chain. -> *Generates 4 combinations.*

### 3. The Ambiguous PID Combinatorial
This occurs when you do not have strict Particle Identification, so you make hypotheses based purely on track charge.
*   **Example:** You require a proton and a $\pi^+$. The detector records 2 positive tracks, but their PID is unclear. The framework will treat Track A as the proton and Track B as the pion (Combination 1), and then swap them, treating Track B as the proton and Track A as the pion (Combination 2).

---
## Dealing with Identical Particles: Symmetric Flags

### The Permutation Problem
When your final state has multiple identical particles (for example, two photons from a $\pi^0$ decay), the default combinatorial engine will treat them as distinct roles: "Photon 1" and "Photon 2". 

If the detector finds Track A and Track B, a standard combinatorial engine will test both permutations:
* **Combination 1:** Photon 1 = Track A, Photon 2 = Track B
* **Combination 2:** Photon 1 = Track B, Photon 2 = Track A

Physically, these are the exact same state. If you calculate the invariant mass of (Photon 1 + Photon 2) for both combinations, you will get the exact same number twice. This wastes CPU cycles and artificially inflates your background statistics (double-counting).

### The Solution: Symmetry Groups
To solve this, the RAD framework allows you to define **Symmetry Groups** (or Symmetric Flags) when setting up your reaction. 

When you tell the framework that Photon 1 and Photon 2 are "symmetric" (indistinguishable), the framework activates a special symmetry filter (`GenerateSymmetricCombinations`) during the event loop. 

**How it works under the hood:**
RAD enforces a strict mathematical ordering rule on the raw detector track indices. It demands that the track index assigned to Photon 1 must be strictly less than the track index assigned to Photon 2. 
* Track A (Index 3) and Track B (Index 7): $3 < 7$, so this combination is **kept**.
* Track B (Index 7) and Track A (Index 3): $7 \not< 3$, so this swapped combination is **instantly rejected**.

### When should I use Symmetric Flags?
You should group identical particles symmetrically **only if they come from the same parent decay** or if their specific assignment does not matter to your physics calculations. 

* **DO use symmetry:** For the two photons decaying from a $\pi^0$. You just need their combined 4-momentum; you don't care which is which. Using symmetry reduces the combinations from 2 to 1.
* **DO NOT use symmetry:** If you are analyzing a reaction with two $\pi^-$ tracks, but one comes from a $\Lambda^0$ decay and the other comes from a $\Xi^-$ decay. Because they have different parent particles, you *must* evaluate both permutations to figure out which pion belongs to which parent!

## How does the RAD framework handle combinatorics?
Older software packages handled combinatorics by literally copying and writing out a brand new, separate event to the output file for every combination. This bloated file sizes and made data incredibly slow to process.

The RAD framework is built on **ROOT's RDataFrame**, which uses modern, vectorized arrays (called `RVec`s). 

When you define your particles and call `MakeCombinations()`, RAD does *not* duplicate the event. Instead, it creates synchronized arrays of data for that single event. If an event has 6 valid combinations, RAD evaluates your physics math across all 6 combinations simultaneously and stores the result as an array of size 6. 

*   `Event 1, Missing Mass Array: [ 0.938, 1.210, 4.500, 2.100, 0.940, 1.111 ]`

By applying a selection mask (like requiring the Missing Mass to be near the proton mass), the framework mathematically filters out the false combinations instantly, leaving only the correct physics state.

---
## Coding Combinatorics & Symmetry

### How do I code combinatorics in my RAD script?
To generate combinatorial events, you must define your particle candidates and then explicitly trigger the engine. You do this by accessing the core reaction object from your `AnalysisManager`.

1. **Access the Reaction:** Retrieve the reaction object from your manager using `auto& reaction = mgr.Reaction();`.
2. **Define Candidates:** Tell the framework what tracks to look for using `SetParticleCandidates`. You must provide a string name, a truth-matching Role ID, an index filter, and the target bank column. 
   * Example: `reaction.SetParticleCandidates("pip", 5, rad::index::FilterIndices(211), {"rec_pid"});`.
   * *API Sugar:* If your specific reaction class supports it, you can use the shortened helper method: `reaction.SetParticleRecPID("pip", 211);`.
3. **Trigger the Engine:** Call `reaction.MakeCombinations();`. This expands the tracks into the Cartesian matrix and **must** be called before you add any processing streams or topology recipes.

### How do I code symmetric particle groups in my RAD script?
If your final state has identical, indistinguishable particles (like two photons from a single $\pi^0$ decay), you must tell the framework to treat them symmetrically to avoid double-counting the permutations. 

1. **Define Candidates:** Define both identical particles normally using `SetParticleCandidates` (e.g., define `"gamma1"` and `"gamma2"`).
2. **Apply Symmetry:** Group them together by calling `mgr.Reaction().SetSymmetryParticles({"gamma1", "gamma2"});`. 
3. **Trigger the Engine:** When you subsequently call `MakeCombinations()`, the engine detects the symmetry group and automatically applies an index-ordering filter during the Cartesian expansion, ensuring only unique physical states are kept.
---
## Why should I analyze data like this?
1. **Zero Nested Loops:** Writing 4-deep nested `for` loops to iterate over particle indices is tedious, unreadable, and highly prone to segmentation faults.
2. **Maximum Efficiency:** By pushing the combinatorics into C++ vector math (`RVec`), the framework analyzes thousands of combinations in a fraction of a millisecond.
3. **Data Preservation:** Instead of throwing away ambiguous events, you preserve your experiment's statistics by systematically evaluating and filtering the background.
