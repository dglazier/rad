
Markdown


# Combinatorial API Refactor for rad/epic-rad

## Motivation
- Old code: Kinematic functions (like Q2, ElectroCMDecay) operate on a single particle assignment per event (using one `config::RVecIndexMap`).
- Desired: Support all valid combinations per event, apply kinematic calculations to each, and return `RVec` of results.

---

## Generic Combinatorial Helper Pattern

### Helper Template

```cpp
namespace rad {
namespace util {

// Generic wrapper for any single-combo kinematic function
template <typename F, typename... Args>
auto ApplyCombinationsGeneric(
    F&& singleComboFunc,
    const ROOT::RVec<config::RVecIndexMap>& combos,
    Args&&... args
) {
    using T = std::invoke_result_t<F, const config::RVecIndexMap&, Args...>;
    ROOT::RVec<T> results(combos.size());
    for (size_t i = 0; i < combos.size(); ++i) {
        results[i] = std::forward<F>(singleComboFunc)(combos[i], std::forward<Args>(args)...);
    }
    return results;
}

} // namespace util
} // namespace rad
Usage Example
Suppose you have:

C++


template<typename Tp, typename Tm>
double Q2(const config::RVecIndexMap& react,
          const ROOT::RVec<Tp>& px,
          const ROOT::RVec<Tp>& py,
          const ROOT::RVec<Tp>& pz,
          const ROOT::RVec<Tm>& m);
To compute Q2 for all combinations:

C++


auto Q2s = rad::util::ApplyCombinationsGeneric(
    rad::electro::Q2<float, float>,
    combos, px, py, pz, m
);
Similarly, for ElectroCMDecay (returns XYZVector):

C++


auto angles = rad::util::ApplyCombinationsGeneric(
    rad::electro::ElectroCMDecay<float, float>,
    combos, px, py, pz, m
);
This pattern works for any kinematic function with the same argument style.

RDataFrame Integration Example
C++


df.Define("Q2s", [](const ROOT::RVec<config::RVecIndexMap>& combos,
                    const ROOT::RVec<float>& px,
                    const ROOT::RVec<float>& py,
                    const ROOT::RVec<float>& pz,
                    const ROOT::RVec<float>& m) {
    return rad::util::ApplyCombinationsGeneric(rad::electro::Q2<float, float>, combos, px, py, pz, m);
}, {"combos", "px", "py", "pz", "m"});
Benefits
No code duplication: Wrap any single-combo function.
Extensible: Any return type, any argument list.
Performance: As efficient as a hand-written loop.
Summary
Use RVec<config::RVecIndexMap> for all event combinations.
Use ApplyCombinationsGeneric() to vectorize any single-combo kinematic function.
Analysis code stays clean, extensible, and efficient.
Code



---

**To copy:**  
Highlight all text inside the box above, copy, and paste into your documentation, README, or wherever you need it.  
If you need it in a file in your repo, let me know!



