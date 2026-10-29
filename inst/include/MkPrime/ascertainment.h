#ifndef MKPRIME_ASCERTAINMENT_H_PUBLIC
#define MKPRIME_ASCERTAINMENT_H_PUBLIC

// Pseudo-data for ascertainment corrections. Score these patterns with a
// Pruner (pruning.h) that shares the model, rates and tree of the real data:
// its sum_site_likelihoods() is then P(constant character), to pass to
// ascertainment_log_correction() (corrections.h). Keeping that Pruner alive
// lets the correction reuse cached CLs across MCMC proposals.
//
// Characters with different missing-data masks need separate pattern sets:
// a masked tip is marginalised, so P(constant) is over the observed tips.

#include <cstddef>
#include <vector>

namespace MkPrime {

// Tip-state matrix (nTip x nPattern, column-major) of the constant patterns:
// column s has every tip in state s, for s = 0..nPattern-1, and -1 (missing)
// at tips where missing[tip - 1] is non-zero. Pass nPattern = k to sum over
// all k states. Under JC with uniform root frequencies every constant
// pattern has the same probability, so nPattern = 1 and a factor of k is
// equivalent and k times cheaper.
inline std::vector<int> constant_patterns(int nTip, int nPattern,
                                          const unsigned char* missing = nullptr) {
  std::vector<int> states(static_cast<std::size_t>(nTip) * nPattern);
  for (int s = 0; s < nPattern; ++s) {
    for (int t = 0; t < nTip; ++t) {
      states[static_cast<std::size_t>(s) * nTip + t] =
        (missing && missing[t]) ? -1 : s;
    }
  }
  return states;
}

}  // namespace MkPrime

#endif  // MKPRIME_ASCERTAINMENT_H_PUBLIC
