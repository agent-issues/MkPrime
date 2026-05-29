#ifndef MKPRIME_CHAIN_RNG_H
#define MKPRIME_CHAIN_RNG_H

#include <cstdint>
#include <cmath>
#include <limits>
#include <random>

// Per-chain RNG. One instance per MCMC chain, owned by run_mcmc_batch_cpp.
// All distribution parameter orderings match Rmath (R::r*) so call-site
// rewrites are mechanical: R::rgamma(shape, scale) -> rng.rgamma(shape, scale).
//
// Phase 2a: code remains serial — chains are dispatched sequentially.
// Phase 2b: the outer chain loop gains #pragma omp parallel for; ChainRng
//           instances carry no R-API dependencies so they are safe in
//           non-master threads.
struct ChainRng {
  std::mt19937_64 gen;

  explicit ChainRng(uint64_t seed) : gen(seed) {}

  // U(0, 1) -- matches R::unif_rand() / unif_rand()
  double unif() {
    return std::generate_canonical<double,
                                   std::numeric_limits<double>::digits>(gen);
  }

  // Normal(mu, sd) -- matches R::rnorm(mu, sd)
  double rnorm(double mu, double sd) {
    std::normal_distribution<double> d(mu, sd);
    return d(gen);
  }

  // Gamma(shape, scale) -- matches R::rgamma(shape, scale).
  // NOTE: std::gamma_distribution(alpha, beta) uses beta as SCALE
  // (same convention as Rmath), confirmed by cppreference: "beta = scale".
  double rgamma(double shape, double scale) {
    std::gamma_distribution<double> d(shape, scale);
    return d(gen);
  }

  // Beta(a, b) -- matches R::rbeta(a, b). std lacks beta; build via
  // two gammas: X ~ Gamma(a, 1), Y ~ Gamma(b, 1), Beta = X / (X + Y).
  double rbeta(double a, double b) {
    double x = rgamma(a, 1.0);
    double y = rgamma(b, 1.0);
    double s = x + y;
    if (s <= 0.0) return 0.5;  // numerical safety; should not occur for valid a, b > 0
    return x / s;
  }
};

#endif // MKPRIME_CHAIN_RNG_H
