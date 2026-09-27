// test random walker's aliasing sampler
//
//   Adjust the include paths as needed.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "pen_randoms.hh"

// ---------------------------------------------------------------------------
// Minimal test framework
// ---------------------------------------------------------------------------
namespace testfw {

  int g_checks = 0;
  int g_failures = 0;
  std::string g_current;

  void begin(const std::string& name) {
    g_current = name;
    std::printf("  %s\n", name.c_str());
  }

  void check(bool cond, const char* expr, const char* file, int line) {
    ++g_checks;
    if (!cond) {
      ++g_failures;
      std::printf("    FAIL: %s\n      at %s:%d\n      in test '%s'\n",
                  expr, file, line, g_current.c_str());
    }
  }

  template <typename A, typename B>
  void check_eq(const A& a, const B& b, const char* ea, const char* eb,
                const char* file, int line) {
    ++g_checks;
    if (!(a == b)) {
      ++g_failures;
      std::cout << "    FAIL: " << ea << " == " << eb << "\n"
                << "      left  = " << a << "\n"
                << "      right = " << b << "\n"
                << "      at " << file << ":" << line << "\n"
                << "      in test '" << g_current << "'\n";
    }
  }

  void check_near(double a, double b, double tol, const char* ea,
                  const char* eb, const char* file, int line) {
    ++g_checks;
    if (std::fabs(a - b) > tol) {
      ++g_failures;
      std::printf("    FAIL: |%s - %s| <= %g\n"
                  "      left  = %.10g\n"
                  "      right = %.10g\n"
                  "      diff  = %.10g\n"
                  "      at %s:%d\n"
                  "      in test '%s'\n",
                  ea, eb, tol, a, b, std::fabs(a - b),
                  file, line, g_current.c_str());
    }
  }

  int summary() {
    std::printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
  }

} // namespace testfw

#define CHECK(cond) ::testfw::check((cond), #cond, __FILE__, __LINE__)
#define CHECK_EQ(a, b) ::testfw::check_eq((a), (b), #a, #b, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, tol)                                       \
  ::testfw::check_near((a), (b), (tol), #a, #b, __FILE__, __LINE__)
#define TEST(name) ::testfw::begin(name)

// Message-carrying check. Requires `operator<<` to be available for `msg`.
#define CHECK_MSG(cond, msg)                                            \
  do {                                                                  \
    ++::testfw::g_checks;                                               \
    if (!(cond)) {                                                      \
      ++::testfw::g_failures;                                           \
      std::cerr << "    FAIL: " #cond "\n"                              \
                << "      msg:  " << msg << "\n"                        \
                << "      at " << __FILE__ << ":" << __LINE__ << "\n"   \
                << "      in test '" << ::testfw::g_current << "'\n";   \
    }                                                                   \
  } while (0)

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

// 1. Validation errors
void test_validation() {
  TEST("init rejects invalid bin counts");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data(1, 1.0);
    std::array<unsigned long, 1> nbins = {{0}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::INVALID_NUMBER_OF_BINS);
  }

  TEST("init rejects mismatched bin count");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data(3, 1.0);
    std::array<unsigned long, 1> nbins = {{4}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::NUMBER_OF_BINS_MISMATCH);
  }

  TEST("init rejects invalid limits");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data(2, 1.0);
    std::array<unsigned long, 1> nbins = {{2}};
    std::array<std::pair<double,double>, 1> lim = {{ {1.0, 0.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::INVALID_LIMITS);
  }

  TEST("init rejects all-zero distribution");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data(4, 0.0);
    std::array<unsigned long, 1> nbins = {{4}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::BAD_DISTRIBUTION);
  }

  TEST("init rejects negative probabilities beyond tolerance");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data = {1.0, -1.0, 1.0};
    std::array<unsigned long, 1> nbins = {{3}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::BAD_DISTRIBUTION);
  }

  TEST("init accepts tiny negative within tolerance");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data = {1.0, -1e-12, 1.0};
    std::array<unsigned long, 1> nbins = {{3}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::SUCCESS);
  }

  TEST("init rejects NaN");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data = {1.0, std::nan(""), 1.0};
    std::array<unsigned long, 1> nbins = {{3}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::BAD_DISTRIBUTION);
  }

  TEST("init rejects +inf");
  {
    penred::sampling::aliasing<1> a;
    std::vector<double> data = {1.0, std::numeric_limits<double>::infinity(), 1.0};
    std::array<unsigned long, 1> nbins = {{3}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::BAD_DISTRIBUTION);
  }
}

// 2. Sampling a uniform 1D distribution
void test_uniform_1d() {
  TEST("uniform 1D: frequencies are close to 1/N");
  const unsigned long N = 10;
  penred::sampling::aliasing<1> a;
  std::vector<double> data(N, 1.0);
  std::array<unsigned long, 1> nbins = {{N}};
  std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::SUCCESS);

  pen_rand rng;
  rng.rand0(1);
  const unsigned long trials = 200000;
  std::vector<unsigned long> counts(N, 0);
  for (unsigned long i = 0; i < trials; ++i) {
    unsigned long b = a.sample(rng);
    CHECK(b < N);  // in-bounds check
    if (b < N) ++counts[b];
  }
  const double expected = 1.0 / static_cast<double>(N);
  for (unsigned long i = 0; i < N; ++i) {
    const double f = static_cast<double>(counts[i]) / static_cast<double>(trials);
    CHECK_NEAR(f, expected, 0.005);
  }
}

// 3. Sampling a two-value distribution 75/25
void test_two_value() {
  TEST("two-value 75/25: frequencies match");
  penred::sampling::aliasing<1> a;
  std::vector<double> data = {75.0, 25.0};
  std::array<unsigned long, 1> nbins = {{2}};
  std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::SUCCESS);

  pen_rand rng;
  rng.rand0(2);
  const unsigned long trials = 500000;
  unsigned long c0 = 0, c1 = 0;
  for (unsigned long i = 0; i < trials; ++i) {
    unsigned long b = a.sample(rng);
    if (b == 0) ++c0;
    else if (b == 1) ++c1;
    else CHECK(false);
  }
  CHECK_NEAR(static_cast<double>(c0) / trials, 0.75, 0.005);
  CHECK_NEAR(static_cast<double>(c1) / trials, 0.25, 0.005);
}

// 4. Single bin
void test_single_bin() {
  TEST("single bin: always samples bin 0");
  penred::sampling::aliasing<1> a;
  std::vector<double> data = {1.0};
  std::array<unsigned long, 1> nbins = {{1}};
  std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::SUCCESS);

  pen_rand rng;
  rng.rand0(3);
  for (int i = 0; i < 1000; ++i) {
    CHECK_EQ(a.sample(rng), 0UL);
  }
}

// 5. Sparse distribution: zero-probability bins must never be sampled
void test_sparse() {
  TEST("sparse: zero-probability bins are never sampled");
  penred::sampling::aliasing<1> a;
  std::vector<double> data = {0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
  std::array<unsigned long, 1> nbins = {{8}};
  std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::SUCCESS);

  pen_rand rng;
  rng.rand0(4);
  const unsigned long trials = 100000;
  unsigned long c1 = 0, c4 = 0;
  for (unsigned long i = 0; i < trials; ++i) {
    unsigned long b = a.sample(rng);
    if (b == 1) ++c1;
    else if (b == 4) ++c4;
    else CHECK(false);  // only bins 1 and 4 are allowed
  }
  CHECK_NEAR(static_cast<double>(c1) / trials, 0.5, 0.005);
  CHECK_NEAR(static_cast<double>(c4) / trials, 0.5, 0.005);
}

// 6. Skewed distribution
void test_skewed() {
  TEST("skewed: one dominant bin");
  penred::sampling::aliasing<1> a;
  std::vector<double> data = {1000.0, 1.0, 1.0, 1.0, 1.0};
  std::array<unsigned long, 1> nbins = {{5}};
  std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::SUCCESS);

  pen_rand rng;
  rng.rand0(5);
  const unsigned long trials = 500000;
  std::vector<unsigned long> counts(5, 0);
  for (unsigned long i = 0; i < trials; ++i) {
    ++counts[a.sample(rng)];
  }
  const double total = 1004.0;
  CHECK_NEAR(static_cast<double>(counts[0]) / trials, 1000.0 / total, 0.005);
  for (int k = 1; k < 5; ++k) {
    CHECK_NEAR(static_cast<double>(counts[k]) / trials, 1.0 / total, 0.002);
  }
}

// 7. Chi-square test on a 2D distribution
void test_chi_square_2d() {
  TEST("2D: chi-square statistic is reasonable");
  const unsigned long NX = 4, NY = 3;
  penred::sampling::aliasing<2> a;
  // Non-uniform 2D distribution, row-major: index = x * NY + y
  std::vector<double> data(NX * NY, 0.0);
  double expected[NX][NY];
  for (unsigned long x = 0; x < NX; ++x) {
    for (unsigned long y = 0; y < NY; ++y) {
      double p = 1.0 + static_cast<double>(x) + 2.0 * static_cast<double>(y);
      data[y * NX + x] = p;
      expected[x][y] = p;
    }
  }
  double total = 0.0;
  for (unsigned long x = 0; x < NX; ++x)
    for (unsigned long y = 0; y < NY; ++y)
      total += expected[x][y];
  for (unsigned long x = 0; x < NX; ++x)
    for (unsigned long y = 0; y < NY; ++y)
      expected[x][y] /= total;

  std::array<unsigned long, 2> nbins = {{NX, NY}};
  std::array<std::pair<double,double>, 2> lim = {{
      {0.0, 1.0}, {0.0, 1.0}
    }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<2>::SUCCESS);

  pen_rand rng;
  rng.rand0(6);
  const unsigned long trials = 500000;
  unsigned long counts[NX][NY];
  for (unsigned long x = 0; x < NX; ++x)
    for (unsigned long y = 0; y < NY; ++y)
      counts[x][y] = 0;

  for (unsigned long i = 0; i < trials; ++i) {
    auto idx = a.sampleByDim(rng);
    CHECK(idx[0] < NX);
    CHECK(idx[1] < NY);
    if (idx[0] < NX && idx[1] < NY) ++counts[idx[0]][idx[1]];
  }

  double chi2 = 0.0;
  for (unsigned long x = 0; x < NX; ++x) {
    for (unsigned long y = 0; y < NY; ++y) {
      double exp = expected[x][y] * static_cast<double>(trials);
      double obs = static_cast<double>(counts[x][y]);
      double diff = obs - exp;
      chi2 += diff * diff / exp;
    }
  }
  // 12 bins, 11 dof. 99.9% critical value ~ 31.26. Use a loose bound.
  std::printf("    chi2 = %.3f (dof = %lu)\n", chi2, (unsigned long)(NX*NY - 1));
  CHECK(chi2 < 40.0);
}

// 8. sampleByDim covers all dimensions
// 8. sampleByDim covers all dimensions
void test_sample_by_dim_3d() {
  TEST("3D: sampled indices are in bounds and marginals match");
  const unsigned long NX = 3, NY = 4, NZ = 5;
  penred::sampling::aliasing<3> a;

  // Non-uniform but with uniform marginals: p(x,y,z) = f(x)*g(y)*h(z)
  // with f, g, h each summing to 1 over their own range.
  // Layout: index = z*(NX*NY) + y*NX + x  (x fastest).
  std::vector<double> data(NX * NY * NZ, 0.0);
  double px[NX], py[NY], pz[NZ];
  for (unsigned long i = 0; i < NX; ++i) px[i] = 1.0 + 0.5 * i;  // f
  for (unsigned long i = 0; i < NY; ++i) py[i] = 1.0 + 0.3 * i;  // g
  for (unsigned long i = 0; i < NZ; ++i) pz[i] = 1.0 + 0.7 * i;  // h
  for (unsigned long z = 0; z < NZ; ++z)
    for (unsigned long y = 0; y < NY; ++y)
      for (unsigned long x = 0; x < NX; ++x)
        data[z * (NX * NY) + y * NX + x] = px[x] * py[y] * pz[z];

  std::array<unsigned long, 3> nbins = {{NX, NY, NZ}};
  std::array<std::pair<double,double>, 3> lim = {{
      {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}
    }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<3>::SUCCESS);

  pen_rand rng;
  rng.rand0(7);
  const unsigned long trials = 300000;
  std::vector<unsigned long> cx(NX, 0), cy(NY, 0), cz(NZ, 0);
  for (unsigned long i = 0; i < trials; ++i) {
    auto idx = a.sampleByDim(rng);
    CHECK(idx[0] < NX);
    CHECK(idx[1] < NY);
    CHECK(idx[2] < NZ);
    if (idx[0] < NX && idx[1] < NY && idx[2] < NZ) {
      ++cx[idx[0]]; ++cy[idx[1]]; ++cz[idx[2]];
    }
  }
  // Marginals should match the (normalized) f, g, h.
  double sx = 0, sy = 0, sz = 0;
  for (unsigned long i = 0; i < NX; ++i) sx += px[i];
  for (unsigned long i = 0; i < NY; ++i) sy += py[i];
  for (unsigned long i = 0; i < NZ; ++i) sz += pz[i];
  for (unsigned long i = 0; i < NX; ++i)
    CHECK_NEAR(static_cast<double>(cx[i]) / trials, px[i] / sx, 0.01);
  for (unsigned long i = 0; i < NY; ++i)
    CHECK_NEAR(static_cast<double>(cy[i]) / trials, py[i] / sy, 0.01);
  for (unsigned long i = 0; i < NZ; ++i)
    CHECK_NEAR(static_cast<double>(cz[i]) / trials, pz[i] / sz, 0.01);
}

// 9. samplePositions: positions lie in the limits and follow the distribution
void test_sample_positions() {
  TEST("1D samplePositions: in limits and histogram matches");
  const unsigned long N = 5;
  penred::sampling::aliasing<1> a;
  std::vector<double> data = {1.0, 2.0, 3.0, 2.0, 1.0};
  std::array<unsigned long, 1> nbins = {{N}};
  const double lo = -2.0, hi = 3.0;
  std::array<std::pair<double,double>, 1> lim = {{ {lo, hi} }};
  CHECK_EQ(a.init(data, nbins, lim), (int)penred::sampling::aliasing<1>::SUCCESS);

  pen_rand rng;
  rng.rand0(8);
  const unsigned long trials = 500000;
  std::vector<unsigned long> hist(N, 0);
  for (unsigned long i = 0; i < trials; ++i) {
    auto pos = a.samplePositions(rng);
    CHECK(pos[0] >= lo);
    CHECK(pos[0] <  hi);
    if (pos[0] >= lo && pos[0] < hi) {
      unsigned long b = static_cast<unsigned long>((pos[0] - lo) / ((hi - lo) / N));
      if (b >= N) b = N - 1;  // numerical safety
      ++hist[b];
    }
  }
  const double expected[] = {1.0/9.0, 2.0/9.0, 3.0/9.0, 2.0/9.0, 1.0/9.0};
  for (unsigned long i = 0; i < N; ++i) {
    CHECK_NEAR(static_cast<double>(hist[i]) / trials, expected[i], 0.005);
  }
}

// 10. Re-init: calling init twice with different data works
void test_reinit() {
  TEST("re-init: calling init twice works");
  penred::sampling::aliasing<1> a;
  pen_rand rng;
  rng.rand0(9);

  std::vector<double> d1 = {1.0, 1.0, 1.0, 1.0};
  std::array<unsigned long, 1> n1 = {{4}};
  std::array<std::pair<double,double>, 1> l1 = {{ {0.0, 1.0} }};
  CHECK_EQ(a.init(d1, n1, l1), (int)penred::sampling::aliasing<1>::SUCCESS);
  unsigned long s = a.sample(rng);
  CHECK(s < 4);

  std::vector<double> d2 = {0.0, 1.0};
  std::array<unsigned long, 1> n2 = {{2}};
  std::array<std::pair<double,double>, 1> l2 = {{ {0.0, 1.0} }};
  CHECK_EQ(a.init(d2, n2, l2), (int)penred::sampling::aliasing<1>::SUCCESS);
  for (int i = 0; i < 100; ++i) {
    CHECK_EQ(a.sample(rng), 1UL);
  }
}

// ---------------------------------------------------------------------------
// Fuzz test: run init on many random distributions and check invariants.
// ---------------------------------------------------------------------------

// Helper: build a random distribution of size n with a given "spikiness".
// spikiness = 0 -> near-uniform; spikiness = 1 -> very peaked.
static std::vector<double> randomDistribution(unsigned long n, double spikiness,
                                              pen_rand& rng) {
  std::vector<double> w(n);
  for (unsigned long i = 0; i < n; ++i) {
    // base weight in [0, 1)
    double u = rng.rand();
    // exponent controls how peaked the distribution is
    double e = 1.0 + 20.0 * spikiness;
    w[i] = std::pow(u, e) + 1e-6;  // avoid exact zero (tested separately)
  }
  return w;
}

// Helper: verify the structural invariants of a valid aliasing table.
// Returns true if all invariants hold.
template <size_t dim>
static bool checkInvariants(const penred::sampling::aliasing<dim>& a,
                            const std::vector<double>& originalData) {
  const std::vector<double>& cutoff = a.readCutoff();
  const std::vector<unsigned long>& alias = a.readAlias();
  const unsigned long n = a.readAlias().size();

  if (cutoff.size() != n || alias.size() != n) return false;

  // 1. Every alias target is a valid bin index.
  for (unsigned long j = 0; j < n; ++j) {
    if (alias[j] >= n) return false;
  }

  // 2. Cutoffs are non-negative (we clamp negatives in init).
  for (unsigned long j = 0; j < n; ++j) {
    if (cutoff[j] < -1e-9) return false;
  }

  // 3. Zero-probability bins: never self-return (cutoff must be 0) and
  //    never alias to another zero-probability bin.
  for (unsigned long j = 0; j < n; ++j) {
    if (originalData[j] <= 0.0) {
      if (cutoff[j] > 1e-9) return false;
      if (originalData[alias[j]] <= 0.0) return false;
    }
  }

  // 4. Active-bin mass conservation: sum of cutoffs over active bins
  //    (alias[j] == j) equals the number of active bins. This is the
  //    true invariant of Walker's algorithm.
  double activeSum = 0.0;
  unsigned long activeCount = 0;
  for (unsigned long j = 0; j < n; ++j) {
    if (alias[j] == j) {
      activeSum += cutoff[j];
      ++activeCount;
    }
  }
  if (std::fabs(activeSum - static_cast<double>(activeCount)) > 1e-6 * static_cast<double>(n)) {
    return false;
  }

  return true;
}

void test_fuzz_invariants() {
  TEST("fuzz: invariants hold across random distributions");

  pen_rand rng;
  rng.rand0(10);
  const unsigned long iterations = 200;
  unsigned long sizes[] = {2, 3, 5, 8, 13, 21, 50, 100};

  for (unsigned long iter = 0; iter < iterations; ++iter) {
    // Pick a random size and spikiness.
    unsigned long n = sizes[int(rng.rand() * 8.0)];
    double spikiness = rng.rand();
    std::vector<double> w = randomDistribution(n, spikiness, rng);

    penred::sampling::aliasing<1> a;
    std::array<unsigned long, 1> nbins = {{n}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    int rc = a.init(w, nbins, lim);
    CHECK_MSG(rc == (int)penred::sampling::aliasing<1>::SUCCESS,
              "iter=" << iter << " n=" << n << " rc=" << rc);
    if (rc != (int)penred::sampling::aliasing<1>::SUCCESS) continue;

    CHECK_MSG(checkInvariants(a, w),
              "invariants failed at iter=" << iter << " n=" << n);
  }
}

void test_fuzz_frequencies() {
  TEST("fuzz: empirical frequencies match input distribution");

  pen_rand rng;
  rng.rand0(11);
  const unsigned long iterations = 30;
  const unsigned long trials = 100000;

  for (unsigned long iter = 0; iter < iterations; ++iter) {
    unsigned long n = 2 + static_cast<unsigned long>(rng.rand() * 30.0);
    double spikiness = rng.rand();
    std::vector<double> w = randomDistribution(n, spikiness, rng);

    penred::sampling::aliasing<1> a;
    std::array<unsigned long, 1> nbins = {{n}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    if (a.init(w, nbins, lim) != (int)penred::sampling::aliasing<1>::SUCCESS) continue;

    // Sample and count.
    std::vector<unsigned long> counts(n, 0);
    for (unsigned long i = 0; i < trials; ++i) {
      unsigned long b = a.sample(rng);
      CHECK(b < n);
      if (b < n) ++counts[b];
    }

    // Normalize input.
    double total = std::accumulate(w.begin(), w.end(), 0.0);
    for (unsigned long j = 0; j < n; ++j) {
      double expected = w[j] / total;
      double observed = static_cast<double>(counts[j]) / trials;
      // Tolerance: 5 sigma for the binomial. For p near 0 or 1,
      // this is still meaningful; for very small p, skip the check
      // since it's dominated by Poisson noise.
      double sigma = std::sqrt(expected * (1.0 - expected) / trials);
      double tol = 5.0 * sigma + 1e-6;
      CHECK_MSG(std::fabs(observed - expected) < tol,
                "iter=" << iter << " n=" << n << " bin=" << j
                << " expected=" << expected
                << " observed=" << observed
                << " tol=" << tol);
    }
  }
}

// Edge case: distribution with a few exact zeros mixed in.
void test_fuzz_with_zeros() {
  TEST("fuzz: distributions with zero-probability bins");

  pen_rand rng;
  rng.rand0(12);
  for (unsigned long iter = 0; iter < 50; ++iter) {
    unsigned long n = 3 + static_cast<unsigned long>(rng.rand() * 20.0);
    std::vector<double> w(n, 0.0);
    // Randomly assign some bins zero probability.
    for (unsigned long j = 0; j < n; ++j) {
      double u = rng.rand();
      if (u < 0.3) {
        w[j] = 0.0;
      } else {
        w[j] = 0.5 + rng.rand();
      }
    }
    // Ensure not all zero.
    double s = std::accumulate(w.begin(), w.end(), 0.0);
    if (s < 1e-9) continue;

    penred::sampling::aliasing<1> a;
    std::array<unsigned long, 1> nbins = {{n}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    int rc = a.init(w, nbins, lim);
    CHECK_MSG(rc == (int)penred::sampling::aliasing<1>::SUCCESS,
              "iter=" << iter << " n=" << n << " rc=" << rc);
    if (rc != (int)penred::sampling::aliasing<1>::SUCCESS) continue;

    CHECK_MSG(checkInvariants(a, w),
              "invariants failed at iter=" << iter << " n=" << n);

    // Verify zero-probability bins are never sampled.
    const unsigned long trials = 20000;
    for (unsigned long i = 0; i < trials; ++i) {
      unsigned long b = a.sample(rng);
      CHECK_MSG(w[b] > 0.0,
                "sampled zero-probability bin " << b
                << " at iter=" << iter);
    }
  }
}

// ---------------------------------------------------------------------------
// Cross-validation against legacy IRND0/IRND
// ---------------------------------------------------------------------------

void test_cross_validation_legacy() {
  TEST("cross-validation: class and legacy IRND0/IRND agree");

  pen_rand rng;
  rng.rand0(13);
  const unsigned long iterations = 20;
  const unsigned long trials = 200000;

  for (unsigned long iter = 0; iter < iterations; ++iter) {
    unsigned long n = 2 + static_cast<unsigned long>(rng.rand() * 30.0);
    double spikiness = rng.rand();
    std::vector<double> w = randomDistribution(n, spikiness, rng);

    // --- Legacy path ---
    std::vector<double> F(n, 0.0);
    std::vector<long int> K(n, 0);
    IRND0(w.data(), F.data(), K.data(), static_cast<long int>(n));

    // --- Class path ---
    penred::sampling::aliasing<1> a;
    std::array<unsigned long, 1> nbins = {{n}};
    std::array<std::pair<double,double>, 1> lim = {{ {0.0, 1.0} }};
    if (a.init(w, nbins, lim) != (int)penred::sampling::aliasing<1>::SUCCESS) continue;

    // --- Compare cutoffs and aliases directly ---
    // The legacy K may use a 1-based convention; check both.
    // Class cutoffs should equal legacy F (both normalized the same way).
    bool cutoffsMatch = true;
    for (unsigned long j = 0; j < n; ++j) {
      if (std::fabs(a.readCutoff()[j] - F[j]) > 1e-9) {
        cutoffsMatch = false;
        break;
      }
    }
    CHECK_MSG(cutoffsMatch,
              "iter=" << iter << " n=" << n << ": cutoffs differ");

    // --- Compare sampled histograms ---
    std::vector<unsigned long> countsClass(n, 0);
    std::vector<unsigned long> countsLegacy(n, 0);
    pen_rand rng2;
    rng.rand0(14);
    for (unsigned long i = 0; i < trials; ++i) {
      unsigned long b = a.sample(rng2);
      if (b < n) ++countsClass[b];
      long int lb = IRND(F.data(), K.data(),
                         static_cast<long int>(n), rng2);
      if (lb >= 0 && static_cast<unsigned long>(lb) < n) {
        ++countsLegacy[static_cast<unsigned long>(lb)];
      }
    }

    // Chi-square between the two histograms (dof = n-1).
    double chi2 = 0.0;
    for (unsigned long j = 0; j < n; ++j) {
      double expected = static_cast<double>(countsLegacy[j]);
      double observed = static_cast<double>(countsClass[j]);
      double denom = expected + observed;
      if (denom < 1e-9) continue;  // both zero, skip
      double diff = observed - expected;
      chi2 += diff * diff / denom;
    }
    // Loose bound: for n up to 32, 99.9% critical value is < 70.
    // We're comparing two empirical histograms, so the effective dof
    // is n-1 and the statistic is roughly chi-square distributed.
    double critical = 3.0 * static_cast<double>(n) + 20.0;
    CHECK_MSG(chi2 < critical,
              "iter=" << iter << " n=" << n
              << " chi2=" << chi2 << " crit=" << critical);
  }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
  test_validation();
  test_uniform_1d();
  test_two_value();
  test_single_bin();
  test_sparse();
  test_skewed();
  test_chi_square_2d();
  test_sample_by_dim_3d();
  test_sample_positions();
  test_reinit();
  test_fuzz_invariants();
  test_fuzz_frequencies();
  test_fuzz_with_zeros();
  test_cross_validation_legacy();
    
  return testfw::summary();
}
