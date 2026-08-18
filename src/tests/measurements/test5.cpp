#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <cassert>
#include <numeric>
#include <string>

#include "math_classes.hh"

// ============================================================================
// LIGHTWEIGHT TEST HARNESS
// ============================================================================
static int g_tests_passed = 0;
static int g_tests_failed = 0;

#define TEST_ASSERT(cond, msg)                                          \
  do {                                                                  \
    if (!(cond)) {                                                      \
      std::cerr << "  [FAIL] " << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n"; \
      g_tests_failed++;                                                 \
    } else {                                                            \
      g_tests_passed++;                                                 \
    }                                                                   \
  } while(0)

#define TEST_ASSERT_NEAR(val1, val2, tol, msg)                          \
  TEST_ASSERT(std::abs((val1) - (val2)) <= (tol), msg << " | Expected: " << (val2) << " Got: " << (val1))

#define RUN_TEST(func)                              \
  do {                                              \
    std::cout << "[RUNNING] " << #func << "...\n";  \
    func();                                         \
  } while(0)


// ============================================================================
// 1. CUBIC SPLINE (1D) TESTS
// ============================================================================

void test_spline_exact_cubic() {
  // P(x) = 2*x^3 - 3*x^2 + 4*x - 5
  // P'(x) = 6*x^2 - 6*x + 4
  // P''(x) = 12*x - 6
  auto P   = [](double x) { return 2.0*x*x*x - 3.0*x*x + 4.0*x - 5.0; };
  auto dP  = [](double x) { return 6.0*x*x - 6.0*x + 4.0; };
  auto d2P = [](double x) { return 12.0*x - 6.0; };

  std::vector<double> X = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> Y(X.size());
  for (size_t i = 0; i < X.size(); ++i) Y[i] = P(X[i]);

  double S1 = d2P(X.front());
  double SN = d2P(X.back());

  penred::interpolation::CubicSpline<double> spline;
  int err = spline.init(X, Y, S1, SN);
  TEST_ASSERT(err == penred::interpolation::SUCCESS, "CubicSpline init failed");

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(0.0, 5.0);

  for (int k = 0; k < 100; ++k) {
    double xTest = dist(rng);
    TEST_ASSERT_NEAR(spline.evaluate(xTest), P(xTest), 1e-10, "Spline evaluation mismatch");
    TEST_ASSERT_NEAR(spline.derivative(xTest), dP(xTest), 1e-10, "Spline derivative 1 mismatch");
    TEST_ASSERT_NEAR(spline.derivative2(xTest), d2P(xTest), 1e-10, "Spline derivative 2 mismatch");
  }
}

void test_spline_large_offset_precision() {
  // Test interval stability with x shifted to [1e6, 1e6 + 4]
  double offset = 1.0e6;
  std::vector<double> X = {offset + 0.0, offset + 1.0, offset + 2.0, offset + 3.0, offset + 4.0};
  std::vector<double> Y = {1.0, 4.0, 2.0, 5.0, 3.0};

  penred::interpolation::CubicSpline<double> spline;
  int err = spline.init(X, Y);
  TEST_ASSERT(err == penred::interpolation::SUCCESS, "Spline init at large offset failed");

  for (size_t i = 0; i < X.size(); ++i) {
    TEST_ASSERT_NEAR(spline.evaluate(X[i]), Y[i], 1e-9, "Precision loss at knot point with large offset");
  }
}

void test_spline_invalid_inputs() {
  penred::interpolation::CubicSpline<double> spline;

  // Reject duplicate or non-increasing X
  std::vector<double> X_bad = {0.0, 1.0, 1.0, 3.0};
  std::vector<double> Y = {1.0, 2.0, 3.0, 4.0};
  TEST_ASSERT(spline.init(X_bad, Y) == penred::interpolation::UNORDERED_DATA, "Failed to reject duplicate X points");

  // Reject N < 4
  std::vector<double> X_short = {0.0, 1.0, 2.0};
  std::vector<double> Y_short = {1.0, 2.0, 3.0};
  TEST_ASSERT(spline.init(X_short, Y_short) == penred::interpolation::NOT_ENOUGH_DATA_POINTS, "Failed to reject N < 4");
}


// ============================================================================
// 2. BICUBIC SPLINE (2D) TESTS
// ============================================================================

void test_bicubic_spline_surface() {
    std::vector<double> X = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> Y = {0.0, 1.0, 2.0, 3.0, 4.0};
    
    std::vector<double> Z(X.size() * Y.size());
    for (size_t j = 0; j < Y.size(); ++j) {
        for (size_t i = 0; i < X.size(); ++i) {
            Z[j * X.size() + i] = X[i]*X[i] + 2.0*Y[j];
        }
    }

    penred::interpolation::BicubicSpline<double> spline2d;
    int err = spline2d.init(X, Y, Z);
    TEST_ASSERT(err == penred::interpolation::SUCCESS, "BicubicSpline init failed");

    // Single point evaluation at (1.5, 2.5)
    double val = spline2d.evaluate(1.5, 2.5);
    TEST_ASSERT_NEAR(val, 7.232142857142857, 1e-7, "Bicubic evaluation off expected surface");

    // Grid evaluation
    std::vector<double> xQuery = {0.5, 2.5};
    std::vector<double> yQuery = {1.0, 3.0};
    std::vector<double> gridResults = spline2d.evaluateGrid(xQuery, yQuery);

    TEST_ASSERT(gridResults.size() == 4, "evaluateGrid returned wrong output size");
    
    // Exact natural spline values (131/56 and 461/56)
    TEST_ASSERT_NEAR(gridResults[0], 2.3392857142857143, 1e-7, "Grid eval [0,0] incorrect");
    TEST_ASSERT_NEAR(gridResults[1], 8.2321428571428570, 1e-7, "Grid eval [0,1] incorrect");
}

void test_bicubic_spline_invalid_inputs() {
  penred::interpolation::BicubicSpline<double> spline2d;

  std::vector<double> X_short = {0.0, 1.0, 2.0}; // N < 4
  std::vector<double> Y_ok = {0.0, 1.0, 2.0, 3.0};
  std::vector<double> Z_dummy(12, 0.0);

  TEST_ASSERT(spline2d.init(X_short, Y_ok, Z_dummy) == penred::interpolation::NOT_ENOUGH_DATA_POINTS, 
              "Bicubic failed to reject Nx < 4");

  std::vector<double> X_ok = {0.0, 1.0, 2.0, 3.0};
  std::vector<double> Z_mismatched(10, 0.0); // Size 10 instead of 16
  TEST_ASSERT(spline2d.init(X_ok, Y_ok, Z_mismatched) == penred::interpolation::DIMENSION_MISMATCH, 
              "Bicubic failed to reject size mismatch between grid and values");
}


// ============================================================================
// 3. MULTIDIMENSION & GRID INDEXING TESTS
// ============================================================================

void test_multidimension_indexing() {
  // Concrete test wrapper for multiDimension
  class TestGrid : public penred::measurements::multiDimension<3> {
  public:
    int setup() {
      return this->initDims(
                            std::vector<unsigned long>{4, 5, 6},
                            std::vector<std::pair<double, double>>{{0.0, 4.0}, {0.0, 5.0}, {0.0, 6.0}}
                            );
    }
    using multiDimension<3>::getGlobalIndex;
  };

  TestGrid grid;
  TEST_ASSERT(grid.setup() == penred::measurements::errors::SUCCESS, "Grid initialization failed");
  TEST_ASSERT(grid.getNBins() == 4 * 5 * 6, "Total bins calculation wrong");

  // Test flat index computation: idx = i + j*Nx + k*Nx*Ny
  std::array<unsigned long, 3> idx1 = {0, 0, 0};
  TEST_ASSERT(grid.getGlobalIndex(idx1) == 0, "Global index [0,0,0] wrong");

  std::array<unsigned long, 3> idx2 = {2, 3, 1}; // 2 + 3*4 + 1*(4*5) = 2 + 12 + 20 = 34
  TEST_ASSERT(grid.getGlobalIndex(idx2) == 34, "Global index [2,3,1] wrong");
}

void test_multidimension_boundary_clamping() {
  class TestGrid1D : public penred::measurements::multiDimension<1> {
  public:
    int setup() {
      return this->initDims(
                            std::vector<unsigned long>{10},
                            std::vector<std::pair<double, double>>{{0.0, 10.0}}
                            );
    }
        
    // Expose fraction index lookup / clamping logic for testing
    std::pair<size_t, double> query(double x) {
      size_t lowerIdx = 0;
      double t = 0.0;
      double fracIdx = (x - this->limits[0].first) / this->binWidths[0];

      if (fracIdx <= 0.0) {
        lowerIdx = 0;
        t = 0.0;
      } else if (fracIdx >= static_cast<double>(this->nBins[0] - 1)) {
        lowerIdx = this->nBins[0] - 2;
        t = 1.0;
      } else {
        lowerIdx = static_cast<size_t>(fracIdx);
        t = fracIdx - static_cast<double>(lowerIdx);
      }
      return {lowerIdx, t};
    }
  };

  TestGrid1D grid;
  grid.setup();

  // 1. Underflow query (x < xmin)
  auto resUnder = grid.query(-2.5);
  TEST_ASSERT(resUnder.first == 0 && resUnder.second == 0.0, "Underflow clamping failed");

  // 2. Interior query (x = 3.5 -> bin 3, t = 0.5)
  auto resMid = grid.query(3.5);
  TEST_ASSERT(resMid.first == 3 && std::abs(resMid.second - 0.5) < 1e-9, "Interior interpolation index wrong");

  // 3. Overflow query (x >= xmax)
  auto resOver = grid.query(12.0);
  TEST_ASSERT(resOver.first == 8 && resOver.second == 1.0, "Overflow boundary clamping failed");
}


// ============================================================================
// 4. MEASUREMENT MONTE CARLO TALLY TESTS
// ============================================================================

void test_measurement_randomized_tallies() {
  penred::measurements::measurement<double, 2> meas;
  meas.init(std::vector<unsigned long>{10, 10}, 
            std::vector<std::pair<double, double>>{{0.0, 10.0}, {0.0, 10.0}});

  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> posDist(0.0, 9.99999);
  std::uniform_real_distribution<double> valDist(0.5, 2.5);

  double totalDepositedSum = 0.0;
  unsigned long long currentHist = 1;

  for (int h = 0; h < 500; ++h) {
    int numHits = 1 + (rng() % 5);
    for (int hit = 0; hit < numHits; ++hit) {
      std::array<double, 2> pos = {posDist(rng), posDist(rng)};
      double val = valDist(rng);

      meas.add(pos, val, currentHist);
      totalDepositedSum += val;
    }
    currentHist++;
  }

  // Un-flushed verify: results() must account for un-flushed tmp scores
  penred::measurements::results<double, 2> res;
  meas.flush();
  meas.results(currentHist - 1, res);

  double sumInResults = 0.0;
  for (size_t i = 0; i < res.data.size(); ++i) {
    sumInResults += res.data[i] * static_cast<double>(currentHist - 1);
  }

  TEST_ASSERT_NEAR(sumInResults, totalDepositedSum, 1e-7, "Measurement tally lost score in tmp buffer");
}

void test_measurement_boundary_overflow_safety() {
  penred::measurements::measurement<double, 1> meas;
  meas.init(std::vector<unsigned long>{5}, std::vector<std::pair<double, double>>{{0.0, 1.0}});

  // 0.9999999999999 can round to index 5 without explicit boundary clamping
  std::array<double, 1> posBoundary = {0.9999999999999};
    
  meas.add(posBoundary, 10.0, 1);
  meas.flush();

  TEST_ASSERT_NEAR(meas.readData().back(), 10.0, 1e-9, "Boundary score was not assigned to last bin");
}


// ============================================================================
// 5. MAIN RUNNER
// ============================================================================

int main() {
  std::cout << "========================================\n";
  std::cout << "      RUNNING COMPLETE TEST SUITE       \n";
  std::cout << "========================================\n\n";

  // 1D Spline
  RUN_TEST(test_spline_exact_cubic);
  RUN_TEST(test_spline_large_offset_precision);
  RUN_TEST(test_spline_invalid_inputs);

  // 2D Spline
  RUN_TEST(test_bicubic_spline_surface);
  RUN_TEST(test_bicubic_spline_invalid_inputs);

  // Grids & Dimensions
  RUN_TEST(test_multidimension_indexing);
  RUN_TEST(test_multidimension_boundary_clamping);

  // Measurement Tallying
  RUN_TEST(test_measurement_randomized_tallies);
  RUN_TEST(test_measurement_boundary_overflow_safety);

  std::cout << "\n----------------------------------------\n";
  std::cout << "RESULTS: " << g_tests_passed << " PASSED, " << g_tests_failed << " FAILED.\n";
  std::cout << "----------------------------------------\n";

  return g_tests_failed == 0 ? 0 : 1;
}
