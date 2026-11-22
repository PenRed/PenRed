#include "test_utils.hh"
#include "math_classes.hh"   // or wherever penred::measurements lives
#include <sstream>
#include <cmath>

using namespace penred::measurements;

TEST(results_init_and_global_index){
  results<double, 2> r;
  std::array<unsigned long,2> nb = {3, 2};
  std::array<std::vector<double>,2> lim = {{ {0,3}, {0,2} }};
  CHECK_OK(r.init(nb, lim));

  CHECK_EQ(r.totalBins(), 6ul);
  CHECK_EQ(r.getGlobalIndex({0,0}), 0ul);
  CHECK_EQ(r.getGlobalIndex({1,0}), 1ul);
  CHECK_EQ(r.getGlobalIndex({2,0}), 2ul);
  CHECK_EQ(r.getGlobalIndex({0,1}), 3ul);
  CHECK_EQ(r.getGlobalIndex({2,1}), 5ul);
}

TEST(results_extractValue_regular_1d_linear){
  results<double, 1> r;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0.0, 2.0} }};
  CHECK_OK(r.init(nb, lim));

  // Put values at the bin left edges: v(0)=0, v(1)=10.
  // Piecewise-linear interpolation with edge-attached values gives
  //   x in [0,1): v = 0 + 10*x
  //   x in [1,2): v = 10 (constant, last bin)
  r.data[0] = 0.0;
  r.data[1] = 10.0;

  double v, s;
  CHECK_OK(r.extractValue(std::array<double,1>{0.0}, v, s)); CHECK_NEAR(v,  0.0, 1e-12);
  CHECK_OK(r.extractValue(std::array<double,1>{0.5}, v, s)); CHECK_NEAR(v,  5.0, 1e-12);
  CHECK_OK(r.extractValue(std::array<double,1>{1.0}, v, s)); CHECK_NEAR(v, 10.0, 1e-12);
  CHECK_OK(r.extractValue(std::array<double,1>{1.5}, v, s)); CHECK_NEAR(v, 10.0, 1e-12);   // last bin constant
  CHECK_OK(r.extractValue(std::array<double,1>{2.0}, v, s)); CHECK_NEAR(v, 10.0, 1e-12);   // right edge
  CHECK_OK(r.extractValue(std::array<double,1>{-1.0}, v, s)); CHECK_NEAR(v, 0.0, 1e-12);   // clamped low
  CHECK_OK(r.extractValue(std::array<double,1>{ 3.0}, v, s)); CHECK_NEAR(v, 10.0, 1e-12);  // clamped high
}

TEST(results_extractValue_variable_1d_linear){
  results<double, 1, penred::grids::VariableGrid<1>> r;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> ed = {{ {0.0, 1.0, 4.0} }};
  CHECK_OK(r.init(nb, ed));

  r.data[0] = 0.0;   // attached at x=0
  r.data[1] = 6.0;   // attached at x=1

  double v, s;
  // Interpolate on [0,1]: v = 0 + 6*x
  CHECK_OK(r.extractValue(std::array<double,1>{0.5}, v, s)); CHECK_NEAR(v, 3.0, 1e-12);
  // On [1,4]: last bin, constant 6
  CHECK_OK(r.extractValue(std::array<double,1>{2.5}, v, s)); CHECK_NEAR(v, 6.0, 1e-12);
  CHECK_OK(r.extractValue(std::array<double,1>{4.0}, v, s)); CHECK_NEAR(v, 6.0, 1e-12);
}

TEST(results_extractValue_2d_bilinear){
  // 2x2 grid, values 0, 1, 10, 11 at the four lower-left corners.
  // Bilinear interpolation should hit these at the corners and average
  // them at the centre of the first "cell".
  results<double, 2> r;
  std::array<unsigned long,2> nb = {2,2};
  std::array<std::vector<double>,2> lim = {{ {0,2}, {0,2} }};
  CHECK_OK(r.init(nb, lim));

  // Global indices: (0,0)=0, (1,0)=1, (0,1)=2, (1,1)=3
  r.data[0] = 0.0;
  r.data[1] = 1.0;
  r.data[2] = 10.0;
  r.data[3] = 11.0;

  double v, s;
  CHECK_OK(r.extractValue(std::array<double,2>{0.0,0.0}, v, s)); CHECK_NEAR(v,  0.0, 1e-12);
  CHECK_OK(r.extractValue(std::array<double,2>{1.0,0.0}, v, s)); CHECK_NEAR(v,  1.0, 1e-12);
  CHECK_OK(r.extractValue(std::array<double,2>{0.0,1.0}, v, s)); CHECK_NEAR(v, 10.0, 1e-12);
  CHECK_OK(r.extractValue(std::array<double,2>{1.0,1.0}, v, s)); CHECK_NEAR(v, 11.0, 1e-12);
  // Bilinear at (0.5, 0.5) = average of 0,1,10,11 = 5.5
  CHECK_OK(r.extractValue(std::array<double,2>{0.5,0.5}, v, s)); CHECK_NEAR(v, 5.5, 1e-12);
}

TEST(results_extractValue_sigma_propagation){
  // Two bins with values 0 and 1 and sigma 0 and sqrt(2).
  // At x=0.5 (midpoint of the first bin, frac=0.5), interpolated sigma is
  //   sqrt(0.5² * 0² + 0.5² * (sqrt(2))²) = sqrt(0.5) ≈ 0.7071
  results<double, 1> r;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0.0, 2.0} }};
  CHECK_OK(r.init(nb, lim));
  r.data[0]  = 0.0;  r.sigma[0] = 0.0;
  r.data[1]  = 1.0;  r.sigma[1] = std::sqrt(2.0);

  double v, s;
  CHECK_OK(r.extractValue({0.5}, v, s));
  CHECK_NEAR(s, std::sqrt(0.5), 1e-12);
}

TEST(results_compatibleDims){
  results<double, 2> a, b;
  std::array<unsigned long,2> nb = {3,2};
  std::array<std::vector<double>,2> lim = {{ {0,3}, {0,2} }};
  CHECK_OK(a.init(nb, lim));
  CHECK_OK(b.init(nb, lim));
  CHECK_TRUE(a.compatibleDims(b));

  std::array<std::vector<double>,2> lim2 = {{ {0,6}, {0,2} }};
  results<double, 2> c;
  CHECK_OK(c.init(nb, lim2));
  CHECK_TRUE(!a.compatibleDims(c));   // different range -> different regular grid
}

TEST(results_normalize){
  results<double, 1> r;
  std::array<unsigned long,1> nb = {3};
  std::array<std::vector<double>,1> lim = {{ {0,3} }};
  CHECK_OK(r.init(nb, lim));
  r.data  = {2.0, 3.0, 5.0};
  r.sigma = {0.2, 0.3, 0.5};

  CHECK_OK(r.normalize());
  CHECK_NEAR(r.data[0], 0.2, 1e-12);
  CHECK_NEAR(r.data[1], 0.3, 1e-12);
  CHECK_NEAR(r.data[2], 0.5, 1e-12);
  // Uncertainties: formula has negative correlation; just check positivity
  for(double s : r.sigma) CHECK_TRUE(s >= 0.0);
}

TEST(results_normalize_rejects_negative){
  results<double, 1> r;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0,2} }};
  CHECK_OK(r.init(nb, lim));
  r.data = {1.0, -1.0};
  CHECK_ERR(r.normalize(), errors::NEGATIVE_BIN_VALUE);
}

TEST(results_subtract_hist){
  results<double, 1> a, b, d;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0,2} }};
  CHECK_OK(a.init(nb, lim));
  CHECK_OK(b.init(nb, lim));
  a.data = {1.0, 2.0}; a.sigma = {0.1, 0.2};
  b.data = {0.5, 0.5}; b.sigma = {0.2, 0.1};

  CHECK_OK(a.subtract(d, b));
  CHECK_NEAR(d.data[0],  0.5, 1e-12);
  CHECK_NEAR(d.data[1],  1.5, 1e-12);
  CHECK_NEAR(d.sigma[0], std::sqrt(0.01 + 0.04), 1e-12);
  CHECK_NEAR(d.sigma[1], std::sqrt(0.04 + 0.01), 1e-12);
}

TEST(results_subtract_scalar){
  results<double, 1> a, d;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0,2} }};
  CHECK_OK(a.init(nb, lim));
  a.data = {1.0, 2.0}; a.sigma = {0.3, 0.4};

  CHECK_OK(a.subtract(d, 1.0, 0.5));
  CHECK_NEAR(d.data[0], 0.0, 1e-12);
  CHECK_NEAR(d.data[1], 1.0, 1e-12);
  CHECK_NEAR(d.sigma[0], std::sqrt(0.09 + 0.25), 1e-12);
}

TEST(results_relative_hist){
  results<double, 1> a, b, d;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0,2} }};
  CHECK_OK(a.init(nb, lim));
  CHECK_OK(b.init(nb, lim));
  a.data = {2.0, 4.0}; a.sigma = {0.2, 0.4};
  b.data = {1.0, 2.0}; b.sigma = {0.1, 0.2};

  CHECK_OK(a.relative(d, b));
  CHECK_NEAR(d.data[0], 2.0, 1e-12);
  CHECK_NEAR(d.data[1], 2.0, 1e-12);
}

TEST(results_relative_div_by_zero_sentinel){
  results<double, 1> a, b, d;
  std::array<unsigned long,1> nb = {1};
  std::array<std::vector<double>,1> lim = {{ {0,1} }};
  CHECK_OK(a.init(nb, lim));
  CHECK_OK(b.init(nb, lim));
  a.data = {3.0}; a.sigma = {0.3};
  b.data = {0.0}; b.sigma = {0.0};

  CHECK_OK(a.relative(d, b, 0.0, -1.0));
  CHECK_NEAR(d.data[0], -1.0, 1e-12);   // sentinel
  CHECK_NEAR(d.sigma[0], 0.0, 1e-12);
}
