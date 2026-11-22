#include "test_utils.hh"
#include "grids.hh"
#include <sstream>

using namespace penred;

TEST(grid_regular_init_and_geometry){
  grids::RegularGrid<2> g;
  std::array<unsigned long, 2> nb = {4, 3};
  std::array<std::vector<double>, 2> lim = {{ {0.0, 4.0}, {-1.0, 2.0} }};
  CHECK_OK(g.init(nb, lim));

  CHECK_EQ(g.nBins[0], 4ul);
  CHECK_EQ(g.nBins[1], 3ul);
  CHECK_EQ(g.totalBins, 12ul);
  CHECK_EQ(g.binsPerIncrement[0], 1ul);
  CHECK_EQ(g.binsPerIncrement[1], 4ul);
  CHECK_NEAR(g.binWidth(0,0), 1.0, 1e-12);
  CHECK_NEAR(g.binWidth(1,0), 1.0, 1e-12);
  CHECK_NEAR(g.binLower(0,2), 2.0, 1e-12);
  CHECK_NEAR(g.binUpper(0,3), 4.0, 1e-12);
  CHECK_NEAR(g.lowEdge(1), -1.0, 1e-12);
  CHECK_NEAR(g.highEdge(1), 2.0, 1e-12);
}

TEST(grid_regular_findBin_clamps){
  grids::RegularGrid<1> g;
  std::array<unsigned long,1> nb = {4};
  std::array<std::vector<double>,1> lim = {{ {0.0, 4.0} }};
  CHECK_OK(g.init(nb, lim));

  CHECK_EQ(g.findBin({-1.0})[0], 0ul);   // below
  CHECK_EQ(g.findBin({ 0.0})[0], 0ul);   // left edge
  CHECK_EQ(g.findBin({ 1.5})[0], 1ul);
  CHECK_EQ(g.findBin({ 3.999})[0], 3ul);
  CHECK_EQ(g.findBin({ 4.0})[0], 3ul);   // right edge clamps
  CHECK_EQ(g.findBin({ 5.0})[0], 3ul);   // above clamps
}

TEST(grid_regular_bracket_frac){
  grids::RegularGrid<1> g;
  std::array<unsigned long,1> nb = {4};
  std::array<std::vector<double>,1> lim = {{ {0.0, 4.0} }};
  CHECK_OK(g.init(nb, lim));

  std::array<unsigned long,1> lo, hi;
  std::array<double,1> frac;

  g.bracket({1.5}, lo, hi, frac);
  CHECK_EQ(lo[0], 1ul); CHECK_EQ(hi[0], 2ul); CHECK_NEAR(frac[0], 0.5, 1e-12);

  // Below low
  g.bracket({-1.0}, lo, hi, frac);
  CHECK_EQ(lo[0], 0ul); CHECK_EQ(hi[0], 0ul); CHECK_NEAR(frac[0], 0.0, 1e-12);

  // Exactly high edge -> last bin, frac 1
  g.bracket({4.0}, lo, hi, frac);
  CHECK_EQ(lo[0], 3ul); CHECK_EQ(hi[0], 3ul); CHECK_NEAR(frac[0], 1.0, 1e-12);

  // Inside last bin
  g.bracket({3.25}, lo, hi, frac);
  CHECK_EQ(lo[0], 3ul); CHECK_EQ(hi[0], 3ul); CHECK_NEAR(frac[0], 0.25, 1e-12);
}

TEST(grid_variable_init_and_find){
  grids::VariableGrid<1> g;
  std::array<unsigned long,1> nb = {3};
  std::array<std::vector<double>,1> ed = {{ {0.0, 0.5, 1.5, 4.0} }};
  CHECK_OK(g.init(nb, ed));

  CHECK_EQ(g.totalBins, 3ul);
  CHECK_NEAR(g.binLower(0,1), 0.5, 1e-12);
  CHECK_NEAR(g.binUpper(0,1), 1.5, 1e-12);
  CHECK_NEAR(g.binWidth(0,1), 1.0, 1e-12);

  CHECK_EQ(g.findBin({-0.1})[0], 0ul);
  CHECK_EQ(g.findBin({ 0.0})[0], 0ul);
  CHECK_EQ(g.findBin({ 0.49})[0], 0ul);
  CHECK_EQ(g.findBin({ 0.5})[0], 1ul);
  CHECK_EQ(g.findBin({ 1.49})[0], 1ul);
  CHECK_EQ(g.findBin({ 1.5})[0], 2ul);
  CHECK_EQ(g.findBin({ 3.99})[0], 2ul);
  CHECK_EQ(g.findBin({ 4.0})[0], 2ul);
  CHECK_EQ(g.findBin({ 9.0})[0], 2ul);
}

TEST(grid_variable_init_rejects_non_monotonic){
  grids::VariableGrid<1> g;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> ed = {{ {0.0, 1.0, 0.5} }};  // not increasing
  CHECK_ERR(g.init(nb, ed), grids::errors::INVALID_LIMITS);
}

TEST(grid_sameAs_regular){
  grids::RegularGrid<1> a, b, c;
  std::array<unsigned long,1> nb = {4};
  std::array<std::vector<double>,1> lim = {{ {0.0, 4.0} }};
  std::array<std::vector<double>,1> lim2 = {{ {0.0, 4.0000004} }};
  std::array<std::vector<double>,1> lim3 = {{ {0.0, 8.0} }};
  CHECK_OK(a.init(nb, lim));
  CHECK_OK(b.init(nb, lim2));
  CHECK_OK(c.init(nb, lim3));
  CHECK_TRUE(a.sameAs(b));    // within relative tolerance
  CHECK_TRUE(!a.sameAs(c));   // different range
}

TEST(grid_sameAs_variable_rejects_opposite_sign){
  grids::VariableGrid<1> a, b;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> e1 = {{ {-1.0, 0.0, 1.0} }};
  std::array<std::vector<double>,1> e2 = {{ {-1.0, 0.0, -1.0} }};  // non-monotonic, will be rejected at init
  CHECK_OK(a.init(nb, e1));
  // e2 invalid, so build b as a valid grid with edges mirroring signs
  std::array<std::vector<double>,1> e3 = {{ {-2.0, -1.0, -0.5} }};
  CHECK_OK(b.init(nb, e3));
  CHECK_TRUE(!a.sameAs(b));
}

TEST(grid_sameAs_variable_sign){
  grids::VariableGrid<1> a, b;
  std::array<unsigned long,1> nb = {1};   // single bin
  std::array<std::vector<double>,1> e1 = {{ {-1.0, 1.0} }};
  std::array<std::vector<double>,1> e2 = {{ { 1.0, -1.0} }};  // invalid (non-monotonic)
  CHECK_OK(a.init(nb, e1));
  // Can't init b with e2, so instead build two 1-bin grids that legitimately
  // differ only in sign of their (single) representative edge:
  std::array<std::vector<double>,1> e3 = {{ {-0.5, 1.0} }};
  std::array<std::vector<double>,1> e4 = {{ { 0.5, 1.0} }};
  CHECK_OK(b.init(nb, e3));
  grids::VariableGrid<1> c;
  CHECK_OK(c.init(nb, e4));
  CHECK_TRUE(!b.sameAs(c));  // -0.5 vs +0.5, opposite signs, both > threshold
}
