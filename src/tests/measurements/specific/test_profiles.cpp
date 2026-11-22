#include "test_utils.hh"
#include "math_classes.hh"

using namespace penred::measurements;
using penred::grids::RegularGrid;
using penred::grids::VariableGrid;

TEST(profile1D_regular){
  results<double, 2> r;
  std::array<unsigned long,2> nb = {2,3};
  std::array<std::vector<double>,2> lim = {{ {0,2}, {0,3} }};
  CHECK_OK(r.init(nb, lim));

  // Fill with data = index + 1
  for(unsigned long i = 0; i < r.totalBins(); ++i) r.data[i] = double(i+1);

  // Profile on dimension 1 -> collapse dimension 0
  results<double, 1, RegularGrid<1>> p;
  CHECK_OK(r.profile1D(1, p));

  CHECK_EQ(p.totalBins(), 3ul);
  // profile bin 0 = r[0,0]+r[1,0] = 1+2 = 3
  // profile bin 1 = r[0,1]+r[1,1] = 3+4 = 7
  // profile bin 2 = r[0,2]+r[1,2] = 5+6 = 11
  CHECK_NEAR(p.data[0],  3.0, 1e-12);
  CHECK_NEAR(p.data[1],  7.0, 1e-12);
  CHECK_NEAR(p.data[2], 11.0, 1e-12);
}

TEST(profile1D_variable_keeps_edges){
  results<double, 2, VariableGrid<2>> r;
  std::array<unsigned long,2> nb = {2,2};
  std::array<std::vector<double>,2> ed = {{
    {0.0, 1.0, 3.0},
    {0.0, 2.0, 5.0}
  }};
  CHECK_OK(r.init(nb, ed));

  results<double, 1, VariableGrid<1>> p;
  CHECK_OK(r.profile1D(1, p));

  // Profile edges on dim 1 must be {0,2,5}
  const auto& edges = p.readGrid().getEdges();
  CHECK_EQ(edges[0].size(), 3ul);
  CHECK_NEAR(edges[0][0], 0.0, 1e-12);
  CHECK_NEAR(edges[0][1], 2.0, 1e-12);
  CHECK_NEAR(edges[0][2], 5.0, 1e-12);
}

TEST(profile2D_regular){
  results<double, 3> r;
  std::array<unsigned long,3> nb = {2,2,2};
  std::array<std::vector<double>,3> lim = {{ {0,2}, {0,2}, {0,2} }};
  CHECK_OK(r.init(nb, lim));

  for(unsigned long i = 0; i < r.totalBins(); ++i) r.data[i] = 1.0;

  results<double, 2, RegularGrid<2>> p;
  CHECK_OK(r.profile2D(0, 2, p));

  CHECK_EQ(p.totalBins(), 4ul);
  for(unsigned long i = 0; i < p.totalBins(); ++i) CHECK_NEAR(p.data[i], 2.0, 1e-12);
}

TEST(to1D_requires_one_effective_dim){
  results<double, 3> r;
  std::array<unsigned long,3> nb = {1,4,1};
  std::array<std::vector<double>,3> lim = {{ {0,1}, {0,4}, {0,1} }};
  CHECK_OK(r.init(nb, lim));

  results<double, 1, RegularGrid<1>> p;
  CHECK_OK(r.to1D(p));
  CHECK_EQ(p.totalBins(), 4ul);
}

TEST(to2D_requires_two_effective_dims){
  results<double, 3> r;
  std::array<unsigned long,3> nb = {1,2,3};
  std::array<std::vector<double>,3> lim = {{ {0,1}, {0,2}, {0,3} }};
  CHECK_OK(r.init(nb, lim));

  results<double, 2, RegularGrid<2>> p;
  CHECK_OK(r.to2D(p));
  CHECK_EQ(p.totalBins(), 6ul);
}
