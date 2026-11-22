#include "test_utils.hh"
#include "math_classes.hh"
#include <sstream>
#include <cstdio>
#include <fstream>

using namespace penred::measurements;

template<class R>
static std::string capturePrint(const R& r, bool onlyEff){
  const char* path = "test_io_tmp.txt";
  FILE* f = std::fopen(path, "w");
  r.print(f, 1, true, true, onlyEff);
  std::fclose(f);
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  std::remove(path);
  return ss.str();
}

TEST(io_roundtrip_regular_1d){
  results<double,1> r;
  std::array<unsigned long,1> nb = {3};
  std::array<std::vector<double>,1> lim = {{ {0.0, 3.0} }};
  CHECK_OK(r.init(nb, lim));
  r.data  = {1.5, 2.5, 3.5};
  r.sigma = {0.1, 0.2, 0.3};
  r.description = "round-trip test";
  r.setDimHeader(0, "x");
  r.setValueHeader("v");
  r.setSigmaHeader("s");

  const std::string text = capturePrint(r, false);
  std::stringstream in(text);

  results<double,1> r2;
  CHECK_OK(r2.read(in));

  CHECK_EQ(r2.totalBins(), 3ul);
  for(unsigned long i = 0; i < 3; ++i){
    CHECK_NEAR(r2.data[i],  r.data[i],  1e-12);
    CHECK_NEAR(r2.sigma[i], r.sigma[i], 1e-12);
  }
}

TEST(io_roundtrip_regular_2d_effective){
  results<double,2> r;
  std::array<unsigned long,2> nb = {1,3};   // dim 0 is a singleton
  std::array<std::vector<double>,2> lim = {{ {0,1}, {0,3} }};
  CHECK_OK(r.init(nb, lim));
  r.data = {1.0, 2.0, 3.0};

  const std::string text = capturePrint(r, true);   // printOnlyEffective
  std::stringstream in(text);

  results<double,2> r2;
  CHECK_OK(r2.read(in));
  CHECK_EQ(r2.totalBins(), 3ul);
  for(unsigned long i = 0; i < 3; ++i) CHECK_NEAR(r2.data[i], r.data[i], 1e-12);
}

TEST(io_roundtrip_variable_1d){
  results<double,1, penred::grids::VariableGrid<1>> r;
  std::array<unsigned long,1> nb = {3};
  std::array<std::vector<double>,1> ed = {{ {0.0, 0.5, 2.0, 4.0} }};
  CHECK_OK(r.init(nb, ed));
  r.data  = {1.0, 2.0, 3.0};
  r.sigma = {0.1, 0.2, 0.3};

  const char* path = "test_io_var_tmp.txt";
  FILE* f = std::fopen(path, "w");
  r.print(f, 1, true, true, false);
  std::fclose(f);

  std::ifstream in(path);
  decltype(r) r2;
  CHECK_OK(r2.read(in));
  std::remove(path);

  const auto& edges = r2.readGrid().getEdges();
  CHECK_EQ(edges[0].size(), 4ul);
  CHECK_NEAR(edges[0][0], 0.0, 1e-12);
  CHECK_NEAR(edges[0][1], 0.5, 1e-12);
  CHECK_NEAR(edges[0][2], 2.0, 1e-12);
  CHECK_NEAR(edges[0][3], 4.0, 1e-12);
  for(unsigned long i = 0; i < 3; ++i) CHECK_NEAR(r2.data[i], r.data[i], 1e-12);
}
