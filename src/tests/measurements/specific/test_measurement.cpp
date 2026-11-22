#include "test_utils.hh"
#include "math_classes.hh"
#include <cmath>

using namespace penred::measurements;

TEST(measurement_add_and_results){
  measurement<double, 1> m;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0.0, 2.0} }};
  CHECK_OK(m.init(nb, lim));

  // Two histories, each depositing in bin 0 and bin 1.
  m.add({0.5}, 1.0, 1);
  m.add({1.5}, 2.0, 1);
  m.add({0.5}, 3.0, 2);
  m.add({1.5}, 4.0, 2);
  m.flush();

  results<double, 1> r;
  m.results(2, r);

  // Bin 0: values 1 and 3 -> mean 2
  // Bin 1: values 2 and 4 -> mean 3
  CHECK_NEAR(r.data[0], 2.0, 1e-12);
  CHECK_NEAR(r.data[1], 3.0, 1e-12);
  // Bin 0: sum=4, sum2=1+9=10; mean=2; var = 10/2 - 4 = 1; sigma=sqrt(1/2)
  CHECK_NEAR(r.sigma[0], std::sqrt(0.5), 1e-12);
}

TEST(measurement_add_rejects_out_of_range){
  measurement<double, 1> m;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0.0, 2.0} }};
  CHECK_OK(m.init(nb, lim));

  m.add({-1.0}, 1.0, 1);   // below
  m.add({ 3.0}, 1.0, 1);   // above
  m.add({ 2.0}, 1.0, 1);   // right edge
  m.flush();

  results<double, 1> r;
  m.results(1, r);
  CHECK_NEAR(r.data[0], 0.0, 1e-12);
  CHECK_NEAR(r.data[1], 0.0, 1e-12);
}

TEST(measurement_cumulative){
  measurement<double, 1> m;
  std::array<unsigned long,1> nb = {3};
  std::array<std::vector<double>,1> lim = {{ {0,3} }};
  CHECK_OK(m.init(nb, lim));

  m.add({0.5}, 1.0, 1);
  m.add({1.5}, 2.0, 1);
  m.add({2.5}, 3.0, 1);
  m.flush();

  measurement<double,1> c = m.cumulative();
  // c[0] = m[0]; c[1] = m[0]+m[1]; c[2] = m[0]+m[1]+m[2]
  CHECK_NEAR(c.readData()[0], 1.0, 1e-12);
  CHECK_NEAR(c.readData()[1], 3.0, 1e-12);
  CHECK_NEAR(c.readData()[2], 6.0, 1e-12);
}

TEST(measurement_add_measurement){
  measurement<double,1> a, b;
  std::array<unsigned long,1> nb = {2};
  std::array<std::vector<double>,1> lim = {{ {0,2} }};
  CHECK_OK(a.init(nb, lim));
  CHECK_OK(b.init(nb, lim));

  a.add({0.5}, 1.0, 1); a.flush();
  b.add({0.5}, 2.0, 1); b.flush();

  CHECK_OK(a.add(b));
  results<double,1> r;
  a.results(1, r);
  CHECK_NEAR(r.data[0], 3.0, 1e-12);
}
