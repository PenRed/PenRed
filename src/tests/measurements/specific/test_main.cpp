#include "test_utils.hh"

namespace penred { namespace test {
  int g_checks = 0;
  int g_failures = 0;
  const char* g_currentTest = nullptr;
}}

int main(){
  using namespace penred::test;
  std::printf("\n%d checks, %d failures\n", g_checks, g_failures);
  return g_failures == 0 ? 0 : 1;
}
