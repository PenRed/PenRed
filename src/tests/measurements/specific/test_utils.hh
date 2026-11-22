#ifndef PENRED_TEST_UTILS_
#define PENRED_TEST_UTILS_

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

namespace penred { namespace test {

  // Counters (one translation unit defines them; see test_main.cpp)
  extern int g_checks;
  extern int g_failures;
  extern const char* g_currentTest;

  inline void setTest(const char* name){ g_currentTest = name; }

  inline void reportFailure(const char* file, int line,
                            const std::string& msg){
    ++g_failures;
    std::fprintf(stderr, "  FAIL [%s] %s:%d: %s\n",
                 g_currentTest ? g_currentTest : "?",
                 file, line, msg.c_str());
  }

  // Basic boolean check
  inline void checkTrue(bool cond, const char* expr,
                        const char* file, int line){
    ++g_checks;
    if(!cond) reportFailure(file, line, std::string("expected true: ") + expr);
  }

  // Approximate equality for doubles
  inline void checkNear(double a, double b, double tol,
                        const char* exprA, const char* exprB,
                        const char* file, int line){
    ++g_checks;
    const double d = std::fabs(a - b);
    const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    if(d > tol * scale){
      char buf[256];
      std::snprintf(buf, sizeof(buf),
                    "expected %s ≈ %s (got %.17g vs %.17g, |Δ|=%.3g > %.3g)",
                    exprA, exprB, a, b, d, tol*scale);
      reportFailure(file, line, buf);
    }
  }

  // Exact equality (for ints, enums, sizes)
  template<class A, class B>
  inline void checkEq(const A& a, const B& b,
                      const char* exprA, const char* exprB,
                      const char* file, int line){
    ++g_checks;
    if(!(a == b)){
      char buf[256];
      std::snprintf(buf, sizeof(buf),
                    "expected %s == %s (got %lld vs %lld)",
                    exprA, exprB,
                    static_cast<long long>(a), static_cast<long long>(b));
      reportFailure(file, line, buf);
    }
  }

    // Check an error code equals SUCCESS (or a specific code)
    template<class A>
    inline void checkOk(const A& actual, const char* expr,
                        const char* file, int line){
      ++g_checks;
      if(static_cast<long long>(actual) != 0){
        char buf[128];
        std::snprintf(buf, sizeof(buf),
                      "%s returned non-zero error (%lld)", expr,
                      static_cast<long long>(actual));
        reportFailure(file, line, buf);
      }
    }

    template<class A, class E>
    inline void checkErr(const A& actual, const E& expected, const char* expr,
                         const char* file, int line){
      ++g_checks;
      if(static_cast<long long>(actual) != static_cast<long long>(expected)){
        char buf[128];
        std::snprintf(buf, sizeof(buf),
                      "%s: expected error %lld, got %lld",
                      expr,
                      static_cast<long long>(expected),
                      static_cast<long long>(actual));
        reportFailure(file, line, buf);
      }
    }

}} // namespace penred::test

#define CHECK_TRUE(x)            ::penred::test::checkTrue((x), #x, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, tol)    ::penred::test::checkNear((a), (b), (tol), #a, #b, __FILE__, __LINE__)
#define CHECK_EQ(a, b)           ::penred::test::checkEq((a), (b), #a, #b, __FILE__, __LINE__)
#define CHECK_OK(x)              ::penred::test::checkOk((x), #x, __FILE__, __LINE__)
#define CHECK_ERR(x, e)          ::penred::test::checkErr((x), (e), #x, __FILE__, __LINE__)

#define TEST(name) \
  static void name(); \
  struct name##_registrar { \
    name##_registrar(){ ::penred::test::setTest(#name); name(); } \
  } name##_instance; \
  static void name()

#endif
