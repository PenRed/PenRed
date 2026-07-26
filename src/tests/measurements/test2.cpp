#include "math_classes.hh"
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>

// Simple test framework
namespace test_framework {
    static int tests_run = 0;
    static int tests_passed = 0;
    static int tests_failed = 0;
    static std::string current_test_name;
    
    void start_test(const std::string& name) {
        current_test_name = name;
        tests_run++;
        std::cout << "Running: " << name << "... ";
    }
    
    void test_passed() {
        tests_passed++;
        std::cout << "PASSED" << std::endl;
    }
    
    void test_failed(const std::string& message = "") {
        tests_failed++;
        std::cout << "FAILED";
        if (!message.empty()) {
            std::cout << " - " << message;
        }
        std::cout << std::endl;
    }
    
    void print_summary() {
        std::cout << "\n=== Test Summary ===" << std::endl;
        std::cout << "Total tests: " << tests_run << std::endl;
        std::cout << "Passed: " << tests_passed << std::endl;
        std::cout << "Failed: " << tests_failed << std::endl;
    }
}

// Custom assertion macros
#define TEST_ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            test_framework::test_failed(message); \
            return; \
        } \
    } while(0)

#define TEST_ASSERT_EQ(val1, val2, message) \
    TEST_ASSERT((val1) == (val2), message)

#define TEST_ASSERT_NE(val1, val2, message) \
    TEST_ASSERT((val1) != (val2), message)

#define TEST_ASSERT_NEAR(val1, val2, tolerance, message) \
    TEST_ASSERT(std::abs((val1) - (val2)) <= (tolerance), message)

#define TEST_ASSERT_TRUE(condition, message) \
    TEST_ASSERT(condition, message)

#define TEST_ASSERT_FALSE(condition, message) \
    TEST_ASSERT(!(condition), message)

// Helper function to fill test data
template<size_t dim>
void fillResults(penred::measurements::results<double, dim>& res,
                 const std::function<double(const std::array<unsigned long, dim>&)>& valueFunc) {
    std::array<unsigned long, dim> indexes;
    std::fill(indexes.begin(), indexes.end(), 0ul);
    
    for(unsigned long i = 0; i < res.getNBins(); ++i) {
        // Calculate indexes
        unsigned long temp = i;
        for(size_t d = 0; d < dim; ++d) {
            indexes[d] = temp % res.readDimBins()[d];
            temp /= res.readDimBins()[d];
        }
        
        res.data[i] = valueFunc(indexes);
        res.sigma[i] = 0.1;
    }
}

// Test functions
void test_basic_2d_profile() {
    test_framework::start_test("Basic 2D Profile");
    
    using namespace penred::measurements;
    
    // Create a 3D results object
    results<double, 3> res3D;
    std::vector<unsigned long> bins = {3, 4, 5};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 3.0}, {0.0, 4.0}, {0.0, 5.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 3D results");
    
    // Fill with test data where value = sum of indexes
    fillResults<3>(res3D, [](const std::array<unsigned long, 3>& idx) -> double {
        return static_cast<double>(idx[0] + idx[1] + idx[2]);
    });
    
    // Create 2D profile for dimensions 0 and 1
    results<double, 2> profile2D;
    err = res3D.profile2D(0, 1, profile2D);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D call failed");
    
    // Check dimensions
    TEST_ASSERT_EQ(profile2D.readDimBins()[0], 3ul, "Wrong number of bins for dim 0");
    TEST_ASSERT_EQ(profile2D.readDimBins()[1], 4ul, "Wrong number of bins for dim 1");
    
    // Verify the profile values
    bool values_correct = true;
    for(unsigned long i = 0; i < 3 && values_correct; ++i) {
        for(unsigned long j = 0; j < 4 && values_correct; ++j) {
            std::array<unsigned long, 2> profIdx = {i, j};
            unsigned long globIdx = profile2D.getGlobalIndex(profIdx);
            
            double expected = 5.0 * (i + j) + 10.0; // 10 = sum(0..4)
            if (std::abs(profile2D.data[globIdx] - expected) > 1e-10) {
                values_correct = false;
                std::cout << "\n  Mismatch at (" << i << "," << j << "): expected " 
                         << expected << " got " << profile2D.data[globIdx];
            }
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Profile values incorrect");
    
    // Check headers
    TEST_ASSERT_EQ(profile2D.readDimHeader(0), res3D.readDimHeader(0), "Dim 0 header mismatch");
    TEST_ASSERT_EQ(profile2D.readDimHeader(1), res3D.readDimHeader(1), "Dim 1 header mismatch");
    
    test_framework::test_passed();
}

void test_profile_with_bin_limits() {
  test_framework::start_test("Profile with Bin Limits");
    
  using namespace penred::measurements;
    
  results<double, 3> res3D;
  std::vector<unsigned long> bins = {5, 5, 5};
  std::vector<std::pair<double, double>> limits = {
    {0.0, 5.0}, {0.0, 5.0}, {0.0, 5.0}
  };
    
  int err = res3D.init(bins, limits);
  TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 3D results");
    
  // Fill with index product
  fillResults<3>(res3D, [](const std::array<unsigned long, 3>& idx) -> double {
    return static_cast<double>((idx[0] + 1) * (idx[1] + 1) * (idx[2] + 1));
  });
    
  // Profile dimensions 1 and 2, but only bins 1-3 for dimension 1
  std::array<std::pair<unsigned long, unsigned long>, 3> binLimits;
  binLimits[0] = {0, 5};  // Full range for dim 0 (will be summed over)
  binLimits[1] = {1, 4};  // Bins 1-3 for dim 1
  binLimits[2] = {0, 5};  // Full range for dim 2
    
  results<double, 2> profile2D;
  err = res3D.profile2D(1, 2, binLimits, profile2D);  // Fixed: dimensions first, then limits
  TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D call failed");
    
  // Check dimensions
  TEST_ASSERT_EQ(profile2D.readDimBins()[0], 3ul, "Wrong number of bins for dim 0");
  TEST_ASSERT_EQ(profile2D.readDimBins()[1], 5ul, "Wrong number of bins for dim 1");
    
  // Verify a specific bin value
  std::array<unsigned long, 2> profIdx = {1, 3};  // dim1=2 (index 1 in restricted range), dim2=3
  unsigned long globIdx = profile2D.getGlobalIndex(profIdx);
    
  double expected = 0.0;
  for(unsigned long k = 0; k < 5; ++k) {
    expected += (k + 1) * 3.0 * 4.0;  // (k+1) * (2+1) * (3+1) = (k+1)*3*4
  }
  TEST_ASSERT_NEAR(profile2D.data[globIdx], expected, 1e-10, "Wrong profile value");
    
  test_framework::test_passed();
}

void test_repeated_dimensions_error() {
    test_framework::start_test("Repeated Dimensions Error");
    
    using namespace penred::measurements;
    
    results<double, 3> res3D;
    std::vector<unsigned long> bins = {3, 3, 3};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 3.0}, {0.0, 3.0}, {0.0, 3.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 3D results");
    
    results<double, 2> profile2D;
    err = res3D.profile2D(0, 0, profile2D);
    TEST_ASSERT_EQ(err, errors::DIMENSION_REPEATED, "Should return DIMENSION_REPEATED error");
    
    test_framework::test_passed();
}

void test_out_of_range_dimension_error() {
    test_framework::start_test("Out of Range Dimension Error");
    
    using namespace penred::measurements;
    
    results<double, 3> res3D;
    std::vector<unsigned long> bins = {3, 3, 3};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 3.0}, {0.0, 3.0}, {0.0, 3.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 3D results");
    
    results<double, 2> profile2D;
    err = res3D.profile2D(0, 3, profile2D);
    TEST_ASSERT_EQ(err, errors::DIMENSION_OUT_OF_RANGE, "Should return DIMENSION_OUT_OF_RANGE error");
    
    test_framework::test_passed();
}

void test_vector_bin_limits() {
    test_framework::start_test("Vector Bin Limits");
    
    using namespace penred::measurements;
    
    results<double, 4> res4D;
    std::vector<unsigned long> bins = {3, 4, 5, 6};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 3.0}, {0.0, 4.0}, {0.0, 5.0}, {0.0, 6.0}
    };
    
    int err = res4D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 4D results");
    
    // Fill with constant values
    std::fill(res4D.data.begin(), res4D.data.end(), 1.0);
    std::fill(res4D.sigma.begin(), res4D.sigma.end(), 0.1);
    
    // Profile dimensions 1 and 3 with specific limits
    std::vector<std::array<unsigned long, 3>> vecLimits = {
        {1, 1, 3},  // Dimension 1: bins 1-2
        {3, 2, 5}   // Dimension 3: bins 2-4
    };
    
    results<double, 2> profile2D;
    err = res4D.profile2D(1, 3, vecLimits, profile2D);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D call failed");
    
    // Check dimensions
    TEST_ASSERT_EQ(profile2D.readDimBins()[0], 2ul, "Wrong number of bins for dim 0");
    TEST_ASSERT_EQ(profile2D.readDimBins()[1], 3ul, "Wrong number of bins for dim 1");
    
    // Value should be sum over dimensions 0 and 2: 3 * 5 = 15 bins, each with value 1.0
    double expectedValue = 15.0;
    double expectedSigma = 0.1 * std::sqrt(15.0);
    
    bool values_correct = true;
    for(unsigned long i = 0; i < profile2D.getNBins() && values_correct; ++i) {
        if (std::abs(profile2D.data[i] - expectedValue) > 1e-10) {
            values_correct = false;
            std::cout << "\n  Wrong value at bin " << i << ": expected " 
                     << expectedValue << " got " << profile2D.data[i];
        }
        if (std::abs(profile2D.sigma[i] - expectedSigma) > 1e-10) {
            values_correct = false;
            std::cout << "\n  Wrong sigma at bin " << i << ": expected " 
                     << expectedSigma << " got " << profile2D.sigma[i];
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Profile values incorrect");
    
    test_framework::test_passed();
}

void test_limits_propagation() {
    test_framework::start_test("Limits Propagation");
    
    using namespace penred::measurements;
    
    results<double, 3> res3D;
    std::vector<unsigned long> bins = {5, 5, 5};
    std::vector<std::pair<double, double>> limits = {
        {10.0, 20.0}, {30.0, 40.0}, {50.0, 60.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 3D results");
    
    // Profile dimensions 0 and 2 with restricted range
    std::array<std::pair<unsigned long, unsigned long>, 3> binLimits;
    binLimits[0] = {1, 4};  // bins 1-3
    binLimits[1] = {0, 5};  // full range
    binLimits[2] = {2, 5};  // bins 2-4
    
    results<double, 2> profile2D;
    err = res3D.profile2D(0, 2, binLimits, profile2D);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D call failed");
    
    auto profLimits = profile2D.readLimits();
    
    // Dimension 0 (from original dim 0, bins 1-3, width = 2.0)
    double expected_low = 10.0 + 1 * 2.0;   // 12.0
    double expected_high = 10.0 + 4 * 2.0;  // 18.0
    TEST_ASSERT_NEAR(profLimits[0].first, expected_low, 1e-10, "Wrong low limit for dim 0");
    TEST_ASSERT_NEAR(profLimits[0].second, expected_high, 1e-10, "Wrong high limit for dim 0");
    
    // Dimension 1 (from original dim 2, bins 2-4, width = 2.0)
    expected_low = 50.0 + 2 * 2.0;   // 54.0
    expected_high = 50.0 + 5 * 2.0;  // 60.0
    TEST_ASSERT_NEAR(profLimits[1].first, expected_low, 1e-10, "Wrong low limit for dim 1");
    TEST_ASSERT_NEAR(profLimits[1].second, expected_high, 1e-10, "Wrong high limit for dim 1");
    
    test_framework::test_passed();
}

void test_original_data_unmodified() {
    test_framework::start_test("Original Data Unmodified");
    
    using namespace penred::measurements;
    
    results<double, 3> res3D;
    std::vector<unsigned long> bins = {3, 3, 3};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 3.0}, {0.0, 3.0}, {0.0, 3.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 3D results");
    
    // Fill with known values
    for(unsigned long i = 0; i < res3D.getNBins(); ++i) {
        res3D.data[i] = static_cast<double>(i);
        res3D.sigma[i] = static_cast<double>(i) * 0.1;
    }
    
    // Save original data
    std::vector<double> originalData = res3D.data;
    std::vector<double> originalSigma = res3D.sigma;
    
    results<double, 2> profile2D;
    err = res3D.profile2D(1, 2, profile2D);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D call failed");
    
    // Verify original data is unchanged
    bool data_unchanged = (res3D.data == originalData);
    TEST_ASSERT_TRUE(data_unchanged, "Original data was modified");
    
    bool sigma_unchanged = (res3D.sigma == originalSigma);
    TEST_ASSERT_TRUE(sigma_unchanged, "Original sigma was modified");
    
    TEST_ASSERT_EQ(res3D.getNBins(), 27ul, "Original number of bins changed");
    
    test_framework::test_passed();
}

void test_multiple_profiles() {
    test_framework::start_test("Multiple Profiles");
    
    using namespace penred::measurements;
    
    results<double, 4> res4D;
    std::vector<unsigned long> bins = {2, 3, 4, 5};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 2.0}, {0.0, 3.0}, {0.0, 4.0}, {0.0, 5.0}
    };
    
    int err = res4D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 4D results");
    
    // Fill with known pattern
    fillResults<4>(res4D, [](const std::array<unsigned long, 4>& idx) -> double {
        return static_cast<double>(idx[0] + 2*idx[1] + 3*idx[2] + 4*idx[3]);
    });
    
    // Create two different 2D profiles
    results<double, 2> profile01, profile23;
    
    err = res4D.profile2D(0, 1, profile01);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "First profile2D call failed");
    
    err = res4D.profile2D(2, 3, profile23);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Second profile2D call failed");
    
    // Check dimensions
    TEST_ASSERT_EQ(profile01.readDimBins()[0], 2ul, "Wrong bins for profile01 dim 0");
    TEST_ASSERT_EQ(profile01.readDimBins()[1], 3ul, "Wrong bins for profile01 dim 1");
    TEST_ASSERT_EQ(profile23.readDimBins()[0], 4ul, "Wrong bins for profile23 dim 0");
    TEST_ASSERT_EQ(profile23.readDimBins()[1], 5ul, "Wrong bins for profile23 dim 1");
    
    // The sum of all values in both profiles should be equal
    double sum01 = std::accumulate(profile01.data.begin(), profile01.data.end(), 0.0);
    double sum23 = std::accumulate(profile23.data.begin(), profile23.data.end(), 0.0);
    
    double max_sum = std::max(std::abs(sum01), std::abs(sum23));
    double tolerance = (max_sum > 1.0) ? 1e-10 * max_sum : 1e-10;
    
    TEST_ASSERT_NEAR(sum01, sum23, tolerance, "Profile sums don't match");
    
    test_framework::test_passed();
}

void test_single_bin_dimensions() {
    test_framework::start_test("Single Bin Dimensions");
    
    using namespace penred::measurements;
    
    results<double, 4> res4D;
    std::vector<unsigned long> bins = {1, 5, 5, 1};  // Dims 0 and 3 are single bins
    std::vector<std::pair<double, double>> limits = {
        {0.0, 1.0}, {0.0, 5.0}, {0.0, 5.0}, {0.0, 1.0}
    };
    
    int err = res4D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 4D results");
    
    // Fill with values depending only on dims 1 and 2
    fillResults<4>(res4D, [](const std::array<unsigned long, 4>& idx) -> double {
        return static_cast<double>(idx[1] * 10 + idx[2]);
    });
    
    // Profile dimensions 1 and 2
    results<double, 2> profile2D;
    err = res4D.profile2D(1, 2, profile2D);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D call failed");
    
    // The values should be the same as the original for dims 1,2
    bool values_correct = true;
    for(unsigned long i = 0; i < 5 && values_correct; ++i) {
        for(unsigned long j = 0; j < 5 && values_correct; ++j) {
            std::array<unsigned long, 2> profIdx = {i, j};
            unsigned long globIdx = profile2D.getGlobalIndex(profIdx);
            double expected = static_cast<double>(i * 10 + j);
            if (std::abs(profile2D.data[globIdx] - expected) > 1e-10) {
                values_correct = false;
                std::cout << "\n  Mismatch at (" << i << "," << j << "): expected " 
                         << expected << " got " << profile2D.data[globIdx];
            }
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Profile values incorrect for single bin dims");
    
    test_framework::test_passed();
}

void test_different_data_types() {
    test_framework::start_test("Different Data Types");
    
    using namespace penred::measurements;
    
    // Test with integer type
    results<int, 3> resInt;
    std::vector<unsigned long> bins = {3, 3, 3};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 3.0}, {0.0, 3.0}, {0.0, 3.0}
    };
    
    int err = resInt.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize integer results");
    
    // Fill with values 1-27
    for(unsigned long i = 0; i < resInt.getNBins(); ++i) {
        resInt.data[i] = static_cast<int>(i + 1);
        resInt.sigma[i] = 1.0;
    }
    
    results<int, 2> profileInt;
    err = resInt.profile2D(0, 1, profileInt);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Integer profile2D call failed");
    
    // Check that values are integers (sum of 3 bins)
    for(unsigned long i = 0; i < profileInt.getNBins(); ++i) {
        TEST_ASSERT_TRUE(std::abs(profileInt.data[i] - std::round(profileInt.data[i])) < 1e-10,
                        "Integer profile value is not an integer");
    }
    
    // Test with float type
    results<float, 3> resFloat;
    err = resFloat.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize float results");
    
    std::fill(resFloat.data.begin(), resFloat.data.end(), 2.5f);
    std::fill(resFloat.sigma.begin(), resFloat.sigma.end(), 0.2f);
    
    results<float, 2> profileFloat;
    err = resFloat.profile2D(0, 1, profileFloat);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Float profile2D call failed");
    
    // Each bin should sum over dimension 2 (3 bins)
    float expected = 2.5f * 3.0f;
    bool values_correct = true;
    for(unsigned long i = 0; i < profileFloat.getNBins() && values_correct; ++i) {
        if (std::abs(profileFloat.data[i] - expected) > 1e-5f) {
            values_correct = false;
            std::cout << "\n  Wrong float value at bin " << i << ": expected " 
                     << expected << " got " << profileFloat.data[i];
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Float profile values incorrect");
    
    test_framework::test_passed();
}

void test_full_range_profile() {
    test_framework::start_test("Full Range Profile (no limits specified)");
    
    using namespace penred::measurements;
    
    results<double, 3> res3D;
    std::vector<unsigned long> bins = {4, 4, 4};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 4.0}, {0.0, 4.0}, {0.0, 4.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 3D results");
    
    std::fill(res3D.data.begin(), res3D.data.end(), 1.0);
    std::fill(res3D.sigma.begin(), res3D.sigma.end(), 0.1);
    
    // Create profile without specifying limits (should use full range)
    results<double, 2> profile2D;
    err = res3D.profile2D(1, 2, profile2D);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D call failed");
    
    // Should have all bins from dimensions 1 and 2
    TEST_ASSERT_EQ(profile2D.readDimBins()[0], 4ul, "Wrong number of bins");
    TEST_ASSERT_EQ(profile2D.readDimBins()[1], 4ul, "Wrong number of bins");
    
    // Each bin should have sum of 4 values (from dimension 0)
    double expectedValue = 4.0;
    double expectedSigma = 0.1 * std::sqrt(4.0);
    
    bool values_correct = true;
    for(unsigned long i = 0; i < profile2D.getNBins() && values_correct; ++i) {
        if (std::abs(profile2D.data[i] - expectedValue) > 1e-10) {
            values_correct = false;
            std::cout << "\n  Wrong value at bin " << i;
        }
        if (std::abs(profile2D.sigma[i] - expectedSigma) > 1e-10) {
            values_correct = false;
            std::cout << "\n  Wrong sigma at bin " << i;
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Full range profile values incorrect");
    
    test_framework::test_passed();
}

void test_profile_same_dimension() {
    test_framework::start_test("Profile Same Dimension - 1D to 1D");

    using namespace penred::measurements;
    
    // Create a 1D results object
    results<double, 1> res1D;
    std::vector<unsigned long> bins = {10};
    std::vector<std::pair<double, double>> limits = {{0.0, 10.0}};
    
    int err = res1D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 1D results");
    
    // Fill with known values
    for(unsigned long i = 0; i < res1D.getNBins(); ++i) {
        res1D.data[i] = static_cast<double>(i * i);
        res1D.sigma[i] = 0.1 * i;
    }
    
    // Profile full range
    results<double, 1> profile;
    err = res1D.profile1D(0, profile);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile1D (full range) failed");
    TEST_ASSERT_EQ(profile.getNBins(), res1D.getNBins(), "Full range profile has wrong number of bins");
    
    // Check data matches
    bool data_matches = true;
    for(unsigned long i = 0; i < profile.getNBins() && data_matches; ++i) {
        if(std::abs(profile.data[i] - res1D.data[i]) > 1e-10 ||
           std::abs(profile.sigma[i] - res1D.sigma[i]) > 1e-10) {
            data_matches = false;
        }
    }
    TEST_ASSERT_TRUE(data_matches, "Full range profile data doesn't match original");
    
    // Check limits match
    TEST_ASSERT_NEAR(profile.readLimits()[0].first, limits[0].first, 1e-10, "Full range low limit incorrect");
    TEST_ASSERT_NEAR(profile.readLimits()[0].second, limits[0].second, 1e-10, "Full range high limit incorrect");
    
    // Profile with limits
    std::array<std::pair<unsigned long, unsigned long>, 1> binLimits;
    binLimits[0] = {2, 8};  // Bins 2-7
    
    results<double, 1> profile2;
    err = res1D.profile1D(0, binLimits, profile2);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile1D (with limits) failed");
    TEST_ASSERT_EQ(profile2.getNBins(), 6ul, "Limited profile has wrong number of bins");
    
    // Check data matches original subset
    data_matches = true;
    for(unsigned long i = 0; i < profile2.getNBins() && data_matches; ++i) {
        if(std::abs(profile2.data[i] - res1D.data[i + 2]) > 1e-10 ||
           std::abs(profile2.sigma[i] - res1D.sigma[i + 2]) > 1e-10) {
            data_matches = false;
        }
    }
    TEST_ASSERT_TRUE(data_matches, "Limited profile data doesn't match original subset");
    
    // Check limits of limited profile
    TEST_ASSERT_NEAR(profile2.readLimits()[0].first, 2.0, 1e-10, "Limited profile low limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[0].second, 8.0, 1e-10, "Limited profile high limit incorrect");
    
    test_framework::test_passed();
}

void test_profile_same_dimension_2d() {
    test_framework::start_test("Profile Same Dimension - 2D to 2D");

    using namespace penred::measurements;
    
    // Create a 2D results object
    results<double, 2> res2D;
    std::vector<unsigned long> bins = {5, 6};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 5.0}, {0.0, 6.0}
    };
    
    int err = res2D.init(bins, limits);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "Failed to initialize 2D results");
    
    // Fill with known values: value = global index
    for(unsigned long i = 0; i < res2D.getNBins(); ++i) {
        res2D.data[i] = static_cast<double>(i);
        res2D.sigma[i] = 0.1 * i;
    }
    
    // Profile full range
    results<double, 2> profile;
    err = res2D.profile2D(0, 1, profile);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D (full range) failed");
    TEST_ASSERT_EQ(profile.readDimBins()[0], res2D.readDimBins()[0], "Full range dim 0 wrong");
    TEST_ASSERT_EQ(profile.readDimBins()[1], res2D.readDimBins()[1], "Full range dim 1 wrong");
    TEST_ASSERT_EQ(profile.getNBins(), res2D.getNBins(), "Full range total bins wrong");
    
    // Check data matches
    bool data_matches = true;
    for(unsigned long i = 0; i < profile.getNBins() && data_matches; ++i) {
        if(std::abs(profile.data[i] - res2D.data[i]) > 1e-10 ||
           std::abs(profile.sigma[i] - res2D.sigma[i]) > 1e-10) {
            data_matches = false;
        }
    }
    TEST_ASSERT_TRUE(data_matches, "Full range profile data doesn't match original");
    
    // Check limits match
    TEST_ASSERT_NEAR(profile.readLimits()[0].first, limits[0].first, 1e-10, "Full range x low limit incorrect");
    TEST_ASSERT_NEAR(profile.readLimits()[0].second, limits[0].second, 1e-10, "Full range x high limit incorrect");
    TEST_ASSERT_NEAR(profile.readLimits()[1].first, limits[1].first, 1e-10, "Full range y low limit incorrect");
    TEST_ASSERT_NEAR(profile.readLimits()[1].second, limits[1].second, 1e-10, "Full range y high limit incorrect");
    
    // Profile with limits
    std::array<std::pair<unsigned long, unsigned long>, 2> binLimits;
    binLimits[0] = {1, 4};  // Bins 1-3 in x
    binLimits[1] = {2, 5};  // Bins 2-4 in y
    
    results<double, 2> profile2;
    err = res2D.profile2D(0, 1, binLimits, profile2);
    TEST_ASSERT_EQ(err, errors::SUCCESS, "profile2D (with limits) failed");
    TEST_ASSERT_EQ(profile2.readDimBins()[0], 3ul, "Limited profile dim 0 wrong");
    TEST_ASSERT_EQ(profile2.readDimBins()[1], 3ul, "Limited profile dim 1 wrong");
    TEST_ASSERT_EQ(profile2.getNBins(), 9ul, "Limited profile total bins wrong");
    
    // Check specific values match original subset
    bool values_correct = true;
    for(unsigned long j = 0; j < profile2.readDimBins()[1] && values_correct; ++j) {
        for(unsigned long i = 0; i < profile2.readDimBins()[0] && values_correct; ++i) {
            // Calculate original indices
            unsigned long orig_i = i + binLimits[0].first;  // i + 1
            unsigned long orig_j = j + binLimits[1].first;  // j + 2
            
            // Expected value is the global index in the original array
            double expected = static_cast<double>(orig_j * bins[0] + orig_i);
            double expected_sigma = 0.1 * expected;
            
            // Get profile index
            std::array<unsigned long, 2> idx = {i, j};
            unsigned long globIdx = profile2.getGlobalIndex(idx);
            
            if(std::abs(profile2.data[globIdx] - expected) > 1e-10 ||
               std::abs(profile2.sigma[globIdx] - expected_sigma) > 1e-10) {
                values_correct = false;
            }
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Limited profile data doesn't match original subset");
    
    // Check limits of limited profile
    double expected_x_low = limits[0].first + binLimits[0].first * res2D.readBinWidth(0);
    double expected_x_high = limits[0].first + binLimits[0].second * res2D.readBinWidth(0);
    double expected_y_low = limits[1].first + binLimits[1].first * res2D.readBinWidth(1);
    double expected_y_high = limits[1].first + binLimits[1].second * res2D.readBinWidth(1);
    
    TEST_ASSERT_NEAR(profile2.readLimits()[0].first, expected_x_low, 1e-10, "Limited x low limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[0].second, expected_x_high, 1e-10, "Limited x high limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[1].first, expected_y_low, 1e-10, "Limited y low limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[1].second, expected_y_high, 1e-10, "Limited y high limit incorrect");
    
    // Check headers are preserved
    TEST_ASSERT_EQ(profile2.readDimHeader(0), res2D.readDimHeader(0), "Dim 0 header mismatch");
    TEST_ASSERT_EQ(profile2.readDimHeader(1), res2D.readDimHeader(1), "Dim 1 header mismatch");
    TEST_ASSERT_EQ(profile2.readValueHeader(), res2D.readValueHeader(), "Value header mismatch");
    TEST_ASSERT_EQ(profile2.readSigmaHeader(), res2D.readSigmaHeader(), "Sigma header mismatch");
    
    test_framework::test_passed();
}

int main() {
    std::cout << "Running 2D Profile Unit Tests\n";
    std::cout << "==============================\n\n";
    
    // Run all tests
    test_basic_2d_profile();
    test_profile_with_bin_limits();
    test_repeated_dimensions_error();
    test_out_of_range_dimension_error();
    test_vector_bin_limits();
    test_limits_propagation();
    test_original_data_unmodified();
    test_multiple_profiles();
    test_single_bin_dimensions();
    test_different_data_types();
    test_full_range_profile();
    test_profile_same_dimension();
    test_profile_same_dimension_2d();
    
    // Print summary
    test_framework::print_summary();
    
    // Return 0 on success, 1 on failure
    return (test_framework::tests_failed == 0) ? 0 : 1;
}
