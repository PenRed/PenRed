#include "math_classes.hh"
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <numeric>

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

// Helper function for approximate equality check
bool approx_equal(double a, double b, double tolerance = 1e-6) {
    if (std::abs(a - b) <= tolerance) return true;
    
    // For large values, use relative tolerance
    double max_val = std::max(std::abs(a), std::abs(b));
    if (max_val > 1.0) {
        return std::abs(a - b) <= tolerance * max_val;
    }
    return false;
}

// Namespace aliases for clarity
namespace pmeas = penred::measurements;
namespace pinterp = penred::interpolation;

// Test: Profile 1D to 1D
void test_profile_1d_to_1d() {
    test_framework::start_test("Profile 1D to 1D");
    
    // Create a 1D results object
    pmeas::results<double, 1> res1D;
    std::vector<unsigned long> bins = {10};
    std::vector<std::pair<double, double>> limits = {{0.0, 10.0}};
    
    int err = res1D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 1D results");
    
    // Fill with known values
    for(unsigned long i = 0; i < res1D.getNBins(); ++i) {
        res1D.data[i] = static_cast<double>(i * i);
        res1D.sigma[i] = 0.1 * i;
    }
    
    // Profile full range
    pmeas::results<double, 1> profile;
    err = res1D.profile1D(0, profile);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "profile1D (full range) failed");
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
    
    pmeas::results<double, 1> profile2;
    err = res1D.profile1D(0, binLimits, profile2);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "profile1D (with limits) failed");
    TEST_ASSERT_EQ(profile2.getNBins(), 6ul, "Limited profile has wrong number of bins");
    
    // Check data
    data_matches = true;
    for(unsigned long i = 0; i < profile2.getNBins() && data_matches; ++i) {
        if(std::abs(profile2.data[i] - res1D.data[i + 2]) > 1e-10 ||
           std::abs(profile2.sigma[i] - res1D.sigma[i + 2]) > 1e-10) {
            data_matches = false;
        }
    }
    TEST_ASSERT_TRUE(data_matches, "Limited profile data doesn't match original subset");
    
    // Check limits
    TEST_ASSERT_NEAR(profile2.readLimits()[0].first, 2.0, 1e-10, "Limited profile low limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[0].second, 8.0, 1e-10, "Limited profile high limit incorrect");
    
    test_framework::test_passed();
}

// Test: Profile 2D to 2D
void test_profile_2d_to_2d() {
    test_framework::start_test("Profile 2D to 2D");
    
    // Create a 2D results object
    pmeas::results<double, 2> res2D;
    std::vector<unsigned long> bins = {5, 6};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 5.0}, {0.0, 6.0}
    };
    
    int err = res2D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 2D results");
    
    // Fill with known values
    for(unsigned long i = 0; i < res2D.getNBins(); ++i) {
        res2D.data[i] = static_cast<double>(i);
        res2D.sigma[i] = 0.1 * i;
    }
    
    // Profile full range
    pmeas::results<double, 2> profile;
    err = res2D.profile2D(0, 1, profile);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "profile2D (full range) failed");
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
    
    // Profile with limits
    std::array<std::pair<unsigned long, unsigned long>, 2> binLimits;
    binLimits[0] = {1, 4};  // Bins 1-3 in x
    binLimits[1] = {2, 5};  // Bins 2-4 in y
    
    pmeas::results<double, 2> profile2;
    err = res2D.profile2D(0, 1, binLimits, profile2);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "profile2D (with limits) failed");
    TEST_ASSERT_EQ(profile2.readDimBins()[0], 3ul, "Limited profile dim 0 wrong");
    TEST_ASSERT_EQ(profile2.readDimBins()[1], 3ul, "Limited profile dim 1 wrong");
    
    // Check specific values
    bool values_correct = true;
    for(unsigned long j = 0; j < profile2.readDimBins()[1] && values_correct; ++j) {
        for(unsigned long i = 0; i < profile2.readDimBins()[0] && values_correct; ++i) {
            unsigned long orig_i = i + binLimits[0].first;
            unsigned long orig_j = j + binLimits[1].first;
            double expected = static_cast<double>(orig_j * bins[0] + orig_i);
            
            std::array<unsigned long, 2> idx = {i, j};
            unsigned long globIdx = profile2.getGlobalIndex(idx);
            
            if(std::abs(profile2.data[globIdx] - expected) > 1e-10) {
                values_correct = false;
            }
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Limited profile data doesn't match original subset");
    
    // Check limits
    double expected_x_low = limits[0].first + binLimits[0].first * res2D.readBinWidth(0);
    double expected_x_high = limits[0].first + binLimits[0].second * res2D.readBinWidth(0);
    double expected_y_low = limits[1].first + binLimits[1].first * res2D.readBinWidth(1);
    double expected_y_high = limits[1].first + binLimits[1].second * res2D.readBinWidth(1);
    
    TEST_ASSERT_NEAR(profile2.readLimits()[0].first, expected_x_low, 1e-10, "Limited x low limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[0].second, expected_x_high, 1e-10, "Limited x high limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[1].first, expected_y_low, 1e-10, "Limited y low limit incorrect");
    TEST_ASSERT_NEAR(profile2.readLimits()[1].second, expected_y_high, 1e-10, "Limited y high limit incorrect");
    
    test_framework::test_passed();
}

// Test: 1D Cubic Spline - Basic Interpolation
void test_1d_interpolation_basic() {
    test_framework::start_test("1D Cubic Spline - Basic Interpolation");
    
    // Create a 2D results object
    pmeas::results<double, 2> res2D;
    std::vector<unsigned long> bins = {10, 5};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 10.0}, {0.0, 5.0}
    };
    
    int err = res2D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 2D results");
    
    // Fill with f(x,y) = sin(x) * cos(y)
    for(unsigned long i = 0; i < res2D.getNBins(); ++i) {
        std::array<unsigned long, 2> indexes;
        indexes[0] = i % bins[0];
        indexes[1] = i / bins[0];
        
        double x = limits[0].first + static_cast<double>(indexes[0]) * res2D.readBinWidth(0);
        double y = limits[1].first + static_cast<double>(indexes[1]) * res2D.readBinWidth(1);
        res2D.data[i] = std::sin(x) * std::cos(y);
        res2D.sigma[i] = 0.01;
    }
    
    // Create 1D interpolation along dimension 0 (summing over y)
    pinterp::CubicSpline<double> spline;
    err = res2D.interpolate1D(0, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create 1D spline");
    TEST_ASSERT_TRUE(spline.initialized(), "Spline should be initialized");
    
    // Test interpolation - expected: sin(x) * sum(cos(y_j))
    bool values_correct = true;
    for(double x = 0.5; x < 9.5; x += 0.5) {
        double interpolated = spline.evaluate(x);
        
        double expected = std::sin(x);
        double cos_sum = 0.0;
        for(int j = 0; j < 5; ++j) {
            double y = limits[1].first + static_cast<double>(j) * res2D.readBinWidth(1);
            cos_sum += std::cos(y);
        }
        expected *= cos_sum;
        
        if (!approx_equal(interpolated, expected, 0.2)) {
            values_correct = false;
            std::cout << "\n  Mismatch at x=" << x << ": expected " << expected 
                     << " got " << interpolated;
            break;
        }
    }
    TEST_ASSERT_TRUE(values_correct, "1D interpolation values incorrect");
    
    test_framework::test_passed();
}

// Test: 1D Cubic Spline - With Bin Limits
void test_1d_interpolation_with_limits() {
    test_framework::start_test("1D Cubic Spline - With Bin Limits");
    
    // Create a 3D results object
    pmeas::results<double, 3> res3D;
    std::vector<unsigned long> bins = {20, 10, 5};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 20.0}, {0.0, 10.0}, {0.0, 5.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 3D results");
    
    // Fill with f(x,y,z) = x^2 + y + z
    for(unsigned long i = 0; i < res3D.getNBins(); ++i) {
        unsigned long temp = i;
        std::array<unsigned long, 3> indexes;
        for(size_t d = 0; d < 3; ++d) {
            indexes[d] = temp % bins[d];
            temp /= bins[d];
        }
        
        double x = limits[0].first + static_cast<double>(indexes[0]) * res3D.readBinWidth(0);
        double y = limits[1].first + static_cast<double>(indexes[1]) * res3D.readBinWidth(1);
        double z = limits[2].first + static_cast<double>(indexes[2]) * res3D.readBinWidth(2);
        res3D.data[i] = x * x + y + z;
        res3D.sigma[i] = 0.01;
    }
    
    // Profile dimension 1 with limits on dimension 2
    std::array<std::pair<unsigned long, unsigned long>, 3> binLimits;
    binLimits[0] = {0, 20};   // Full range for dim 0 (summed over)
    binLimits[1] = {0, 10};   // Full range for dim 1 (profiled)
    binLimits[2] = {1, 4};    // Bins 1-3 for dim 2 (summed over)
    
    pinterp::CubicSpline<double> spline;
    err = res3D.interpolate1D(1, binLimits, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create 1D spline with limits");
    
    // Calculate expected: sum(x_i^2) + 20*y + sum(z_k over restricted range)
    double sum_x2 = 0.0;
    for(int i = 0; i < 20; ++i) {
        double x = limits[0].first + static_cast<double>(i) * res3D.readBinWidth(0);
        sum_x2 += x * x;
    }
    
    double sum_z = 0.0;
    for(int k = 1; k < 4; ++k) {
        double z = limits[2].first + static_cast<double>(k) * res3D.readBinWidth(2);
        sum_z += z;
    }
    
    // Test interpolation: expected = sum_x2 + 20*y + sum_z
    bool values_correct = true;
    for(double y = 0.5; y < 9.5; y += 1.0) {
        double interpolated = spline.evaluate(y);
        double expected = sum_x2 + 20.0 * y + sum_z;
        if (!approx_equal(interpolated, expected, 0.01 * std::abs(expected))) {
            values_correct = false;
            std::cout << "\n  Mismatch at y=" << y << ": expected " << expected 
                     << " got " << interpolated;
            break;
        }
    }
    TEST_ASSERT_TRUE(values_correct, "1D interpolation with limits values incorrect");
    
    test_framework::test_passed();
}

// Test: 1D Cubic Spline - Vector Bin Limits
void test_1d_interpolation_vector_limits() {
    test_framework::start_test("1D Cubic Spline - Vector Bin Limits");
    
    // Create a 3D results object
    pmeas::results<double, 3> res3D;
    std::vector<unsigned long> bins = {10, 10, 10};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 10.0}, {0.0, 10.0}, {0.0, 10.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 3D results");
    
    // Fill with f(x,y,z) = x + 2*y + 3*z
    for(unsigned long i = 0; i < res3D.getNBins(); ++i) {
        unsigned long temp = i;
        std::array<unsigned long, 3> idx;
        for(size_t d = 0; d < 3; ++d) {
            idx[d] = temp % bins[d];
            temp /= bins[d];
        }
        
        double x = limits[0].first + static_cast<double>(idx[0]) * res3D.readBinWidth(0);
        double y = limits[1].first + static_cast<double>(idx[1]) * res3D.readBinWidth(1);
        double z = limits[2].first + static_cast<double>(idx[2]) * res3D.readBinWidth(2);
        res3D.data[i] = x + 2.0 * y + 3.0 * z;
        res3D.sigma[i] = 0.01;
    }
    
    // Use vector-based limits for dimension 2
    std::vector<std::array<unsigned long, 3>> vecLimits = {
        {2, 2, 8}  // Dimension 2: bins 2-7
    };
    
    pinterp::CubicSpline<double> spline;
    err = res3D.interpolate1D(1, vecLimits, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create 1D spline with vector limits");
    
    // Expected: sum(x_i) + 2*10*y + 3*sum(z_k over range)
    double sum_x = 0.0;
    for(int i = 0; i < 10; ++i) {
        double x = limits[0].first + static_cast<double>(i) * res3D.readBinWidth(0);
        sum_x += x;
    }
    
    double sum_z = 0.0;
    for(int k = 2; k < 8; ++k) {
        double z = limits[2].first + static_cast<double>(k) * res3D.readBinWidth(2);
        sum_z += z;
    }
    
    bool values_correct = true;
    for(double y = 0.5; y < 9.5; y += 1.0) {
        double interpolated = spline.evaluate(y);
        double expected = sum_x + 20.0 * y + 3.0 * sum_z;
        if (!approx_equal(interpolated, expected, 0.01 * std::abs(expected))) {
            values_correct = false;
            std::cout << "\n  Mismatch at y=" << y << ": expected " << expected 
                     << " got " << interpolated;
            break;
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Vector limits interpolation values incorrect");
    
    test_framework::test_passed();
}

// Test: 1D Cubic Spline - Derivatives
void test_1d_spline_derivatives() {
    test_framework::start_test("1D Cubic Spline - Derivatives");
    
    // Create a 1D profile with f(x) = x^2
    pmeas::results<double, 1> res1D;
    std::vector<unsigned long> bins = {100};
    std::vector<std::pair<double, double>> limits = {{-5.0, 5.0}};
    
    int err = res1D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 1D results");
    
    for(unsigned long i = 0; i < res1D.getNBins(); ++i) {
        double x = limits[0].first + static_cast<double>(i) * res1D.readBinWidth(0);
        res1D.data[i] = x * x;
        res1D.sigma[i] = 0.001;
    }
    
    pinterp::CubicSpline<double> spline;
    err = res1D.interpolate1D(0, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create spline");
    
    // Test function value at x=2: should be 4
    double val = spline.evaluate(2.0);
    TEST_ASSERT_TRUE(approx_equal(val, 4.0, 0.01), "Function value at x=2 incorrect");
    
    // Test first derivative at x=2: should be 2x = 4
    double deriv = spline.derivative(2.0);
    TEST_ASSERT_TRUE(approx_equal(deriv, 4.0, 0.1), "First derivative at x=2 incorrect");
    
    // Test second derivative: should be 2 everywhere
    double deriv2 = spline.derivative2(2.0);
    TEST_ASSERT_TRUE(approx_equal(deriv2, 2.0, 0.2), "Second derivative at x=2 incorrect");
    
    // Test at x=-1: f(-1)=1, f'(-1)=-2
    val = spline.evaluate(-1.0);
    TEST_ASSERT_TRUE(approx_equal(val, 1.0, 0.01), "Function value at x=-1 incorrect");
    
    deriv = spline.derivative(-1.0);
    TEST_ASSERT_TRUE(approx_equal(deriv, -2.0, 0.1), "First derivative at x=-1 incorrect");
    
    // Test at x=0: f(0)=0, f'(0)=0
    val = spline.evaluate(0.0);
    TEST_ASSERT_TRUE(approx_equal(val, 0.0, 0.01), "Function value at x=0 incorrect");
    
    deriv = spline.derivative(0.0);
    TEST_ASSERT_TRUE(approx_equal(deriv, 0.0, 0.2), "First derivative at x=0 incorrect");
    
    test_framework::test_passed();
}

// Test: 1D Cubic Spline - Vector Evaluation
void test_interpolation_vector_evaluation() {
    test_framework::start_test("1D Cubic Spline - Vector Evaluation");
    
    // Create a 1D profile with f(x) = x^2
    pmeas::results<double, 1> res1D;
    std::vector<unsigned long> bins = {20};
    std::vector<std::pair<double, double>> limits = {{0.0, 10.0}};
    
    int err = res1D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 1D results");
    
    for(unsigned long i = 0; i < res1D.getNBins(); ++i) {
        double x = limits[0].first + static_cast<double>(i) * res1D.readBinWidth(0);
        res1D.data[i] = x * x;
        res1D.sigma[i] = 0.001;
    }
    
    pinterp::CubicSpline<double> spline;
    err = res1D.interpolate1D(0, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create spline");
    
    // Test vector evaluation
    std::vector<double> xVals = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    std::vector<double> results = spline.evaluate(xVals);
    
    TEST_ASSERT_EQ(results.size(), xVals.size(), "Vector evaluation size incorrect");
    
    bool values_correct = true;
    for(size_t i = 0; i < xVals.size() && values_correct; ++i) {
        double expected = xVals[i] * xVals[i];
        if (!approx_equal(results[i], expected, 0.1)) {
            values_correct = false;
            std::cout << "\n  Vector evaluation mismatch at x=" << xVals[i] 
                     << ": expected " << expected << " got " << results[i];
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Vector evaluation values incorrect");
    
    test_framework::test_passed();
}

// Test: 2D Bicubic Spline - Basic Interpolation
void test_2d_interpolation_basic() {
    test_framework::start_test("2D Bicubic Spline - Basic Interpolation");
    
    // Create a 2D results object
    pmeas::results<double, 2> res2D;
    std::vector<unsigned long> bins = {8, 6};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 8.0}, {0.0, 6.0}
    };
    
    int err = res2D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 2D results");
    
    // Fill with f(x,y) = sin(x/2) * cos(y/2) + 1
    for(unsigned long i = 0; i < res2D.getNBins(); ++i) {
        std::array<unsigned long, 2> indexes;
        indexes[0] = i % bins[0];
        indexes[1] = i / bins[0];
        
        double x = limits[0].first + static_cast<double>(indexes[0]) * res2D.readBinWidth(0);
        double y = limits[1].first + static_cast<double>(indexes[1]) * res2D.readBinWidth(1);
        res2D.data[i] = std::sin(x/2.0) * std::cos(y/2.0) + 1.0;
        res2D.sigma[i] = 0.01;
    }
    
    // Create 2D interpolation
    pinterp::BicubicSpline<double> spline;
    err = res2D.interpolate2D(0, 1, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create 2D spline");
    
    // Test interpolation - values should match the original function
    bool values_correct = true;
    for(double x = 0.5; x < 7.5; x += 1.0) {
        for(double y = 0.5; y < 5.5; y += 1.0) {
            double interpolated = spline.evaluate(x, y);
            double expected = std::sin(x/2.0) * std::cos(y/2.0) + 1.0;
            if (!approx_equal(interpolated, expected, 0.15)) {
                values_correct = false;
                std::cout << "\n  Mismatch at (" << x << "," << y << "): expected " 
                         << expected << " got " << interpolated;
                goto test_done;
            }
        }
    }
    test_done:
    TEST_ASSERT_TRUE(values_correct, "2D interpolation values incorrect");
    
    test_framework::test_passed();
}

// Test: 2D Bicubic Spline - From 3D Data
void test_2d_interpolation_from_3d() {
    test_framework::start_test("2D Bicubic Spline - From 3D Data");
    
    // Create a 3D results object
    pmeas::results<double, 3> res3D;
    std::vector<unsigned long> bins = {10, 8, 4};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 10.0}, {0.0, 8.0}, {0.0, 4.0}
    };
    
    int err = res3D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 3D results");
    
    // Fill with f(x,y,z) = x*y + z^2
    for(unsigned long i = 0; i < res3D.getNBins(); ++i) {
        unsigned long temp = i;
        std::array<unsigned long, 3> indexes;
        for(size_t d = 0; d < 3; ++d) {
            indexes[d] = temp % bins[d];
            temp /= bins[d];
        }
        
        double x = limits[0].first + static_cast<double>(indexes[0]) * res3D.readBinWidth(0);
        double y = limits[1].first + static_cast<double>(indexes[1]) * res3D.readBinWidth(1);
        double z = limits[2].first + static_cast<double>(indexes[2]) * res3D.readBinWidth(2);
        res3D.data[i] = x * y + z * z;
        res3D.sigma[i] = 0.01;
    }
    
    // Create 2D interpolation of dimensions 0 and 1 (summing over dimension 2)
    pinterp::BicubicSpline<double> spline;
    err = res3D.interpolate2D(0, 1, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create 2D spline from 3D");
    
    // Calculate expected: sum over z of (x*y + z^2) = 4*x*y + sum(z^2)
    double sum_z2 = 0.0;
    for(int k = 0; k < 4; ++k) {
        double z = limits[2].first + static_cast<double>(k) * res3D.readBinWidth(2);
        sum_z2 += z * z;
    }
    
    // Test interpolation
    bool values_correct = true;
    for(double x = 0.5; x < 9.5; x += 2.0) {
        for(double y = 0.5; y < 7.5; y += 2.0) {
            double interpolated = spline.evaluate(x, y);
            double expected = 4.0 * x * y + sum_z2;
            if (!approx_equal(interpolated, expected, 0.01 * std::abs(expected))) {
                values_correct = false;
                std::cout << "\n  Mismatch at (" << x << "," << y << "): expected " 
                         << expected << " got " << interpolated;
                goto test_done2;
            }
        }
    }
    test_done2:
    TEST_ASSERT_TRUE(values_correct, "2D interpolation from 3D values incorrect");
    
    test_framework::test_passed();
}

// Test: 2D Bicubic Spline - With Bin Limits
void test_2d_interpolation_with_limits() {
    test_framework::start_test("2D Bicubic Spline - With Bin Limits");
    
    // Create a 4D results object
    pmeas::results<double, 4> res4D;
    std::vector<unsigned long> bins = {6, 8, 10, 4};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 6.0}, {0.0, 8.0}, {0.0, 10.0}, {0.0, 4.0}
    };
    
    int err = res4D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 4D results");
    
    // Fill with f(w,x,y,z) = w*x + y + z^2
    for(unsigned long i = 0; i < res4D.getNBins(); ++i) {
        unsigned long temp = i;
        std::array<unsigned long, 4> idx;
        for(size_t d = 0; d < 4; ++d) {
            idx[d] = temp % bins[d];
            temp /= bins[d];
        }
        
        double w = limits[0].first + static_cast<double>(idx[0]) * res4D.readBinWidth(0);
        double x = limits[1].first + static_cast<double>(idx[1]) * res4D.readBinWidth(1);
        double y = limits[2].first + static_cast<double>(idx[2]) * res4D.readBinWidth(2);
        double z = limits[3].first + static_cast<double>(idx[3]) * res4D.readBinWidth(3);
        res4D.data[i] = w * x + y + z * z;
        res4D.sigma[i] = 0.01;
    }
    
    // Create 2D interpolation with limits
    std::array<std::pair<unsigned long, unsigned long>, 4> binLimits;
    binLimits[0] = {0, 6};   // Full range for dim 0 (summed over: 6 bins)
    binLimits[1] = {0, 8};   // Full range for dim 1 (profiled: 8 bins)
    binLimits[2] = {2, 8};   // Bins 2-7 for dim 2 (profiled: 6 bins)
    binLimits[3] = {1, 3};   // Bins 1-2 for dim 3 (summed over: 2 bins)
    
    pinterp::BicubicSpline<double> spline;
    err = res4D.interpolate2D(1, 2, binLimits, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create 2D spline with limits");
    
    // Calculate expected values
    // Sum over dim 0 (w): 6 bins, and dim 3 (z): 2 bins
    // f(w,x,y,z) = w*x + y + z^2
    // Sum = x * (nz * sum_w) + (nw * nz) * y + nw * sum_z2
    // where nw = 6, nz = 2
    
    double sum_w = 0.0;
    for(int i = 0; i < 6; ++i) {
        double w = limits[0].first + static_cast<double>(i) * res4D.readBinWidth(0);
        sum_w += w;
    }
    
    double sum_z2 = 0.0;
    for(int k = 1; k < 3; ++k) {
        double z = limits[3].first + static_cast<double>(k) * res4D.readBinWidth(3);
        sum_z2 += z * z;
    }
    
    const unsigned long nw = 6;  // Number of bins in dim 0
    const unsigned long nz = 2;  // Number of bins in dim 3
    
    // Test at several points
    bool values_correct = true;
    double test_points[][2] = {{1.0, 3.0}, {4.0, 6.0}, {7.0, 9.0}};
    
    for(int p = 0; p < 3 && values_correct; ++p) {
        double x_val = test_points[p][0];
        double y_val = test_points[p][1];
        
        double interpolated = spline.evaluate(x_val, y_val);
        double expected = x_val * nz * sum_w + nw * nz * y_val + nw * sum_z2;
        
        if (!approx_equal(interpolated, expected, 0.01 * std::abs(expected))) {
            values_correct = false;
            std::cout << "\n  Mismatch at (" << x_val << "," << y_val 
                     << "): expected " << expected << " got " << interpolated;
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Interpolation with limits values incorrect");
    
    test_framework::test_passed();
}

// Test: 2D Grid Evaluation
void test_2d_grid_evaluation() {
    test_framework::start_test("2D Bicubic Spline - Grid Evaluation");
    
    // Create a 2D results object with f(x,y) = x + y
    pmeas::results<double, 2> res2D;
    std::vector<unsigned long> bins = {5, 5};
    std::vector<std::pair<double, double>> limits = {
        {0.0, 5.0}, {0.0, 5.0}
    };
    
    int err = res2D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 2D results");
    
    for(unsigned long j = 0; j < bins[1]; ++j) {
        for(unsigned long i = 0; i < bins[0]; ++i) {
            double x = limits[0].first + static_cast<double>(i) * res2D.readBinWidth(0);
            double y = limits[1].first + static_cast<double>(j) * res2D.readBinWidth(1);
            res2D.data[j * bins[0] + i] = x + y;
            res2D.sigma[j * bins[0] + i] = 0.01;
        }
    }
    
    // Create 2D interpolation
    pinterp::BicubicSpline<double> spline;
    err = res2D.interpolate2D(0, 1, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create 2D spline");
    
    // Test grid evaluation
    std::vector<double> xVals = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> yVals = {1.0, 2.0, 3.0, 4.0};
    
    std::vector<double> gridResults = spline.evaluateGrid(xVals, yVals);
    TEST_ASSERT_EQ(gridResults.size(), xVals.size() * yVals.size(), "Grid result size incorrect");
    
    // Check grid values
    bool values_correct = true;
    for(size_t j = 0; j < yVals.size() && values_correct; ++j) {
        for(size_t i = 0; i < xVals.size() && values_correct; ++i) {
            double val = gridResults[j * xVals.size() + i];
            double expected = xVals[i] + yVals[j];
            if (!approx_equal(val, expected, 0.15)) {
                values_correct = false;
                std::cout << "\n  Grid mismatch at (" << xVals[i] << "," << yVals[j] 
                         << "): expected " << expected << " got " << val;
            }
        }
    }
    TEST_ASSERT_TRUE(values_correct, "Grid evaluation values incorrect");
    
    test_framework::test_passed();
}

// Test: Error Handling
void test_interpolation_error_handling() {
    test_framework::start_test("Interpolation - Error Handling");
    
    // Test 1: Not enough data points for 1D spline
    pmeas::results<double, 2> res2D_small;
    std::vector<unsigned long> bins = {3, 3};
    std::vector<std::pair<double, double>> limits = {{0.0, 3.0}, {0.0, 3.0}};
    
    int err = res2D_small.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize small 2D results");
    
    for(unsigned long i = 0; i < res2D_small.getNBins(); ++i) {
        res2D_small.data[i] = static_cast<double>(i);
    }
    
    pinterp::CubicSpline<double> spline1D;
    err = res2D_small.interpolate1D(0, spline1D);
    TEST_ASSERT_EQ(err, pinterp::NOT_ENOUGH_DATA_POINTS, 
                   "Should fail with not enough data points");
    TEST_ASSERT_FALSE(spline1D.initialized(), "Spline should not be initialized after error");
    
    // Test 2: Dimension mismatch for bicubic spline
    std::vector<double> xVals = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> yVals = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> values(15);  // Wrong size (should be 16)
    
    pinterp::BicubicSpline<double> spline2D;
    err = spline2D.init(xVals, yVals, values);
    TEST_ASSERT_EQ(err, pinterp::DIMENSION_MISMATCH, 
                   "Should fail with dimension mismatch");
    
    // Test 3: Unordered data points for 1D spline
    std::vector<double> xUnordered = {1.0, 3.0, 2.0, 4.0};  // Not sorted
    std::vector<double> yUnordered = {1.0, 8.0, 4.0, 16.0};
    
    pinterp::CubicSpline<double> spline1D_2;
    err = spline1D_2.init(xUnordered, yUnordered);
    TEST_ASSERT_EQ(err, pinterp::UNORDERED_DATA, 
                   "Should fail with unordered data");
    
    test_framework::test_passed();
}

// Test: Clear and Reuse
void test_interpolation_clear_and_reuse() {
    test_framework::start_test("Interpolation - Clear and Reuse");
    
    // Create first dataset: f(x) = x^2
    pmeas::results<double, 1> res1D_1;
    std::vector<unsigned long> bins1 = {10};
    std::vector<std::pair<double, double>> limits1 = {{0.0, 10.0}};
    
    int err = res1D_1.init(bins1, limits1);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize first dataset");
    
    for(unsigned long i = 0; i < res1D_1.getNBins(); ++i) {
        double x = limits1[0].first + static_cast<double>(i) * res1D_1.readBinWidth(0);
        res1D_1.data[i] = x * x;
    }
    
    pinterp::CubicSpline<double> spline;
    
    // First use
    err = res1D_1.interpolate1D(0, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed first interpolation");
    TEST_ASSERT_TRUE(spline.initialized(), "Spline should be initialized");
    
    double val1 = spline.evaluate(5.0);
    TEST_ASSERT_TRUE(approx_equal(val1, 25.0, 0.5), "First evaluation incorrect");
    
    // Create second dataset: f(x) = x^3
    pmeas::results<double, 1> res1D_2;
    std::vector<unsigned long> bins2 = {10};
    std::vector<std::pair<double, double>> limits2 = {{0.0, 10.0}};
    
    err = res1D_2.init(bins2, limits2);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize second dataset");
    
    for(unsigned long i = 0; i < res1D_2.getNBins(); ++i) {
        double x = limits2[0].first + static_cast<double>(i) * res1D_2.readBinWidth(0);
        res1D_2.data[i] = x * x * x;
    }
    
    // Clear and reinitialize
    spline.clear();
    TEST_ASSERT_FALSE(spline.initialized(), "Spline should not be initialized after clear");
    
    err = res1D_2.interpolate1D(0, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed second interpolation");
    
    // Test second dataset
    double val2 = spline.evaluate(3.0);
    TEST_ASSERT_TRUE(approx_equal(val2, 27.0, 1.0), "Second evaluation at x=3 incorrect");
    
    // Verify first evaluation would be different now
    double val1_new = spline.evaluate(5.0);
    TEST_ASSERT_FALSE(approx_equal(val1_new, 25.0, 0.5), 
                      "Should give different result after reinitialization");
    TEST_ASSERT_TRUE(approx_equal(val1_new, 125.0, 1.0), "Reinitialized evaluation at x=5 incorrect");
    
    test_framework::test_passed();
}

// Test: Out of Range Evaluation
void test_interpolation_out_of_range() {
    test_framework::start_test("Interpolation - Out of Range Evaluation");
    
    // Create a 1D profile with f(x) = x^2 on [0, 10]
    pmeas::results<double, 1> res1D;
    std::vector<unsigned long> bins = {20};
    std::vector<std::pair<double, double>> limits = {{0.0, 10.0}};
    
    int err = res1D.init(bins, limits);
    TEST_ASSERT_EQ(err, pmeas::errors::SUCCESS, "Failed to initialize 1D results");
    
    for(unsigned long i = 0; i < res1D.getNBins(); ++i) {
        double x = limits[0].first + static_cast<double>(i) * res1D.readBinWidth(0);
        res1D.data[i] = x * x;
    }
    
    pinterp::CubicSpline<double> spline;
    err = res1D.interpolate1D(0, spline);
    TEST_ASSERT_EQ(err, pinterp::SUCCESS, "Failed to create spline");
    
    // Test evaluation below the range
    double val_below = spline.evaluate(-1.0);
    TEST_ASSERT_TRUE(std::isfinite(val_below), "Evaluation below range should be finite");
    
    // Test evaluation above the range
    double val_above = spline.evaluate(11.0);
    TEST_ASSERT_TRUE(std::isfinite(val_above), "Evaluation above range should be finite");
    
    test_framework::test_passed();
}

int main() {
    std::cout << "Running Interpolation Unit Tests\n";
    std::cout << "================================\n\n";
    
    // Same-dimension profile tests
    test_profile_1d_to_1d();
    test_profile_2d_to_2d();
    
    // 1D Interpolation tests
    test_1d_interpolation_basic();
    test_1d_interpolation_with_limits();
    test_1d_interpolation_vector_limits();
    test_1d_spline_derivatives();
    test_interpolation_vector_evaluation();
    
    // 2D Interpolation tests
    test_2d_interpolation_basic();
    test_2d_interpolation_from_3d();
    test_2d_interpolation_with_limits();
    test_2d_grid_evaluation();
    
    // General interpolation tests
    test_interpolation_error_handling();
    test_interpolation_clear_and_reuse();
    test_interpolation_out_of_range();
    
    // Print summary
    test_framework::print_summary();
    
    return (test_framework::tests_failed == 0) ? 0 : 1;
}
