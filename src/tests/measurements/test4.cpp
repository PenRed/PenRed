#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <iomanip>
#include <sstream>

// Include your library here
#include "math_classes.hh"

// Simple test framework without external dependencies
class TestSuite {
private:
    int testsPassed = 0;
    int testsFailed = 0;
    std::string currentTest;
    
    static constexpr double EPSILON = 1.0e-10;
    
public:
    void startTest(const std::string& testName) {
        currentTest = testName;
        std::cout << "Running: " << testName << "... ";
    }
    
    void assertTrue(bool condition, const std::string& message = "") {
        if (condition) {
            pass();
        } else {
            fail(message);
        }
    }
    
    void assertEqual(double expected, double actual, const std::string& message = "") {
        if (std::abs(expected - actual) < EPSILON) {
            pass();
        } else {
            std::stringstream ss;
            ss << message << " (Expected: " << expected 
               << ", Got: " << actual << ")";
            fail(ss.str());
        }
    }
    
    void assertEqual(int expected, int actual, const std::string& message = "") {
        if (expected == actual) {
            pass();
        } else {
            std::stringstream ss;
            ss << message << " (Expected: " << expected 
               << ", Got: " << actual << ")";
            fail(ss.str());
        }
    }
    
    void assertEqual(unsigned long expected, unsigned long actual, const std::string& message = "") {
        if (expected == actual) {
            pass();
        } else {
            std::stringstream ss;
            ss << message << " (Expected: " << expected 
               << ", Got: " << actual << ")";
            fail(ss.str());
        }
    }
    
    void assertEqual(const std::string& expected, const std::string& actual, const std::string& message = "") {
        if (expected == actual) {
            pass();
        } else {
            std::stringstream ss;
            ss << message << " (Expected: \"" << expected 
               << "\", Got: \"" << actual << "\")";
            fail(ss.str());
        }
    }
    
    int getTestsFailed() const { return testsFailed; }
    
    void printSummary() {
        std::cout << "\n" << std::string(50, '=') << std::endl;
        std::cout << "Test Summary:" << std::endl;
        std::cout << "  Passed: " << testsPassed << std::endl;
        std::cout << "  Failed: " << testsFailed << std::endl;
        std::cout << "  Total:  " << (testsPassed + testsFailed) << std::endl;
        std::cout << std::string(50, '=') << std::endl;
    }
    
private:
    void pass() {
        std::cout << "PASSED" << std::endl;
        testsPassed++;
    }
    
    void fail(const std::string& message = "") {
        std::cout << "FAILED";
        if (!message.empty()) {
            std::cout << " - " << message;
        }
        std::cout << std::endl;
        testsFailed++;
    }
};

// Helper functions to create test grids
penred::measurements::results<double, 1> createLinearGrid1D(
    const std::array<unsigned long, 1>& bins,
    const std::array<std::pair<double, double>, 1>& limits) {
    
    penred::measurements::results<double, 1> grid;
    grid.init(bins, limits);
    
    grid.forEach([&](unsigned long idx, const std::array<unsigned long, 1>& localIdx) {
        double x = limits[0].first + (localIdx[0] + 0.5) * grid.readBinWidth(0);
        grid.data[idx] = 2.0 * x;
        grid.sigma[idx] = 0.1;
    });
    
    return grid;
}

penred::measurements::results<double, 2> createLinearGrid2D(
    const std::array<unsigned long, 2>& bins,
    const std::array<std::pair<double, double>, 2>& limits) {
    
    penred::measurements::results<double, 2> grid;
    grid.init(bins, limits);
    
    grid.forEach([&](unsigned long idx, const std::array<unsigned long, 2>& localIdx) {
        double x = limits[0].first + (localIdx[0] + 0.5) * grid.readBinWidth(0);
        double y = limits[1].first + (localIdx[1] + 0.5) * grid.readBinWidth(1);
        grid.data[idx] = 2.0 * x + 3.0 * y;
        grid.sigma[idx] = 0.1;
    });
    
    return grid;
}

penred::measurements::results<double, 3> createLinearGrid3D(
    const std::array<unsigned long, 3>& bins,
    const std::array<std::pair<double, double>, 3>& limits) {
    
    penred::measurements::results<double, 3> grid;
    grid.init(bins, limits);
    
    grid.forEach([&](unsigned long idx, const std::array<unsigned long, 3>& localIdx) {
        double x = limits[0].first + (localIdx[0] + 0.5) * grid.readBinWidth(0);
        double y = limits[1].first + (localIdx[1] + 0.5) * grid.readBinWidth(1);
        double z = limits[2].first + (localIdx[2] + 0.5) * grid.readBinWidth(2);
        grid.data[idx] = 2.0 * x + 3.0 * y + 4.0 * z;
        grid.sigma[idx] = 0.1;
    });
    
    return grid;
}

// =============== extractValue Tests ===============

// Test extractValue with 1D grid and all overloads
void testExtractValue1D(TestSuite& ts) {
    ts.startTest("extractValue - 1D grid");
    
    std::array<unsigned long, 1> bins = {{10}};
    std::array<std::pair<double, double>, 1> limits = {{
        std::pair<double, double>(0.0, 10.0)
    }};
    
    auto grid = createLinearGrid1D(bins, limits);
    
    double value, uncertainty;
    int err;
    
    // Test with std::array overload
    // Position 2.5: between bin 2 (left edge 2.0, value 5.0) and bin 3 (left edge 3.0, value 7.0)
    std::array<double, 1> pos_arr = {{2.5}};
    err = grid.extractValue(pos_arr, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Array overload: SUCCESS");
    ts.assertEqual(6.0, value, "Array overload: f(2.5) = 6.0");
    
    // Test with std::vector overload
    std::vector<double> pos_vec = {2.5};
    err = grid.extractValue(pos_vec, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Vector overload: SUCCESS");
    ts.assertEqual(6.0, value, "Vector overload: f(2.5) = 6.0");
    
    // Test with scalar overload (1D)
    err = grid.extractValue(2.5, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Scalar overload: SUCCESS");
    ts.assertEqual(6.0, value, "Scalar overload: f(2.5) = 6.0");
    
    // Test at bin left edge
    err = grid.extractValue(2.0, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "At left edge: SUCCESS");
    ts.assertEqual(5.0, value, "f(2.0) = 5.0");
    
    // Test out of bounds (clamped to edge)
    err = grid.extractValue(-1.0, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Out of bounds lower: SUCCESS");
    ts.assertEqual(1.0, value, "Clamped to first bin value");
    
    err = grid.extractValue(11.0, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Out of bounds upper: SUCCESS");
    ts.assertEqual(19.0, value, "Clamped to last bin value");
    
    // Test vector dimension mismatch
    std::vector<double> bad_vec = {1.0, 2.0};
    err = grid.extractValue(bad_vec, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::DIMENSION_MISMATCH, err, 
                  "Vector size mismatch should return error");
}

// Test extractValue with 2D grid and all overloads
void testExtractValue2D(TestSuite& ts) {
    ts.startTest("extractValue - 2D grid");
    
    std::array<unsigned long, 2> bins = {{10, 10}};
    std::array<std::pair<double, double>, 2> limits = {{
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0)
    }};
    
    auto grid = createLinearGrid2D(bins, limits);
    
    double value, uncertainty;
    int err;
    
    // Test with std::array overload
    // Position (0.5, 0.5): both on bin left edges (5*0.1=0.5), so no interpolation
    std::array<double, 2> pos_arr = {{0.5, 0.5}};
    err = grid.extractValue(pos_arr, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Array overload: SUCCESS");
    ts.assertEqual(2.75, value, "Array overload: f(0.5,0.5) = 2.75");
    
    // Test with std::vector overload
    std::vector<double> pos_vec = {0.5, 0.5};
    err = grid.extractValue(pos_vec, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Vector overload: SUCCESS");
    ts.assertEqual(2.75, value, "Vector overload: f(0.5,0.5) = 2.75");
    
    // Test with two scalars overload (2D)
    err = grid.extractValue(0.5, 0.5, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Two scalars overload: SUCCESS");
    ts.assertEqual(2.75, value, "Two scalars overload: f(0.5,0.5) = 2.75");
    
    // Test arbitrary point
    // Position (0.2, 0.7): both on bin left edges, no interpolation
    err = grid.extractValue(0.2, 0.7, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Arbitrary point: SUCCESS");
    ts.assertEqual(2.75, value, "f(0.2,0.7) = 2.75");
    
    // Test with interpolation
    // Position (0.25, 0.55): x between edges 0.2 and 0.3, y between edges 0.5 and 0.6
    err = grid.extractValue(0.25, 0.55, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Interpolated point: SUCCESS");
    // Expected: x fraction=0.5, y fraction=0.5
    // Value = average of (2,5), (2,6), (3,5), (3,6)
    // = 0.25*(value at 2,5 + value at 2,6 + value at 3,5 + value at 3,6)
    // = 0.25*(0.2*2+0.3*5+0.25 + 0.2*2+0.3*6+0.25 + 0.2*3+0.3*5+0.25 + 0.2*3+0.3*6+0.25)
    // Let's compute: 0.25*(0.4+1.5+0.25 + 0.4+1.8+0.25 + 0.6+1.5+0.25 + 0.6+1.8+0.25)
    // = 0.25*(2.15 + 2.45 + 2.35 + 2.65) = 0.25*9.6 = 2.4
    ts.assertEqual(2.4, value, "f(0.25,0.55) = 2.4");
}

// Test extractValue with 3D grid and all overloads
void testExtractValue3D(TestSuite& ts) {
    ts.startTest("extractValue - 3D grid");
    
    std::array<unsigned long, 3> bins = {{5, 5, 5}};
    std::array<std::pair<double, double>, 3> limits = {{
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0)
    }};
    
    auto grid = createLinearGrid3D(bins, limits);
    
    double value, uncertainty;
    int err;
    
    // Test with std::array overload
    // Position (0.5, 0.5, 0.5): all at bin centers, interpolation between all surrounding bins
    std::array<double, 3> pos_arr = {{0.5, 0.5, 0.5}};
    err = grid.extractValue(pos_arr, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Array overload: SUCCESS");
    ts.assertEqual(5.4, value, "Array overload: f(0.5,0.5,0.5) = 5.4");
    
    // Test with std::vector overload
    std::vector<double> pos_vec = {0.5, 0.5, 0.5};
    err = grid.extractValue(pos_vec, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Vector overload: SUCCESS");
    ts.assertEqual(5.4, value, "Vector overload: f(0.5,0.5,0.5) = 5.4");
    
    // Test with three scalars overload (3D)
    err = grid.extractValue(0.5, 0.5, 0.5, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Three scalars overload: SUCCESS");
    ts.assertEqual(5.4, value, "Three scalars overload: f(0.5,0.5,0.5) = 5.4");
    
    // Test arbitrary point
    // Position (0.2, 0.7, 0.3)
    err = grid.extractValue(0.2, 0.7, 0.3, value, uncertainty);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Arbitrary point: SUCCESS");
    ts.assertEqual(4.6, value, "f(0.2,0.7,0.3) = 4.6");
    
    // Test uncertainty is positive
    ts.assertTrue(uncertainty > 0.0, "Uncertainty should be positive");
}

// =============== extractSpectrum1D Tests ===============

// Test extractSpectrum1D with 3D grid - various overloads
void testExtractSpectrum1D_3D(TestSuite& ts) {
    ts.startTest("extractSpectrum1D - 3D grid overloads");
    
    std::array<unsigned long, 3> bins = {{10, 10, 10}};
    std::array<std::pair<double, double>, 3> limits = {{
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0)
    }};
    
    auto grid = createLinearGrid3D(bins, limits);
    penred::measurements::results<double, 1> spectrum;
    int err;
    
    // Test with std::array overload
    std::array<double, 2> pos_arr = {{0.5, 0.5}};
    err = grid.extractSpectrum1D(0, pos_arr, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Array overload: SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(10), spectrum.getNBins(), 
                  "Array overload: correct number of bins");
    
    // Test with std::vector overload
    std::vector<double> pos_vec = {0.5, 0.5};
    err = grid.extractSpectrum1D(0, pos_vec, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Vector overload: SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(10), spectrum.getNBins(), 
                  "Vector overload: correct number of bins");
    
    // Test with two scalars overload (3D)
    err = grid.extractSpectrum1D(0, 0.5, 0.5, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Two scalars overload: SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(10), spectrum.getNBins(), 
                  "Two scalars overload: correct number of bins");
    
    for(unsigned long i = 0; i < 10; ++i) {
        double expected = 0.2 * i + 3.95;
        ts.assertEqual(expected, spectrum.data[i], 
                      "Spectrum[" + std::to_string(i) + "] = " + std::to_string(expected));
    }
    
    // Test vector dimension mismatch
    std::vector<double> bad_vec = {0.5};  // Should be size 2 for 3D grid
    err = grid.extractSpectrum1D(0, bad_vec, spectrum);
    ts.assertEqual(penred::measurements::errors::DIMENSION_MISMATCH, err, 
                  "Vector size mismatch should return error");
    
    // Test invalid spectrum dimension
    std::array<double, 2> pos_valid = {{0.5, 0.5}};
    err = grid.extractSpectrum1D(5, pos_valid, spectrum);  // 5 >= 3
    ts.assertEqual(penred::measurements::errors::DIMENSION_OUT_OF_RANGE, err, 
                  "Invalid dimension should return error");
}

// Test extractSpectrum1D with 2D grid overloads
void testExtractSpectrum1D_2D(TestSuite& ts) {
    ts.startTest("extractSpectrum1D - 2D grid overloads");
    
    std::array<unsigned long, 2> bins = {{8, 6}};
    std::array<std::pair<double, double>, 2> limits = {{
        std::pair<double, double>(0.0, 8.0),
        std::pair<double, double>(0.0, 6.0)
    }};
    
    auto grid = createLinearGrid2D(bins, limits);
    penred::measurements::results<double, 1> spectrum;
    int err;
    
    // Test with std::array overload (extract along x at y=3.0)
    // y=3.0 is exactly at left edge of bin 3 (bin width = 1.0)
    // So no interpolation in y, just bin 3 values
    std::array<double, 1> pos_arr = {{3.0}};
    err = grid.extractSpectrum1D(0, pos_arr, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Array overload: SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(8), spectrum.getNBins(), 
                  "Correct number of bins for x dimension");
    
    // Test with std::vector overload
    std::vector<double> pos_vec = {3.0};
    err = grid.extractSpectrum1D(0, pos_vec, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Vector overload: SUCCESS");
    
    // Test with single scalar overload (2D)
    err = grid.extractSpectrum1D(0, 3.0, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Single scalar overload: SUCCESS");
    
    // Check values: f(x_i, 3.0) where x_i is bin center
    // Bin i: center at i+0.5, value = 2.0*(i+0.5) + 3.0*(3.0+0.5)
    // = 2.0*i + 1.0 + 10.5 = 2.0*i + 11.5
    for(unsigned long i = 0; i < 8; ++i) {
        double expected = 2.0 * i + 11.5;
        ts.assertEqual(expected, spectrum.data[i], 
                      "Spectrum[" + std::to_string(i) + "] = " + std::to_string(expected));
    }
    
    // Test extracting along y at x=4.0
    err = grid.extractSpectrum1D(1, 4.0, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Extract along y: SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(6), spectrum.getNBins(), 
                  "Correct number of bins for y dimension");
    
    // f(4.0, y_j) = 2.0*(4.0+0.5) + 3.0*(j+0.5) = 9.0 + 3.0*j + 1.5 = 3.0*j + 10.5
    for(unsigned long i = 0; i < 6; ++i) {
        double expected = 3.0 * i + 10.5;
        ts.assertEqual(expected, spectrum.data[i], 
                      "Spectrum y[" + std::to_string(i) + "] = " + std::to_string(expected));
    }
}

// Test extractSpectrum1D with 3D grid - different spectrum dimensions
void testExtractSpectrum1D_DifferentDimensions(TestSuite& ts) {
    ts.startTest("extractSpectrum1D - different spectrum dimensions");
    
    std::array<unsigned long, 3> bins = {{5, 5, 5}};
    std::array<std::pair<double, double>, 3> limits = {{
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0)
    }};
    
    auto grid = createLinearGrid3D(bins, limits);
    penred::measurements::results<double, 1> spectrum;
    
    // Extract along dimension 1 (y) at x=0.5, z=0.3
    // x=0.5: between bin 2 (left edge 0.4) and bin 3 (left edge 0.6), fraction=0.5
    // z=0.3: between bin 1 (left edge 0.2) and bin 2 (left edge 0.4), fraction=0.5
    // For each y bin j: value = average of 4 corners in x,z
    // = 0.25*(value at 2,j,1 + value at 2,j,2 + value at 3,j,1 + value at 3,j,2)
    // = 0.4*2.5 + 0.6*j + 0.8*1.5 + 0.9 = 1.0 + 0.6*j + 1.2 + 0.9 = 0.6*j + 3.1
    int err = grid.extractSpectrum1D(1, 0.5, 0.3, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Extract along y: SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(5), spectrum.getNBins(), "5 bins along y");
    
    for(unsigned long i = 0; i < 5; ++i) {
        double expected = 0.6 * i + 3.1;
        ts.assertEqual(expected, spectrum.data[i], 
                      "Spectrum y[" + std::to_string(i) + "] = " + std::to_string(expected));
    }
    
    // Extract along dimension 2 (z) at x=0.2, y=0.7
    // x=0.2: exactly at left edge of bin 1, fraction=0.0
    // y=0.7: between bin 3 (left edge 0.6) and bin 4 (left edge 0.8), fraction=0.5
    // For each z bin k: value = average of 2 corners in y
    // = 0.4*1 + 0.6*3.5 + 0.8*k + 0.9 = 0.4 + 2.1 + 0.8*k + 0.9 = 0.8*k + 3.4
    err = grid.extractSpectrum1D(2, 0.2, 0.7, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Extract along z: SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(5), spectrum.getNBins(), "5 bins along z");
    
    for(unsigned long i = 0; i < 5; ++i) {
        double expected = 0.8 * i + 3.4;
        ts.assertEqual(expected, spectrum.data[i], 
                      "Spectrum z[" + std::to_string(i) + "] = " + std::to_string(expected));
    }
}

// Test edge cases
void testExtractSpectrum1D_EdgeCases(TestSuite& ts) {
    ts.startTest("extractSpectrum1D - edge cases");
    
    std::array<unsigned long, 3> bins = {{4, 4, 4}};
    std::array<std::pair<double, double>, 3> limits = {{
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0)
    }};
    
    auto grid = createLinearGrid3D(bins, limits);
    penred::measurements::results<double, 1> spectrum;
    
    // Test at position outside lower bounds
    int err = grid.extractSpectrum1D(0, -0.5, -0.5, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Should handle out-of-bounds gracefully");
    ts.assertEqual(static_cast<unsigned long>(4), spectrum.getNBins(), "Should still create spectrum");
    
    // Test at position outside upper bounds
    err = grid.extractSpectrum1D(0, 1.5, 1.5, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "Should handle out-of-bounds gracefully");
    ts.assertEqual(static_cast<unsigned long>(4), spectrum.getNBins(), "Should still create spectrum");
}

// Test uncertainty propagation
void testExtractSpectrum1D_Uncertainty(TestSuite& ts) {
    ts.startTest("extractSpectrum1D - uncertainty propagation");
    
    std::array<unsigned long, 2> bins = {{3, 3}};
    std::array<std::pair<double, double>, 2> limits = {{
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0)
    }};
    
    penred::measurements::results<double, 2> grid;
    grid.init(bins, limits);
    
    // Fill with constant values but varying uncertainties
    for(unsigned long i = 0; i < 9; ++i) {
        grid.data[i] = 1.0;
        grid.sigma[i] = (i % 3) * 0.1 + 0.1;
    }
    
    // Extract spectrum along x at y=0.5
    penred::measurements::results<double, 1> spectrum;
    int err = grid.extractSpectrum1D(0, 0.5, spectrum);
    
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "SUCCESS");
    
    // Check that uncertainties are positive and reasonable
    for(unsigned long i = 0; i < 3; ++i) {
        ts.assertTrue(spectrum.sigma[i] > 0.0, 
                     "Uncertainty[" + std::to_string(i) + "] should be positive");
        ts.assertTrue(spectrum.sigma[i] < 0.5, 
                     "Uncertainty[" + std::to_string(i) + "] should be reasonable");
    }
}

// Test header preservation
void testExtractSpectrum1D_Headers(TestSuite& ts) {
    ts.startTest("extractSpectrum1D - header preservation");
    
    std::array<unsigned long, 3> bins = {{5, 5, 5}};
    std::array<std::pair<double, double>, 3> limits = {{
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0),
        std::pair<double, double>(0.0, 1.0)
    }};
    
    penred::measurements::results<double, 3> grid;
    grid.init(bins, limits);
    
    // Set custom headers
    grid.setDimHeader(0, "Energy");
    grid.setDimHeader(1, "Angle");
    grid.setDimHeader(2, "Time");
    grid.setValueHeader("Cross Section");
    grid.setSigmaHeader("Uncertainty");
    
    // Extract spectrum along Energy (dim 0)
    penred::measurements::results<double, 1> spectrum;
    int err = grid.extractSpectrum1D(0, 0.5, 0.5, spectrum);
    
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "SUCCESS");
    ts.assertEqual(std::string("Energy"), spectrum.readDimHeader(0), "Dimension header preserved");
    ts.assertEqual(std::string("Cross Section"), spectrum.readValueHeader(), "Value header preserved");
    ts.assertEqual(std::string("Uncertainty"), spectrum.readSigmaHeader(), "Sigma header preserved");
    
    // Extract spectrum along Angle (dim 1)
    err = grid.extractSpectrum1D(1, 0.5, 0.5, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "SUCCESS");
    ts.assertEqual(std::string("Angle"), spectrum.readDimHeader(0), "Dimension header preserved for Angle");
}

// Test binning consistency
void testExtractSpectrum1D_Binning(TestSuite& ts) {
    ts.startTest("extractSpectrum1D - binning consistency");
    
    std::array<unsigned long, 3> bins = {{8, 6, 4}};
    std::array<std::pair<double, double>, 3> limits = {{
        std::pair<double, double>(0.0, 8.0),
        std::pair<double, double>(0.0, 6.0),
        std::pair<double, double>(0.0, 4.0)
    }};
    
    auto grid = createLinearGrid3D(bins, limits);
    penred::measurements::results<double, 1> spectrum;
    
    int err = grid.extractSpectrum1D(0, 3.0, 2.0, spectrum);
    ts.assertEqual(penred::measurements::errors::SUCCESS, err, "SUCCESS");
    ts.assertEqual(static_cast<unsigned long>(8), spectrum.getNBins(), "8 bins");
    ts.assertEqual(0.0, spectrum.readLimits()[0].first, "Lower limit 0.0");
    ts.assertEqual(8.0, spectrum.readLimits()[0].second, "Upper limit 8.0");
    ts.assertEqual(1.0, spectrum.readBinWidth(0), "Bin width 1.0");
}

int main() {
    TestSuite ts;
    
    std::cout << "Running extraction tests (C++14)..." << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    // extractValue tests
    testExtractValue1D(ts);
    testExtractValue2D(ts);
    testExtractValue3D(ts);
    
    // extractSpectrum1D tests
    testExtractSpectrum1D_3D(ts);
    testExtractSpectrum1D_2D(ts);
    testExtractSpectrum1D_DifferentDimensions(ts);
    testExtractSpectrum1D_EdgeCases(ts);
    testExtractSpectrum1D_Uncertainty(ts);
    testExtractSpectrum1D_Headers(ts);
    testExtractSpectrum1D_Binning(ts);
    
    ts.printSummary();
    
    return ts.getTestsFailed() > 0 ? 1 : 0;
}
