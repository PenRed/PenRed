//
// Self-contained test suite for PenRed Animation Library
//

#include "math_classes.hh"
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

using namespace penred;

// Test result structure
struct TestResult {
    std::string name;
    bool passed;
    std::string message;
    
    TestResult(const std::string& n, bool p, const std::string& m = "")
        : name(n), passed(p), message(m) {}
};

// Test runner class
class TestRunner {
private:
    std::vector<TestResult> results;
    int testsRun;
    int testsPassed;
    
public:
    TestRunner() : testsRun(0), testsPassed(0) {}
    
    void runTest(const std::string& name, bool (*testFunc)()) {
        testsRun++;
        bool passed = testFunc();
        if (passed) testsPassed++;
        
        results.emplace_back(name, passed);
        std::cout << (passed ? "PASS" : "FAIL") << " - " << name << std::endl;
    }
    
    template<typename Func>
    void runTest(const std::string& name, Func testFunc) {
        testsRun++;
        bool passed = testFunc();
        if (passed) testsPassed++;
        
        results.emplace_back(name, passed);
        std::cout << (passed ? "PASS" : "FAIL") << " - " << name << std::endl;
    }
    
    void printSummary() const {
        std::cout << "\n=== TEST SUMMARY ===" << std::endl;
        std::cout << "Tests run: " << testsRun << std::endl;
        std::cout << "Tests passed: " << testsPassed << std::endl;
        std::cout << "Tests failed: " << (testsRun - testsPassed) << std::endl;
        std::cout << "Success rate: " << (testsPassed * 100.0 / testsRun) << "%" << std::endl;
        
        if (testsRun != testsPassed) {
            std::cout << "\nFailed tests:" << std::endl;
            for (const auto& result : results) {
                if (!result.passed) {
                    std::cout << "  - " << result.name;
                    if (!result.message.empty()) {
                        std::cout << " (" << result.message << ")";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
    
    bool allPassed() const { return testsRun == testsPassed; }
};

// Assertion macros for tests
#define ASSERT_TRUE(condition) \
    if (!(condition)) { \
        std::cout << "    Assertion failed: " << #condition << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    }

#define ASSERT_FALSE(condition) \
    if ((condition)) { \
        std::cout << "    Assertion failed: " << #condition << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    }

#define ASSERT_EQUAL(a, b) \
    if ((a) != (b)) { \
        std::cout << "    Assertion failed: " << #a << " == " << #b << " (" << (a) << " != " << (b) << ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    }

#define ASSERT_NEAR(a, b, tolerance) \
    if (std::abs((a) - (b)) > (tolerance)) { \
        std::cout << "    Assertion failed: " << #a << " ~= " << #b << " (" << (a) << " != " << (b) << ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return false; \
    }

// Use more appropriate tolerances for different operations
#define ASSERT_APPROX_EQUAL(a, b) ASSERT_NEAR(a, b, 1e-10)
#define ASSERT_ROTATION_EQUAL(a, b) ASSERT_NEAR(a, b, 1e-6)  // More tolerant for rotations

// Test functions
bool testVector3D_BasicOperations() {
    vector3D<double> v1(1.0, 2.0, 3.0);
    vector3D<double> v2(4.0, 5.0, 6.0);
    
    // Addition
    auto v3 = v1 + v2;
    ASSERT_EQUAL(v3.x, 5.0);
    ASSERT_EQUAL(v3.y, 7.0);
    ASSERT_EQUAL(v3.z, 9.0);
    
    // Subtraction
    auto v4 = v2 - v1;
    ASSERT_EQUAL(v4.x, 3.0);
    ASSERT_EQUAL(v4.y, 3.0);
    ASSERT_EQUAL(v4.z, 3.0);
    
    // Dot product
    double dot = v1 * v2;
    ASSERT_EQUAL(dot, 32.0);
    
    // Cross product
    auto cross = v1 ^ v2;
    ASSERT_EQUAL(cross.x, -3.0);
    ASSERT_EQUAL(cross.y, 6.0);
    ASSERT_EQUAL(cross.z, -3.0);
    
    return true;
}

bool testVector3D_Normalization() {
    vector3D<double> v(3.0, 0.0, 0.0);
    v.normalize();
    ASSERT_APPROX_EQUAL(v.x, 1.0);
    ASSERT_APPROX_EQUAL(v.y, 0.0);
    ASSERT_APPROX_EQUAL(v.z, 0.0);
    ASSERT_APPROX_EQUAL(v.mod(), 1.0);
    
    vector3D<double> v2(1.0, 1.0, 1.0);
    v2.normalize();
    double expected = 1.0 / std::sqrt(3.0);
    ASSERT_APPROX_EQUAL(v2.x, expected);
    ASSERT_APPROX_EQUAL(v2.y, expected);
    ASSERT_APPROX_EQUAL(v2.z, expected);
    
    return true;
}

bool testVector3D_LinearInterpolation() {
    vector3D<double> v1(0.0, 0.0, 0.0);
    vector3D<double> v2(10.0, 20.0, 30.0);
    
    auto result = v1.lerp(v2, 0.5);
    ASSERT_APPROX_EQUAL(result.x, 5.0);
    ASSERT_APPROX_EQUAL(result.y, 10.0);
    ASSERT_APPROX_EQUAL(result.z, 15.0);
    
    // Test clamping
    auto below = v1.lerp(v2, -0.5);
    ASSERT_APPROX_EQUAL(below.x, 0.0);
    ASSERT_APPROX_EQUAL(below.y, 0.0);
    ASSERT_APPROX_EQUAL(below.z, 0.0);
    
    auto above = v1.lerp(v2, 1.5);
    ASSERT_APPROX_EQUAL(above.x, 10.0);
    ASSERT_APPROX_EQUAL(above.y, 20.0);
    ASSERT_APPROX_EQUAL(above.z, 30.0);
    
    return true;
}

bool testQuaternion_BasicOperations() {
    Quaternion<double> q1; // Identity
    ASSERT_EQUAL(q1.w, 1.0);
    ASSERT_EQUAL(q1.x, 0.0);
    ASSERT_EQUAL(q1.y, 0.0);
    ASSERT_EQUAL(q1.z, 0.0);
    
    // Axis-angle constructor
    std::array<double, 3> axis = {0.0, 0.0, 1.0};
    Quaternion<double> q2(axis, M_PI/2.0);
    
    // Should be normalized
    double norm = q2.w*q2.w + q2.x*q2.x + q2.y*q2.y + q2.z*q2.z;
    ASSERT_ROTATION_EQUAL(norm, 1.0);
    
    // Conjugate
    auto conj = q2.conjugate();
    ASSERT_APPROX_EQUAL(conj.w, q2.w);
    ASSERT_APPROX_EQUAL(conj.x, -q2.x);
    ASSERT_APPROX_EQUAL(conj.y, -q2.y);
    ASSERT_APPROX_EQUAL(conj.z, -q2.z);
    
    return true;
}

bool testQuaternion_Rotation() {
    // 90-degree rotation around Z-axis
    std::array<double, 3> axis = {0.0, 0.0, 1.0};
    Quaternion<double> q(axis, M_PI/2.0);
    
    // Rotate a vector along X-axis
    vector3D<double> v(1.0, 0.0, 0.0);
    q.rotate(v);
    
    ASSERT_ROTATION_EQUAL(v.x, 0.0);
    ASSERT_ROTATION_EQUAL(v.y, 1.0);
    ASSERT_ROTATION_EQUAL(v.z, 0.0);
    
    // Test with array
    double arr[3] = {1.0, 0.0, 0.0};
    q.rotate(arr);
    ASSERT_ROTATION_EQUAL(arr[0], 0.0);
    ASSERT_ROTATION_EQUAL(arr[1], 1.0);
    ASSERT_ROTATION_EQUAL(arr[2], 0.0);
    
    return true;
}

bool testQuaternion_Slerp() {
    // From 0 to 180 degrees around Z-axis
    std::array<double, 3> axis = {0.0, 0.0, 1.0};
    Quaternion<double> q1(axis, 0.0);
    Quaternion<double> q2(axis, M_PI);
    
    // 50% interpolation should give 90 degrees
    auto q3 = q1.slerp(q2, 0.5);
    
    vector3D<double> v(1.0, 0.0, 0.0);
    q3.rotate(v);
    
    ASSERT_ROTATION_EQUAL(v.x, 0.0);
    ASSERT_ROTATION_EQUAL(v.y, 1.0);
    ASSERT_ROTATION_EQUAL(v.z, 0.0);
    
    return true;
}

bool testTranslation_Basic() {
    transforms::Translation<double> t1({1.0, 2.0, 3.0});
    
    vector3D<double> v(0.0, 0.0, 0.0);
    t1.apply(v);
    ASSERT_EQUAL(v.x, 1.0);
    ASSERT_EQUAL(v.y, 2.0);
    ASSERT_EQUAL(v.z, 3.0);
    
    // Inverse
    t1.applyInv(v);
    ASSERT_EQUAL(v.x, 0.0);
    ASSERT_EQUAL(v.y, 0.0);
    ASSERT_EQUAL(v.z, 0.0);
    
    return true;
}

bool testTranslation_Interpolation() {
    transforms::Translation<double> t1({0.0, 0.0, 0.0});
    transforms::Translation<double> t2({10.0, 20.0, 30.0});
    
    auto t3 = t1.lerp(t2, 0.5);
    vector3D<double> v(0.0, 0.0, 0.0);
    t3.apply(v);
    
    ASSERT_EQUAL(v.x, 5.0);
    ASSERT_EQUAL(v.y, 10.0);
    ASSERT_EQUAL(v.z, 15.0);
    
    return true;
}

bool testRotation_Basic() {
    // 90-degree rotation around Z-axis at origin
    std::array<double, 3> axis = {0.0, 0.0, 1.0};
    Quaternion<double> q(axis, M_PI/2.0);
    transforms::Rotation<double> r(q);
    
    vector3D<double> v(1.0, 0.0, 0.0);
    r.apply(v);
    
    ASSERT_ROTATION_EQUAL(v.x, 0.0);
    ASSERT_ROTATION_EQUAL(v.y, 1.0);
    ASSERT_ROTATION_EQUAL(v.z, 0.0);
    
    // Inverse rotation
    r.applyInv(v);
    ASSERT_ROTATION_EQUAL(v.x, 1.0);
    ASSERT_ROTATION_EQUAL(v.y, 0.0);
    ASSERT_ROTATION_EQUAL(v.z, 0.0);
    
    return true;
}

bool testKeyframe_Basic() {
    transforms::Keyframe<double, double> kf1;
    kf1.frame = 0.0;
    kf1.translation.vector = vector3D<double>(0.0, 0.0, 0.0);
    kf1.rotation.quaternion = Quaternion<double>::identity();
    
    transforms::Keyframe<double, double> kf2;
    kf2.frame = 10.0;
    kf2.translation.vector = vector3D<double>(10.0, 0.0, 0.0);
    std::array<double, 3> axis = {0.0, 0.0, 1.0};
    kf2.rotation.quaternion = Quaternion<double>(axis, M_PI/2.0);
    
    // Test comparison operators
    ASSERT_TRUE(kf1 < kf2);
    ASSERT_TRUE(kf2 > kf1);
    ASSERT_TRUE(kf1 <= kf2);
    ASSERT_TRUE(kf2 >= kf1);
    ASSERT_TRUE(kf1 < 5.0);
    ASSERT_TRUE(kf2 > 5.0);
    
    return true;
}

bool testKeyframe_LinearInterpolation() {
    transforms::Keyframe<double, double> kf1;
    kf1.frame = 0.0;
    kf1.translation.vector = vector3D<double>(0.0, 0.0, 0.0);
    kf1.rotation.quaternion = Quaternion<double>::identity();
    
    transforms::Keyframe<double, double> kf2;
    kf2.frame = 10.0;
    kf2.translation.vector = vector3D<double>(10.0, 0.0, 0.0);
    std::array<double, 3> axis = {0.0, 0.0, 1.0};
    kf2.rotation.quaternion = Quaternion<double>(axis, M_PI/2.0);
    
    auto kf_mid = kf1.linearInterpolate(kf2, 5.0);
    ASSERT_EQUAL(kf_mid.frame, 5.0);
    
    // Test application
    vector3D<double> v(1.0, 0.0, 0.0);
    kf_mid.apply(v);
    
    // At 50% interpolation:
    // - Translation: (5,0,0)
    // - Rotation: 45 degrees
    // Point (1,0,0) rotated 45deg = (cos45, sin45, 0) = (~0.707, ~0.707, 0)
    // Then translated by (5,0,0) = (~5.707, ~0.707, 0)
    ASSERT_ROTATION_EQUAL(v.x, 5.0 + std::cos(M_PI/4.0));
    ASSERT_ROTATION_EQUAL(v.y, std::sin(M_PI/4.0));
    ASSERT_ROTATION_EQUAL(v.z, 0.0);
    
    return true;
}

bool testKeyframe_BezierInterpolation() {
    // Test bezier value calculation
    double v0 = 0.0, v0r = 2.0, v1l = 8.0, v1 = 10.0;
    
    double val0 = transforms::Keyframe<double, double>::bezierValue(0.0, v0, v0r, v1l, v1);
    ASSERT_EQUAL(val0, v0);
    
    double val1 = transforms::Keyframe<double, double>::bezierValue(1.0, v0, v0r, v1l, v1);
    ASSERT_EQUAL(val1, v1);
    
    // Test with linear handlers (should be linear interpolation)
    double val_mid = transforms::Keyframe<double, double>::bezierValue(0.5, 0.0, 2.5, 7.5, 10.0);
    ASSERT_APPROX_EQUAL(val_mid, 5.0);
    
    return true;
}

bool testAnimation_Basic() {
    transforms::Animation<double, double> animation;
    
    // Parse from string stream - simple translation animation
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "0.0 0 0 0 1 0 0 0 \n";          // Start at origin, identity rotation
    ss << "10.0 10 0 0 1 0 0 0\n";        // End at (10,0,0), identity rotation
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 2);
    
    // Test keyframe retrieval and frame values
    auto kf_start = animation.getKeyframe(0.0);
    ASSERT_EQUAL(kf_start.frame, 0.0);
    
    auto kf_mid = animation.getKeyframe(5.0);
    ASSERT_EQUAL(kf_mid.frame, 5.0);
    
    auto kf_end = animation.getKeyframe(10.0);
    ASSERT_EQUAL(kf_end.frame, 10.0);
    
    // Test transformations at keyframes
    vector3D<double> v1(1.0, 2.0, 3.0);
    kf_start.apply(v1);
    ASSERT_EQUAL(v1.x, 1.0);  // No translation at start
    ASSERT_EQUAL(v1.y, 2.0);
    ASSERT_EQUAL(v1.z, 3.0);
    
    vector3D<double> v2(1.0, 2.0, 3.0);
    kf_end.apply(v2);
    ASSERT_EQUAL(v2.x, 11.0); // Full translation at end: 1 + 10
    ASSERT_EQUAL(v2.y, 2.0);
    ASSERT_EQUAL(v2.z, 3.0);
    
    // Test interpolation at midpoint
    vector3D<double> v3(1.0, 2.0, 3.0);
    kf_mid.apply(v3);
    ASSERT_EQUAL(v3.x, 6.0);  // Half translation: 1 + 5
    ASSERT_EQUAL(v3.y, 2.0);
    ASSERT_EQUAL(v3.z, 3.0);
    
    // Test direct transformation function
    vector3D<double> v4(1.0, 2.0, 3.0);
    animation.apply(5.0, v4);
    ASSERT_EQUAL(v4.x, 6.0);  // Should match kf_mid transformation
    ASSERT_EQUAL(v4.y, 2.0);
    ASSERT_EQUAL(v4.z, 3.0);
    
    // Test interpolation at 25% and 75%
    vector3D<double> v5(1.0, 2.0, 3.0);
    animation.apply(2.5, v5);
    ASSERT_EQUAL(v5.x, 3.5);  // 1 + 2.5
    ASSERT_EQUAL(v5.y, 2.0);
    ASSERT_EQUAL(v5.z, 3.0);
    
    vector3D<double> v6(1.0, 2.0, 3.0);
    animation.apply(7.5, v6);
    ASSERT_EQUAL(v6.x, 8.5);  // 1 + 7.5
    ASSERT_EQUAL(v6.y, 2.0);
    ASSERT_EQUAL(v6.z, 3.0);
    
    return true;
}

bool testAnimation_BasicWithRotation() {
    transforms::Animation<double, double> animation;
    
    // Use higher precision quaternion values
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "0.0 0 0 0 1 0 0 0 \n";                    // Start: origin, 0° rotation
    ss << "10.0 10 0 0 0.7071067811865475 0 0 0.7071067811865475\n"; // 90° rotation with higher precision
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 2);
    
    // Test transformations at start
    auto kf_start = animation.getKeyframe(0.0);
    vector3D<double> v_start(1.0, 0.0, 0.0);
    kf_start.apply(v_start);
    ASSERT_EQUAL(v_start.x, 1.0);
    ASSERT_EQUAL(v_start.y, 0.0);
    ASSERT_EQUAL(v_start.z, 0.0);
    
    // Test transformations at end
    auto kf_end = animation.getKeyframe(10.0);
    vector3D<double> v_end(1.0, 0.0, 0.0);
    kf_end.apply(v_end);
    ASSERT_NEAR(v_end.x, 10.0, 1e-6);  // Translated by 10 in X
    ASSERT_NEAR(v_end.y, 1.0, 1e-6);   // Rotated 90°: (1,0,0) -> (0,1,0) but then translated
    ASSERT_NEAR(v_end.z, 0.0, 1e-6);
    
    return true;
}

bool testAnimation_EdgeCases() {
  transforms::Animation<double, double> empty_animation;
  ASSERT_TRUE(empty_animation.empty());
  ASSERT_EQUAL(empty_animation.size(), 0);
    
  // Single keyframe
  transforms::Animation<double, double> single_animation;
  std::stringstream ss;
  ss << "0.0 0.0 0.0\n";
  ss << "5.0 1 2 3 1 0 0 0\n"; // Keyframe at time 5.0 with translation (1,2,3)
  single_animation.parse(ss);
    
  ASSERT_FALSE(single_animation.empty());
  ASSERT_EQUAL(single_animation.size(), 1);
    
  auto kf1 = single_animation.getKeyframe(0.0);  // Before keyframe time
  ASSERT_EQUAL(kf1.frame, 0.0);  // Frame should be the requested value (0.0)
    
  // The transformation should be the identity
  vector3D<double> v1(10.0, 20.0, 30.0);
  kf1.apply(v1);
  ASSERT_EQUAL(v1.x, 10.0); 
  ASSERT_EQUAL(v1.y, 20.0);   
  ASSERT_EQUAL(v1.z, 30.0); 
    
  auto kf2 = single_animation.getKeyframe(5.0);  // At keyframe time
  ASSERT_EQUAL(kf2.frame, 5.0);  // Frame should be the requested value (5.0)
    
  vector3D<double> v2(10.0, 20.0, 30.0);
  kf2.apply(v2);
  ASSERT_EQUAL(v2.x, 11.0); // Same transformation applied
  ASSERT_EQUAL(v2.y, 22.0);
  ASSERT_EQUAL(v2.z, 33.0);
    
  auto kf3 = single_animation.getKeyframe(10.0); // After keyframe time
  ASSERT_EQUAL(kf3.frame, 10.0); // Frame should be the requested value (10.0)
    
  vector3D<double> v3(10.0, 20.0, 30.0);
  kf3.apply(v3);
  ASSERT_EQUAL(v3.x, 11.0); // Same transformation applied
  ASSERT_EQUAL(v3.y, 22.0);
  ASSERT_EQUAL(v3.z, 33.0);
    
  return true;
}

bool testAnimation_SingleKeyframe_Identity() {
    transforms::Animation<double, double> animation;
    
    // Parse single keyframe
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "5.0 1 2 3 1 0 0 0\n"; // Identity rotation, translation (1,2,3)
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 1);
    
    // Test before the keyframe - should return identity with requested frame
    auto kf_before = animation.getKeyframe(0.0);
    ASSERT_EQUAL(kf_before.frame, 0.0);
    
    // Apply identity transformation - should not change the vector
    vector3D<double> v_before(10.0, 20.0, 30.0);
    kf_before.apply(v_before);
    ASSERT_EQUAL(v_before.x, 10.0); // No translation applied
    ASSERT_EQUAL(v_before.y, 20.0);
    ASSERT_EQUAL(v_before.z, 30.0);
    
    return true;
}

bool testAnimation_SingleKeyframe_ExactFrame() {
    transforms::Animation<double, double> animation;
    
    // Parse single keyframe with specific transformation
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "5.0 10 20 30 1 0 0 0\n"; // Translation (10,20,30), identity rotation
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 1);
    
    // Test at the exact keyframe time
    auto kf_exact = animation.getKeyframe(5.0);
    ASSERT_EQUAL(kf_exact.frame, 5.0);
    
    // Apply transformation - should apply the full translation
    vector3D<double> v_exact(1.0, 2.0, 3.0);
    kf_exact.apply(v_exact);
    ASSERT_EQUAL(v_exact.x, 11.0); // 1 + 10
    ASSERT_EQUAL(v_exact.y, 22.0); // 2 + 20
    ASSERT_EQUAL(v_exact.z, 33.0); // 3 + 30
    
    return true;
}

bool testAnimation_SingleKeyframe_AfterFrame() {
    transforms::Animation<double, double> animation;
    
    // Use higher precision quaternion
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "5.0 0 0 0 0.7071067811865475 0 0 0.7071067811865475\n"; // 90deg rotation with higher precision
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 1);
    
    // Test after the keyframe
    auto kf_after = animation.getKeyframe(10.0);
    ASSERT_EQUAL(kf_after.frame, 10.0);
    
    // Apply transformation
    vector3D<double> v_after(1.0, 0.0, 0.0);
    kf_after.apply(v_after);
    ASSERT_NEAR(v_after.x, 0.0, 1e-6);
    ASSERT_NEAR(v_after.y, 1.0, 1e-6);
    ASSERT_NEAR(v_after.z, 0.0, 1e-6);
    
    return true;
}

bool testAnimation_EmptyAnimation() {
    transforms::Animation<double, double> animation;
    ASSERT_TRUE(animation.empty());
    
    // Test various frames with empty animation
    for(double frame : {-10.0, 0.0, 10.0, 100.0}) {
        auto kf = animation.getKeyframe(frame);
        ASSERT_EQUAL(kf.frame, frame); // Should return identity with requested frame
        
        // Identity transformation should not change vectors
        vector3D<double> v(7.0, 8.0, 9.0);
        kf.apply(v);
        ASSERT_EQUAL(v.x, 7.0);
        ASSERT_EQUAL(v.y, 8.0);
        ASSERT_EQUAL(v.z, 9.0);
        
        // Test inverse transformation also does nothing
        vector3D<double> v2(7.0, 8.0, 9.0);
        kf.applyInv(v2);
        ASSERT_EQUAL(v2.x, 7.0);
        ASSERT_EQUAL(v2.y, 8.0);
        ASSERT_EQUAL(v2.z, 9.0);
    }
    
    return true;
}

bool testAnimation_BeforeFirstKeyframe() {
    transforms::Animation<double, double> animation;
    
    // Add keyframes starting at frame 10.0
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "10.0 5 0 0 1 0 0 0\n"; // Translation (5,0,0)
    ss << "20.0 10 0 0 1 0 0 0\n"; // Translation (10,0,0)
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 2);
    
    // Test before first keyframe
    auto kf_before = animation.getKeyframe(5.0);
    ASSERT_EQUAL(kf_before.frame, 5.0);
    
    // Should return identity transformation for frame 5.0
    vector3D<double> v(1.0, 2.0, 3.0);
    kf_before.apply(v);
    ASSERT_EQUAL(v.x, 1.0); // No translation applied
    ASSERT_EQUAL(v.y, 2.0);
    ASSERT_EQUAL(v.z, 3.0);
    
    return true;
}

bool testAnimation_AfterLastKeyframe() {
    transforms::Animation<double, double> animation;
    
    // Add keyframes ending at frame 20.0
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "0.0 0 0 0 1 0 0 0 \n"; // Identity
    ss << "10.0 0 5 0 1 0 0 0 \n"; // Translation (0,5,0)
    ss << "20.0 0 10 0 1 0 0 0 \n"; // Translation (0,10,0)
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 3);
    
    // Test after last keyframe
    auto kf_after = animation.getKeyframe(30.0);
    ASSERT_EQUAL(kf_after.frame, 30.0); // Should return last keyframe with frame value 30
    
    // Should apply the transformation from the last keyframe
    vector3D<double> v(1.0, 2.0, 3.0);
    kf_after.apply(v);
    ASSERT_EQUAL(v.x, 1.0);  // X unchanged
    ASSERT_EQUAL(v.y, 12.0); // 2 + 10 translation in Y
    ASSERT_EQUAL(v.z, 3.0);  // Z unchanged
    
    return true;
}

bool testAnimation_ExactKeyframeMatch() {
    transforms::Animation<double, double> animation;
    
    // Add multiple keyframes
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "0.0 1 0 0 1 0 0 0\n";  // Translation (1,0,0)
    ss << "10.0 0 2 0 1 0 0 0\n"; // Translation (0,2,0)  
    ss << "20.0 0 0 3 1 0 0 0\n"; // Translation (0,0,3)
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 3);
    
    // Test exact keyframe matches
    auto kf0 = animation.getKeyframe(0.0);
    ASSERT_EQUAL(kf0.frame, 0.0);
    vector3D<double> v0(0.0, 0.0, 0.0);
    kf0.apply(v0);
    ASSERT_EQUAL(v0.x, 1.0); // Applied translation from frame 0.0
    
    auto kf10 = animation.getKeyframe(10.0);
    ASSERT_EQUAL(kf10.frame, 10.0);
    vector3D<double> v10(0.0, 0.0, 0.0);
    kf10.apply(v10);
    ASSERT_EQUAL(v10.y, 2.0); // Applied translation from frame 10.0
    
    auto kf20 = animation.getKeyframe(20.0);
    ASSERT_EQUAL(kf20.frame, 20.0);
    vector3D<double> v20(0.0, 0.0, 0.0);
    kf20.apply(v20);
    ASSERT_EQUAL(v20.z, 3.0); // Applied translation from frame 20.0
    
    return true;
}

bool testAnimation_InterpolationAccuracy() {
    transforms::Animation<double, double> animation;
    
    // Simple linear translation animation
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "0.0 0 0 0 1 0 0 0\n";    // Start at origin
    ss << "100.0 100 0 0 1 0 0 0\n"; // End at (100,0,0)
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 2);
    
    // Test multiple interpolation points
    vector3D<double> v(0.0, 0.0, 0.0);
    
    // At 25% - should be at (25,0,0)
    animation.apply(25.0, v);
    ASSERT_APPROX_EQUAL(v.x, 25.0);
    ASSERT_APPROX_EQUAL(v.y, 0.0);
    ASSERT_APPROX_EQUAL(v.z, 0.0);
    
    // At 50% - should be at (50,0,0)
    v = vector3D<double>(0.0, 0.0, 0.0);
    animation.apply(50.0, v);
    ASSERT_APPROX_EQUAL(v.x, 50.0);
    ASSERT_APPROX_EQUAL(v.y, 0.0);
    ASSERT_APPROX_EQUAL(v.z, 0.0);
    
    // At 75% - should be at (75,0,0)
    v = vector3D<double>(0.0, 0.0, 0.0);
    animation.apply(75.0, v);
    ASSERT_APPROX_EQUAL(v.x, 75.0);
    ASSERT_APPROX_EQUAL(v.y, 0.0);
    ASSERT_APPROX_EQUAL(v.z, 0.0);
    
    return true;
}

bool testAnimation_RotationInterpolation() {
    transforms::Animation<double, double> animation;
    
    // Use higher precision quaternions
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    ss << "0.0 0 0 0 1 0 0 0\n";                              // 0 degrees
    ss << "10.0 0 0 0 0.7071067811865475 0 0 0.7071067811865475\n"; // 90 degrees
    ss << "20.0 0 0 0 0 0 0 1\n";                            // 180 degrees
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 3);
    
    vector3D<double> v(1.0, 0.0, 0.0);
    
    // At frame 5.0 - should be 45 degrees
    animation.apply(5.0, v);
    double expected_x_45 = std::cos(M_PI/4.0);
    double expected_y_45 = std::sin(M_PI/4.0);
    ASSERT_NEAR(v.x, expected_x_45, 1e-6);
    ASSERT_NEAR(v.y, expected_y_45, 1e-6);
    ASSERT_NEAR(v.z, 0.0, 1e-6);
    
    // At frame 15.0 - should be 135 degrees  
    v = vector3D<double>(1.0, 0.0, 0.0);
    animation.apply(15.0, v);
    double expected_x_135 = std::cos(3.0*M_PI/4.0);
    double expected_y_135 = std::sin(3.0*M_PI/4.0);
    ASSERT_NEAR(v.x, expected_x_135, 1e-6);
    ASSERT_NEAR(v.y, expected_y_135, 1e-6);
    ASSERT_NEAR(v.z, 0.0, 1e-6);
    
    return true;
}

bool testTriangle_Operations() {
    triangle<double> tri;
    tri.v1 = vector3D<double>(0.0, 0.0, 0.0);
    tri.v2 = vector3D<double>(2.0, 0.0, 0.0);
    tri.v3 = vector3D<double>(0.0, 3.0, 0.0);
    
    ASSERT_EQUAL(tri.minx(), 0.0);
    ASSERT_EQUAL(tri.maxx(), 2.0);
    ASSERT_EQUAL(tri.miny(), 0.0);
    ASSERT_EQUAL(tri.maxy(), 3.0);
    ASSERT_EQUAL(tri.minz(), 0.0);
    ASSERT_EQUAL(tri.maxz(), 0.0);
    
    auto centroid = tri.centroid();
    ASSERT_APPROX_EQUAL(centroid.x, 2.0/3.0);
    ASSERT_APPROX_EQUAL(centroid.y, 1.0);
    ASSERT_APPROX_EQUAL(centroid.z, 0.0);
    
    return true;
}

bool testBox_Operations() {
    box<double> b(0.0, 0.0, 0.0, 2.0, 3.0, 4.0);
    
    ASSERT_EQUAL(b.minx(), 0.0);
    ASSERT_EQUAL(b.maxx(), 2.0);
    ASSERT_EQUAL(b.miny(), 0.0);
    ASSERT_EQUAL(b.maxy(), 3.0);
    ASSERT_EQUAL(b.minz(), 0.0);
    ASSERT_EQUAL(b.maxz(), 4.0);
    
    ASSERT_EQUAL(b.dx(), 2.0);
    ASSERT_EQUAL(b.dy(), 3.0);
    ASSERT_EQUAL(b.dz(), 4.0);
    ASSERT_EQUAL(b.volume(), 24.0);
    
    // Test point containment
    vector3D<double> inside(1.0, 1.0, 1.0);
    vector3D<double> outside(5.0, 5.0, 5.0);
    
    ASSERT_TRUE(b.in(inside, 1e-10));
    ASSERT_FALSE(b.in(outside, 1e-10));
    
    return true;
}

bool testDirectionTransformation_Orthogonality() {
    // Test that direction transformations preserve orthogonality and normalization
    transforms::Keyframe<double, double> kf;
    kf.frame = 0.0;
    kf.translation.vector = vector3D<double>(10.0, 20.0, 30.0); // Translation should not affect directions
    kf.rotation.quaternion = Quaternion<double>(vector3D<double>(0.0, 1.0, 1.0), M_PI/3.0); // Complex rotation

    vector3D<double> dir1(1.0, 0.0, 0.0);
    vector3D<double> dir2(0.0, 1.0, 0.0);
    vector3D<double> dir3(0.0, 0.0, 1.0);

    // Store original properties
    //double original_dot12 = dir1 * dir2;
    //double original_dot13 = dir1 * dir3;
    //double original_dot23 = dir2 * dir3;
    //double original_len1 = dir1.mod();
    //double original_len2 = dir2.mod();
    //double original_len3 = dir3.mod();

    //printf("=== Direction Transformation Test ===\n");
    //printf("Original directions are orthogonal and normalized\n");

    // Apply inverse rotation to directions
    kf.rotation.applyInv(dir1);
    kf.rotation.applyInv(dir2);
    kf.rotation.applyInv(dir3);

    //printf("After inverse rotation:\n");
    //printf("dir1: %s (length: %f)\n", dir1.stringify().c_str(), dir1.mod());
    //printf("dir2: %s (length: %f)\n", dir2.stringify().c_str(), dir2.mod());
    //printf("dir3: %s (length: %f)\n", dir3.stringify().c_str(), dir3.mod());

    // Directions should remain orthogonal and normalized
    ASSERT_NEAR(dir1.mod(), 1.0, 1e-10);
    ASSERT_NEAR(dir2.mod(), 1.0, 1e-10);
    ASSERT_NEAR(dir3.mod(), 1.0, 1e-10);
    
    ASSERT_NEAR(dir1 * dir2, 0.0, 1e-10);
    ASSERT_NEAR(dir1 * dir3, 0.0, 1e-10);
    ASSERT_NEAR(dir2 * dir3, 0.0, 1e-10);

    return true;
}

bool testPerformance_ManyKeyframes() {
    transforms::Animation<double, double> animation;
    
    // Create many keyframes
    std::stringstream ss;
    ss << "0.0 0.0 0.0\n";
    for(int i = 0; i < 100; ++i) {
        ss << static_cast<double>(i) << " " << i << " 0 0 1 0 0 0 0 0 0\n";
    }
    
    size_t num_frames = animation.parse(ss);
    ASSERT_EQUAL(num_frames, 100);
    
    // Test interpolation at multiple points
    for(int i = 0; i <= 99; i += 10) {
        auto kf = animation.getKeyframe(static_cast<double>(i));
        ASSERT_EQUAL(kf.frame, static_cast<double>(i));
        
        vector3D<double> v(1.0, 0.0, 0.0);
        animation.apply(static_cast<double>(i), v);
        ASSERT_EQUAL(v.x, 1.0 + static_cast<double>(i));
    }
    
    return true;
}

int main() {
    std::cout << "Running PenRed Animation Library Tests\n";
    std::cout << "=====================================\n\n";
    
    TestRunner runner;
    
    // Vector3D tests
    runner.runTest("Vector3D Basic Operations", testVector3D_BasicOperations);
    runner.runTest("Vector3D Normalization", testVector3D_Normalization);
    runner.runTest("Vector3D Linear Interpolation", testVector3D_LinearInterpolation);
    
    // Quaternion tests
    runner.runTest("Quaternion Basic Operations", testQuaternion_BasicOperations);
    runner.runTest("Quaternion Rotation", testQuaternion_Rotation);
    runner.runTest("Quaternion SLERP", testQuaternion_Slerp);
    
    // Transform tests
    runner.runTest("Translation Basic", testTranslation_Basic);
    runner.runTest("Translation Interpolation", testTranslation_Interpolation);
    runner.runTest("Rotation Basic", testRotation_Basic);
    
    // Keyframe tests
    runner.runTest("Keyframe Basic", testKeyframe_Basic);
    runner.runTest("Keyframe Linear Interpolation", testKeyframe_LinearInterpolation);
    runner.runTest("Keyframe Bezier Interpolation", testKeyframe_BezierInterpolation);
    
    // Animation tests
    runner.runTest("Animation Basic", testAnimation_Basic);
    runner.runTest("Animation Basic With Rotation", testAnimation_BasicWithRotation);
    runner.runTest("Animation Edge Cases", testAnimation_EdgeCases);
    runner.runTest("Animation Single Keyframe Identity", testAnimation_SingleKeyframe_Identity);
    runner.runTest("Animation Single Keyframe Exact Frame", testAnimation_SingleKeyframe_ExactFrame);
    runner.runTest("Animation Single Keyframe After Frame", testAnimation_SingleKeyframe_AfterFrame);
    runner.runTest("Animation Empty Animation", testAnimation_EmptyAnimation);
    runner.runTest("Animation Before First Keyframe", testAnimation_BeforeFirstKeyframe);
    runner.runTest("Animation After Last Keyframe", testAnimation_AfterLastKeyframe);
    runner.runTest("Animation Exact Keyframe Match", testAnimation_ExactKeyframeMatch);
    runner.runTest("Animation Interpolation Accuracy", testAnimation_InterpolationAccuracy);
    runner.runTest("Animation Rotation Interpolation", testAnimation_RotationInterpolation);
    
    // Geometry tests
    runner.runTest("Triangle Operations", testTriangle_Operations);
    runner.runTest("Box Operations", testBox_Operations);
    
    // Performance tests
    runner.runTest("Performance Many Keyframes", testPerformance_ManyKeyframes);

    // Hierarchical transformation tests
    runner.runTest("Direction Transformation - Orthogonality", testDirectionTransformation_Orthogonality);
    
    runner.printSummary();
    
    return runner.allPassed() ? 0 : 1;
}
