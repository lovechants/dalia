#include <math/special/erf.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::special;

TEST(erf_values) {
    assert_near(erf(0.0), 0.0, 1e-10);
    assert_near(erf(1.0), 0.8427007929, 1e-8);
    assert_near(erf(-1.0), -0.8427007929, 1e-8);
    assert_near(erf(2.0), 0.9953222650, 1e-8);
}

TEST(erfc_values) {
    assert_near(erfc(0.0), 1.0, 1e-10);
    assert_near(erfc(1.0), 1.0 - 0.8427007929, 1e-8);
}

TEST(erf_inverse) {
    double result_0 = erf_inv(0.0);
    std::cout << "erf_inv(0.0) = " << result_0 << " (expected 0.0)\n";
    
    double result_05 = erf_inv(0.5);
    std::cout << "erf_inv(0.5) = " << result_05 << " (expected 0.4769362762)\n";
    
    double test_val = 0.8;
    double inv = erf_inv(test_val);
    double forward = erf(inv);
    std::cout << "erf_inv(0.8) = " << inv << ", erf(erf_inv(0.8)) = " << forward << "\n";
    
    assert_near(erf_inv(0.0), 0.0, 1e-6);
    assert_near(erf_inv(0.5), 0.4769362762, 1e-3);
    assert_near(erf(erf_inv(0.8)), 0.8, 1e-3);
}

TEST(normal_cdf) {
    assert_near(normal_cdf(0.0, 0.0, 1.0), 0.5, 1e-8);
    assert_near(normal_cdf(1.0, 0.0, 1.0), 0.8413447460, 1e-8);
    assert_near(normal_cdf(-1.0, 0.0, 1.0), 0.1586552540, 1e-8);
}

TEST(normal_pdf) {
    assert_near(normal_pdf(0.0, 0.0, 1.0), 0.3989422804, 1e-8);
    assert_near(normal_pdf(1.0, 0.0, 1.0), 0.2419707245, 1e-8);
}

RUN_ALL_TESTS()
