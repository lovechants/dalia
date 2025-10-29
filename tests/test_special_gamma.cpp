#include <math/special/gamma.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::special;

TEST(gamma_integers) {
    assert_near(gamma(1.0), 1.0, 1e-10);
    assert_near(gamma(2.0), 1.0, 1e-10);
    assert_near(gamma(3.0), 2.0, 1e-10);
    assert_near(gamma(4.0), 6.0, 1e-10);
    assert_near(gamma(5.0), 24.0, 1e-10);
}

TEST(gamma_half_integers) {
    assert_near(gamma(0.5), std::sqrt(std::numbers::pi), 1e-10);
    assert_near(gamma(1.5), 0.5 * std::sqrt(std::numbers::pi), 1e-10);
}

TEST(log_gamma) {
    assert_near(log_gamma(1.0), 0.0, 1e-10);
    assert_near(log_gamma(2.0), 0.0, 1e-10);
    assert_near(log_gamma(10.0), std::log(362880.0), 1e-8);
}

TEST(factorial) {
    assert_near(factorial<double>(0), 1.0, 1e-10);
    assert_near(factorial<double>(5), 120.0, 1e-10);
    assert_near(factorial<double>(10), 3628800.0, 1e-6);
}

TEST(binomial_coefficient) {
    assert_near(binomial_coefficient<double>(5, 2), 10.0, 1e-10);
    assert_near(binomial_coefficient<double>(10, 3), 120.0, 1e-10);
    assert_near(binomial_coefficient<double>(20, 10), 184756.0, 1e-6);
}

RUN_ALL_TESTS()
