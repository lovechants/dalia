#include <math/stats/descriptive/central.hpp>
#include <math/stats/descriptive/dispersion.hpp>
#include <math/stats/descriptive/correlation.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::stats::descriptive;

TEST(mean_vector) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    assert_near(mean(data), 3.0, 1e-10);
}

TEST(median_odd) {
    std::vector<double> data = {1.0, 3.0, 5.0, 7.0, 9.0};
    assert_near(median(data), 5.0, 1e-10);
}

TEST(median_even) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
    assert_near(median(data), 2.5, 1e-10);
}

TEST(variance_sample) {
    std::vector<double> data = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
    assert_near(variance(data, true), 4.571428571, 1e-8);
}

TEST(std_dev_sample) {
    std::vector<double> data = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
    assert_near(std_dev(data, true), 2.138089935, 1e-8);
}

TEST(range_test) {
    std::vector<double> data = {1.0, 5.0, 3.0, 9.0, 2.0};
    assert_near(range(data), 8.0, 1e-10);
}

TEST(covariance_test) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {2.0, 4.0, 5.0, 4.0, 5.0};
    assert_near(covariance(x, y, true), 1.5, 1e-10);
}

TEST(correlation_test) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {2.0, 4.0, 6.0, 8.0, 10.0};
    assert_near(correlation(x, y), 1.0, 1e-10);
}

TEST(quantile_test) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    assert_near(quantile(data, 0.5), 3.0, 1e-10);
    assert_near(quantile(data, 0.25), 2.0, 1e-10);
    assert_near(quantile(data, 0.75), 4.0, 1e-10);
}

RUN_ALL_TESTS()
