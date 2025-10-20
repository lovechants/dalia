#include <math/core/vector.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;

TEST(vector_construction) {
    Vec3<double> v1(1.0, 2.0, 3.0);
    assert_eq(v1[0], 1.0);
    assert_eq(v1[1], 2.0);
    assert_eq(v1[2], 3.0);
}

TEST(vector_addition) {
    Vec3<double> v1(1.0, 2.0, 3.0);
    Vec3<double> v2(4.0, 5.0, 6.0);
    auto v3 = v1 + v2;
    assert_eq(v3[0], 5.0);
    assert_eq(v3[1], 7.0);
    assert_eq(v3[2], 9.0);
}

TEST(vector_dot_product) {
    Vec3<double> v1(1.0, 2.0, 3.0);
    Vec3<double> v2(4.0, 5.0, 6.0);
    double result = dot(v1, v2);
    assert_eq(result, 32.0);
}

TEST(vector_cross_product) {
    Vec3<double> v1(1.0, 0.0, 0.0);
    Vec3<double> v2(0.0, 1.0, 0.0);
    auto v3 = cross(v1, v2);
    assert_eq(v3[0], 0.0);
    assert_eq(v3[1], 0.0);
    assert_eq(v3[2], 1.0);
}

TEST(vector_norm) {
    Vec3<double> v(3.0, 4.0, 0.0);
    double n = norm(v);
    assert_near(n, 5.0, 1e-10);
}

RUN_ALL_TESTS()
