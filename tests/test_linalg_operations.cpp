#include <math/linalg/operations.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::linalg;

TEST(determinant_2x2) {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;
    double det = determinant(m);
    assert_near(det, -2.0, 1e-10);
}

TEST(determinant_3x3) {
    Matrix<double, 3, 3> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0; m(0, 2) = 3.0;
    m(1, 0) = 0.0; m(1, 1) = 1.0; m(1, 2) = 4.0;
    m(2, 0) = 5.0; m(2, 1) = 6.0; m(2, 2) = 0.0;
    double det = determinant(m);
    assert_near(det, 1.0, 1e-10);
}

TEST(inverse_2x2) {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;
    auto inv = inverse(m);
    auto identity = m * inv;
    assert_near(identity(0, 0), 1.0, 1e-10);
    assert_near(identity(1, 1), 1.0, 1e-10);
    assert_near(identity(0, 1), 0.0, 1e-10);
    assert_near(identity(1, 0), 0.0, 1e-10);
}

TEST(inverse_3x3) {
    auto m = Matrix<double, 3, 3>::identity();
    m(0, 1) = 2.0;
    m(1, 2) = 3.0;
    auto inv = inverse(m);
    auto identity = m * inv;
    assert_near(identity(0, 0), 1.0, 1e-10);
    assert_near(identity(1, 1), 1.0, 1e-10);
    assert_near(identity(2, 2), 1.0, 1e-10);
}

RUN_ALL_TESTS()
