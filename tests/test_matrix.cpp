#include <math/core/matrix.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;

TEST(matrix_construction) {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0;
    m(0, 1) = 2.0;
    m(1, 0) = 3.0;
    m(1, 1) = 4.0;
    assert_eq(m(0, 0), 1.0);
    assert_eq(m(1, 1), 4.0);
}

TEST(matrix_identity) {
    auto m = Matrix<double, 3, 3>::identity();
    assert_eq(m(0, 0), 1.0);
    assert_eq(m(1, 1), 1.0);
    assert_eq(m(2, 2), 1.0);
    assert_eq(m(0, 1), 0.0);
}

TEST(matrix_addition) {
    Matrix<double, 2, 2> m1, m2;
    m1(0, 0) = 1.0; m1(0, 1) = 2.0;
    m1(1, 0) = 3.0; m1(1, 1) = 4.0;
    m2(0, 0) = 5.0; m2(0, 1) = 6.0;
    m2(1, 0) = 7.0; m2(1, 1) = 8.0;
    auto m3 = m1 + m2;
    assert_eq(m3(0, 0), 6.0);
    assert_eq(m3(1, 1), 12.0);
}

TEST(matrix_multiplication) {
    Matrix<double, 2, 2> m1, m2;
    m1(0, 0) = 1.0; m1(0, 1) = 2.0;
    m1(1, 0) = 3.0; m1(1, 1) = 4.0;
    m2(0, 0) = 5.0; m2(0, 1) = 6.0;
    m2(1, 0) = 7.0; m2(1, 1) = 8.0;
    auto m3 = m1 * m2;
    assert_eq(m3(0, 0), 19.0);
    assert_eq(m3(0, 1), 22.0);
    assert_eq(m3(1, 0), 43.0);
    assert_eq(m3(1, 1), 50.0);
}

TEST(matrix_vector_multiply) {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;
    Vec2<double> v(5.0, 6.0);
    auto result = m * v;
    assert_eq(result[0], 17.0);
    assert_eq(result[1], 39.0);
}

TEST(matrix_transpose) {
    Matrix<double, 2, 3> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0; m(0, 2) = 3.0;
    m(1, 0) = 4.0; m(1, 1) = 5.0; m(1, 2) = 6.0;
    auto mt = m.transpose();
    assert_eq(mt(0, 0), 1.0);
    assert_eq(mt(1, 0), 2.0);
    assert_eq(mt(2, 1), 6.0);
}

RUN_ALL_TESTS()
