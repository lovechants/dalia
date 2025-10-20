#include <math/linalg/norm.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::linalg;

TEST(l1_norm_vector) {
    Vec3<double> v(1.0, -2.0, 3.0);
    assert_near(l1_norm(v), 6.0, 1e-10);
}

TEST(l2_norm_vector) {
    Vec3<double> v(3.0, 4.0, 0.0);
    assert_near(l2_norm(v), 5.0, 1e-10);
}

TEST(linf_norm_vector) {
    Vec3<double> v(1.0, -5.0, 3.0);
    assert_near(linf_norm(v), 5.0, 1e-10);
}

TEST(frobenius_norm_matrix) {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 2.0; m(1, 1) = 1.0;
    assert_near(frobenius_norm(m), std::sqrt(10.0), 1e-10);
}

TEST(matrix_1_norm) {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 3.0;
    m(1, 0) = 2.0; m(1, 1) = 4.0;
    assert_near(matrix_1_norm(m), 7.0, 1e-10);
}

TEST(matrix_inf_norm) {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;
    assert_near(matrix_inf_norm(m), 7.0, 1e-10);
}

RUN_ALL_TESTS()
