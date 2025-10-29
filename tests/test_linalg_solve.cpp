#include <math/linalg/solve.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::linalg;

TEST(forward_substitution) {
    Matrix<double, 3, 3> L;
    L(0, 0) = 1.0; L(0, 1) = 0.0; L(0, 2) = 0.0;
    L(1, 0) = 2.0; L(1, 1) = 1.0; L(1, 2) = 0.0;
    L(2, 0) = 3.0; L(2, 1) = 4.0; L(2, 2) = 1.0;
    
    Vec3<double> b(1.0, 5.0, 18.0);
    auto x = forward_substitution(L, b);
    
    auto result = L * x;
    assert_near(result[0], b[0], 1e-10);
    assert_near(result[1], b[1], 1e-10);
    assert_near(result[2], b[2], 1e-10);
}

TEST(backward_substitution) {
    Matrix<double, 3, 3> U;
    U(0, 0) = 1.0; U(0, 1) = 2.0; U(0, 2) = 3.0;
    U(1, 0) = 0.0; U(1, 1) = 1.0; U(1, 2) = 4.0;
    U(2, 0) = 0.0; U(2, 1) = 0.0; U(2, 2) = 1.0;
    
    Vec3<double> b(14.0, 11.0, 2.0);
    auto x = backward_substitution(U, b);
    
    auto result = U * x;
    assert_near(result[0], b[0], 1e-10);
    assert_near(result[1], b[1], 1e-10);
    assert_near(result[2], b[2], 1e-10);
}

TEST(solve_linear_system) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 2.0; A(0, 1) = 1.0; A(0, 2) = 1.0;
    A(1, 0) = 4.0; A(1, 1) = 3.0; A(1, 2) = 3.0;
    A(2, 0) = 8.0; A(2, 1) = 7.0; A(2, 2) = 9.0;
    
    Vec3<double> b(4.0, 10.0, 24.0);
    
    auto x = solve(A, b);
    assert_true(x.has_value());
    
    auto result = A * (*x);
    assert_near(result[0], b[0], 1e-10);
    assert_near(result[1], b[1], 1e-10);
    assert_near(result[2], b[2], 1e-10);
}

TEST(least_squares) {
    Matrix<double, 4, 2> A;
    A(0, 0) = 1.0; A(0, 1) = 1.0;
    A(1, 0) = 1.0; A(1, 1) = 2.0;
    A(2, 0) = 1.0; A(2, 1) = 3.0;
    A(3, 0) = 1.0; A(3, 1) = 4.0;
    
    Vec4<double> b(2.0, 3.0, 5.0, 6.0);
    
    auto x = least_squares(A, b);
    
    Matrix<double, 4, 2> A_copy = A;
    auto residual = A_copy * x;
    
    double error = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        error += (residual[i] - b[i]) * (residual[i] - b[i]);
    }
    
    assert_true(error < 1.0);
}

TEST(matrix_inverse) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;
    A(1, 0) = 0.0; A(1, 1) = 1.0; A(1, 2) = 4.0;
    A(2, 0) = 5.0; A(2, 1) = 6.0; A(2, 2) = 0.0;
    
    auto inv = matrix_inverse(A);
    assert_true(inv.has_value());
    
    auto identity = A * (*inv);
    
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            assert_near(identity(i, j), expected, 1e-10);
        }
    }
}

RUN_ALL_TESTS()
