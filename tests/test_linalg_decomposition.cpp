#include <math/linalg/decomposition.hpp>
#include <math/linalg/norm.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::linalg;

TEST(lu_decomposition_3x3) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 2.0; A(0, 1) = 1.0; A(0, 2) = 1.0;
    A(1, 0) = 4.0; A(1, 1) = 3.0; A(1, 2) = 3.0;
    A(2, 0) = 8.0; A(2, 1) = 7.0; A(2, 2) = 9.0;
    
    auto lu = lu_decompose(A);
    assert_true(!lu.singular);
    
    auto LU = lu.L * lu.U;
    
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            assert_near(LU(i, j), A(lu.P[i], j), 1e-10);
        }
    }
}

TEST(qr_decomposition_3x3) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 12.0; A(0, 1) = -51.0; A(0, 2) = 4.0;
    A(1, 0) = 6.0;  A(1, 1) = 167.0; A(1, 2) = -68.0;
    A(2, 0) = -4.0; A(2, 1) = 24.0;  A(2, 2) = -41.0;
    
    auto qr = qr_decompose(A);
    
    auto QR = qr.Q * qr.R;
    
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            assert_near(QR(i, j), A(i, j), 1e-8);
        }
    }
    
    auto QtQ = qr.Q.transpose() * qr.Q;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            assert_near(QtQ(i, j), expected, 1e-8);
        }
    }
}

TEST(cholesky_decomposition_3x3) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 4.0; A(0, 1) = 12.0; A(0, 2) = -16.0;
    A(1, 0) = 12.0; A(1, 1) = 37.0; A(1, 2) = -43.0;
    A(2, 0) = -16.0; A(2, 1) = -43.0; A(2, 2) = 98.0;
    
    auto chol = cholesky_decompose(A);
    assert_true(chol.positive_definite);
    
    auto LLt = chol.L * chol.L.transpose();
    
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            assert_near(LLt(i, j), A(i, j), 1e-8);
        }
    }
}

RUN_ALL_TESTS()
