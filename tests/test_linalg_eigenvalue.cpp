#include <math/linalg/eigenvalue.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::linalg;

TEST(power_iteration_symmetric) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 2.0; A(0, 1) = 1.0; A(0, 2) = 0.0;
    A(1, 0) = 1.0; A(1, 1) = 2.0; A(1, 2) = 1.0;
    A(2, 0) = 0.0; A(2, 1) = 1.0; A(2, 2) = 2.0;
    
    auto [eigenvalue, eigenvector] = power_iteration(A);
    
    auto Av = A * eigenvector;
    auto lambda_v = eigenvalue * eigenvector;
    
    for (std::size_t i = 0; i < 3; ++i) {
        assert_near(Av[i], lambda_v[i], 1e-6);
    }
}

TEST(qr_algorithm_symmetric) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 4.0; A(0, 1) = 1.0; A(0, 2) = 0.0;
    A(1, 0) = 1.0; A(1, 1) = 3.0; A(1, 2) = 1.0;
    A(2, 0) = 0.0; A(2, 1) = 1.0; A(2, 2) = 2.0;
    
    auto result = qr_algorithm(A);
    assert_true(result.converged);
    
    for (std::size_t i = 0; i < 3; ++i) {
        Vector<double, 3> v;
        for (std::size_t j = 0; j < 3; ++j) {
            v[j] = result.eigenvectors(j, i);
        }
        
        auto Av = A * v;
        auto lambda_v = result.eigenvalues[i] * v;
        
        for (std::size_t j = 0; j < 3; ++j) {
            assert_near(Av[j], lambda_v[j], 1e-6);
        }
    }
}

TEST(rayleigh_quotient) {
    Matrix<double, 3, 3> A;
    A(0, 0) = 2.0; A(0, 1) = 0.0; A(0, 2) = 0.0;
    A(1, 0) = 0.0; A(1, 1) = 3.0; A(1, 2) = 0.0;
    A(2, 0) = 0.0; A(2, 1) = 0.0; A(2, 2) = 5.0;
    
    Vec3<double> v(0.0, 1.0, 0.0);
    double rq = rayleigh_quotient(A, v);
    
    assert_near(rq, 3.0, 1e-10);
}

RUN_ALL_TESTS()
