#ifndef MATH_LINALG_EIGENVALUE_HPP
#define MATH_LINALG_EIGENVALUE_HPP

#include "../core/matrix.hpp"
#include "../core/vector.hpp"
#include "decomposition.hpp"
#include "norm.hpp"
#include <cmath>
#include <limits>

namespace math::linalg {

template<concepts::Arithmetic T, std::size_t N>
struct EigenResult {
    std::array<T, N> eigenvalues;
    Matrix<T, N, N> eigenvectors;
    bool converged;
};

template<concepts::Arithmetic T, std::size_t N>
std::pair<T, Vector<T, N>> power_iteration(const Matrix<T, N, N>& A, 
                                            std::size_t max_iter = 1000,
                                            T tolerance = T{1e-10}) {
    Vector<T, N> v;
    for (std::size_t i = 0; i < N; ++i) {
        v[i] = T{1};
    }
    v = normalize(v);
    
    T eigenvalue = T{0};
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        auto v_new = A * v;
        T eigenvalue_new = dot(v, v_new);
        
        v_new = normalize(v_new);
        
        if (std::abs(eigenvalue_new - eigenvalue) < tolerance) {
            return {eigenvalue_new, v_new};
        }
        
        eigenvalue = eigenvalue_new;
        v = v_new;
    }
    
    return {eigenvalue, v};
}

template<concepts::Arithmetic T, std::size_t N>
EigenResult<T, N> qr_algorithm(const Matrix<T, N, N>& A,
                                 std::size_t max_iter = 1000,
                                 T tolerance = T{1e-10}) {
    EigenResult<T, N> result;
    result.converged = false;
    
    Matrix<T, N, N> Ak = A;
    Matrix<T, N, N> Q_total = Matrix<T, N, N>::identity();
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        auto qr = qr_decompose(Ak);
        Ak = qr.R * qr.Q;
        Q_total = Q_total * qr.Q;
        
        T off_diag_norm = T{0};
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                if (i != j) {
                    off_diag_norm += Ak(i, j) * Ak(i, j);
                }
            }
        }
        off_diag_norm = std::sqrt(off_diag_norm);
        
        if (off_diag_norm < tolerance) {
            result.converged = true;
            break;
        }
    }
    
    for (std::size_t i = 0; i < N; ++i) {
        result.eigenvalues[i] = Ak(i, i);
    }
    
    result.eigenvectors = Q_total;
    
    return result;
}

template<concepts::Arithmetic T, std::size_t N>
T rayleigh_quotient(const Matrix<T, N, N>& A, const Vector<T, N>& x) {
    auto Ax = A * x;
    return dot(x, Ax) / dot(x, x);
}

template<concepts::Arithmetic T, std::size_t N>
std::pair<T, Vector<T, N>> inverse_iteration(const Matrix<T, N, N>& A,
                                               T sigma,
                                               std::size_t max_iter = 1000,
                                               T tolerance = T{1e-10}) {
    auto A_shifted = A;
    for (std::size_t i = 0; i < N; ++i) {
        A_shifted(i, i) -= sigma;
    }
    
    Vector<T, N> v;
    for (std::size_t i = 0; i < N; ++i) {
        v[i] = T{1};
    }
    v = normalize(v);
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        auto v_new_opt = solve_lu(A_shifted, v);
        if (!v_new_opt) {
            break;
        }
        
        auto v_new = normalize(*v_new_opt);
        
        if (l2_norm(v_new - v) < tolerance) {
            T eigenvalue = rayleigh_quotient(A, v_new);
            return {eigenvalue, v_new};
        }
        
        v = v_new;
    }
    
    T eigenvalue = rayleigh_quotient(A, v);
    return {eigenvalue, v};
}

}

#endif
