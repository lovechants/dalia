#ifndef MATH_LINALG_SOLVE_HPP
#define MATH_LINALG_SOLVE_HPP

#include "../core/matrix.hpp"
#include "../core/vector.hpp"
#include "decomposition.hpp"
#include <optional>
#include <limits>

namespace math::linalg {

template<concepts::Arithmetic T, std::size_t N>
Vector<T, N> forward_substitution(const Matrix<T, N, N>& L, const Vector<T, N>& b) {
    Vector<T, N> x;
    
    for (std::size_t i = 0; i < N; ++i) {
        T sum = b[i];
        for (std::size_t j = 0; j < i; ++j) {
            sum -= L(i, j) * x[j];
        }
        x[i] = sum / L(i, i);
    }
    
    return x;
}

template<concepts::Arithmetic T, std::size_t N>
Vector<T, N> backward_substitution(const Matrix<T, N, N>& U, const Vector<T, N>& b) {
    Vector<T, N> x;
    
    for (std::size_t i = N; i-- > 0; ) {
        T sum = b[i];
        for (std::size_t j = i + 1; j < N; ++j) {
            sum -= U(i, j) * x[j];
        }
        x[i] = sum / U(i, i);
    }
    
    return x;
}

template<concepts::Arithmetic T, std::size_t N>
std::optional<Vector<T, N>> solve_lu(const Matrix<T, N, N>& A, const Vector<T, N>& b) {
    auto lu = lu_decompose(A);
    
    if (lu.singular) {
        return std::nullopt;
    }
    
    Vector<T, N> b_permuted;
    for (std::size_t i = 0; i < N; ++i) {
        b_permuted[i] = b[lu.P[i]];
    }
    
    auto y = forward_substitution(lu.L, b_permuted);
    auto x = backward_substitution(lu.U, y);
    
    return x;
}

template<concepts::Arithmetic T, std::size_t N>
std::optional<Vector<T, N>> solve_cholesky(const Matrix<T, N, N>& A, const Vector<T, N>& b) {
    auto chol = cholesky_decompose(A);
    
    if (!chol.positive_definite) {
        return std::nullopt;
    }
    
    auto y = forward_substitution(chol.L, b);
    auto x = backward_substitution(chol.L.transpose(), y);
    
    return x;
}

template<concepts::Arithmetic T, std::size_t N>
std::optional<Vector<T, N>> solve(const Matrix<T, N, N>& A, const Vector<T, N>& b) {
    return solve_lu(A, b);
}

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
Vector<T, Cols> least_squares(const Matrix<T, Rows, Cols>& A, const Vector<T, Rows>& b) {
    auto At = A.transpose();
    auto AtA = At * A;
    auto Atb = At * b;
    
    auto result = solve(AtA, Atb);
    
    if (result) {
        return *result;
    }
    
    Vector<T, Cols> zero;
    for (std::size_t i = 0; i < Cols; ++i) {
        zero[i] = T{0};
    }
    return zero;
}

template<concepts::Arithmetic T, std::size_t N>
std::optional<Matrix<T, N, N>> matrix_inverse(const Matrix<T, N, N>& A) {
    auto lu = lu_decompose(A);
    
    if (lu.singular) {
        return std::nullopt;
    }
    
    Matrix<T, N, N> inv;
    
    for (std::size_t j = 0; j < N; ++j) {
        Vector<T, N> e;
        for (std::size_t i = 0; i < N; ++i) {
            e[i] = (lu.P[i] == j) ? T{1} : T{0};
        }
        
        auto y = forward_substitution(lu.L, e);
        auto x = backward_substitution(lu.U, y);
        
        for (std::size_t i = 0; i < N; ++i) {
            inv(i, j) = x[i];
        }
    }
    
    return inv;
}

}

#endif
