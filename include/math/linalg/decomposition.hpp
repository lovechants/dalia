#ifndef MATH_LINALG_DECOMPOSITION_HPP
#define MATH_LINALG_DECOMPOSITION_HPP

#include "../core/matrix.hpp"
#include "../core/vector.hpp"
#include <cmath>
#include <utility>
#include <optional>

namespace math::linalg {

template<concepts::Arithmetic T, std::size_t N>
struct LUDecomposition {
    Matrix<T, N, N> L;
    Matrix<T, N, N> U;
    std::array<std::size_t, N> P;
    bool singular;
};

template<concepts::Arithmetic T, std::size_t N>
LUDecomposition<T, N> lu_decompose(const Matrix<T, N, N>& A) {
    LUDecomposition<T, N> result;
    result.L = Matrix<T, N, N>::identity();
    result.U = A;
    result.singular = false;
    
    for (std::size_t i = 0; i < N; ++i) {
        result.P[i] = i;
    }
    
    constexpr T epsilon = std::numeric_limits<T>::epsilon() * T{100};
    
    for (std::size_t k = 0; k < N; ++k) {
        T max_val = std::abs(result.U(k, k));
        std::size_t pivot_row = k;
        
        for (std::size_t i = k + 1; i < N; ++i) {
            T val = std::abs(result.U(i, k));
            if (val > max_val) {
                max_val = val;
                pivot_row = i;
            }
        }
        
        if (max_val < epsilon) {
            result.singular = true;
            return result;
        }
        
        if (pivot_row != k) {
            for (std::size_t j = 0; j < N; ++j) {
                std::swap(result.U(k, j), result.U(pivot_row, j));
                if (j < k) {
                    std::swap(result.L(k, j), result.L(pivot_row, j));
                }
            }
            std::swap(result.P[k], result.P[pivot_row]);
        }
        
        for (std::size_t i = k + 1; i < N; ++i) {
            T factor = result.U(i, k) / result.U(k, k);
            result.L(i, k) = factor;
            
            for (std::size_t j = k; j < N; ++j) {
                result.U(i, j) -= factor * result.U(k, j);
            }
        }
    }
    
    return result;
}

template<concepts::Arithmetic T, std::size_t N>
struct QRDecomposition {
    Matrix<T, N, N> Q;
    Matrix<T, N, N> R;
};

template<concepts::Arithmetic T, std::size_t N>
QRDecomposition<T, N> qr_decompose(const Matrix<T, N, N>& A) {
    QRDecomposition<T, N> result;
    result.Q = Matrix<T, N, N>::zeros();
    result.R = Matrix<T, N, N>::zeros();
    
    std::array<Vector<T, N>, N> q;
    
    for (std::size_t j = 0; j < N; ++j) {
        Vector<T, N> a_j;
        for (std::size_t i = 0; i < N; ++i) {
            a_j[i] = A(i, j);
        }
        
        Vector<T, N> u_j = a_j;
        
        for (std::size_t i = 0; i < j; ++i) {
            T proj = dot(a_j, q[i]);
            result.R(i, j) = proj;
            for (std::size_t k = 0; k < N; ++k) {
                u_j[k] -= proj * q[i][k];
            }
        }
        
        T norm_u = norm(u_j);
        result.R(j, j) = norm_u;
        
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T{100};
        if (norm_u > epsilon) {
            q[j] = u_j / norm_u;
        } else {
            for (std::size_t i = 0; i < N; ++i) {
                q[j][i] = T{0};
            }
        }
    }
    
    for (std::size_t j = 0; j < N; ++j) {
        for (std::size_t i = 0; i < N; ++i) {
            result.Q(i, j) = q[j][i];
        }
    }
    
    return result;
}

template<concepts::Arithmetic T, std::size_t N>
struct CholeskyDecomposition {
    Matrix<T, N, N> L;
    bool positive_definite;
};

template<concepts::Arithmetic T, std::size_t N>
CholeskyDecomposition<T, N> cholesky_decompose(const Matrix<T, N, N>& A) {
    CholeskyDecomposition<T, N> result;
    result.L = Matrix<T, N, N>::zeros();
    result.positive_definite = true;
    
    constexpr T epsilon = std::numeric_limits<T>::epsilon() * T{100};
    
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            T sum = T{0};
            
            if (j == i) {
                for (std::size_t k = 0; k < j; ++k) {
                    sum += result.L(j, k) * result.L(j, k);
                }
                T diag = A(j, j) - sum;
                
                if (diag <= epsilon) {
                    result.positive_definite = false;
                    return result;
                }
                
                result.L(j, j) = std::sqrt(diag);
            } else {
                for (std::size_t k = 0; k < j; ++k) {
                    sum += result.L(i, k) * result.L(j, k);
                }
                result.L(i, j) = (A(i, j) - sum) / result.L(j, j);
            }
        }
    }
    
    return result;
}

}

#endif
