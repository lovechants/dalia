#ifndef MATH_LINALG_OPERATIONS_HPP
#define MATH_LINALG_OPERATIONS_HPP

#include "../core/matrix.hpp"
#include "../core/vector.hpp"

namespace math::linalg {

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
constexpr Matrix<T, Cols, Rows> transpose(const Matrix<T, Rows, Cols>& m) {
    return m.transpose();
}

template<concepts::Arithmetic T, std::size_t N>
constexpr T trace(const Matrix<T, N, N>& m) {
    return m.trace();
}

template<concepts::Arithmetic T, std::size_t N>
constexpr T determinant(const Matrix<T, N, N>& m) {
    if constexpr (N == 1) {
        return m(0, 0);
    } else if constexpr (N == 2) {
        return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
    } else if constexpr (N == 3) {
        return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1))
             - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
             + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    } else {
        T det = T{0};
        for (std::size_t j = 0; j < N; ++j) {
            Matrix<T, N-1, N-1> submatrix;
            for (std::size_t i = 1; i < N; ++i) {
                std::size_t col_idx = 0;
                for (std::size_t k = 0; k < N; ++k) {
                    if (k == j) continue;
                    submatrix(i-1, col_idx) = m(i, k);
                    col_idx++;
                }
            }
            T sign = (j % 2 == 0) ? T{1} : T{-1};
            det += sign * m(0, j) * determinant(submatrix);
        }
        return det;
    }
}

template<concepts::Arithmetic T>
constexpr Matrix<T, 2, 2> inverse(const Matrix<T, 2, 2>& m) {
    T det = determinant(m);
    Matrix<T, 2, 2> result;
    result(0, 0) =  m(1, 1) / det;
    result(0, 1) = -m(0, 1) / det;
    result(1, 0) = -m(1, 0) / det;
    result(1, 1) =  m(0, 0) / det;
    return result;
}

template<concepts::Arithmetic T>
constexpr Matrix<T, 3, 3> inverse(const Matrix<T, 3, 3>& m) {
    T det = determinant(m);
    Matrix<T, 3, 3> result;
    
    result(0, 0) = (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) / det;
    result(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) / det;
    result(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) / det;
    
    result(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) / det;
    result(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) / det;
    result(1, 2) = (m(0, 2) * m(1, 0) - m(0, 0) * m(1, 2)) / det;
    
    result(2, 0) = (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0)) / det;
    result(2, 1) = (m(0, 1) * m(2, 0) - m(0, 0) * m(2, 1)) / det;
    result(2, 2) = (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0)) / det;
    
    return result;
}

}

#endif
