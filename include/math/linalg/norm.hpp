#ifndef MATH_LINALG_NORM_HPP
#define MATH_LINALG_NORM_HPP

#include "../core/vector.hpp"
#include "../core/matrix.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace math::linalg {

template<concepts::Arithmetic T, std::size_t N>
T l1_norm(const Vector<T, N>& v) {
    T sum = T{0};
    for (std::size_t i = 0; i < N; ++i) {
        sum += std::abs(v[i]);
    }
    return sum;
}

template<concepts::Arithmetic T, std::size_t N>
T l2_norm(const Vector<T, N>& v) {
    return v.norm();
}

template<concepts::Arithmetic T, std::size_t N>
T linf_norm(const Vector<T, N>& v) {
    T max_val = std::abs(v[0]);
    for (std::size_t i = 1; i < N; ++i) {
        max_val = std::max(max_val, std::abs(v[i]));
    }
    return max_val;
}

template<concepts::Arithmetic T, std::size_t N, int P>
T lp_norm(const Vector<T, N>& v) {
    static_assert(P > 0, "p must be positive");
    T sum = T{0};
    for (std::size_t i = 0; i < N; ++i) {
        sum += std::pow(std::abs(v[i]), P);
    }
    return std::pow(sum, T{1} / P);
}

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
T frobenius_norm(const Matrix<T, Rows, Cols>& m) {
    T sum = T{0};
    for (std::size_t i = 0; i < Rows; ++i) {
        for (std::size_t j = 0; j < Cols; ++j) {
            sum += m(i, j) * m(i, j);
        }
    }
    return std::sqrt(sum);
}

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
T max_norm(const Matrix<T, Rows, Cols>& m) {
    T max_val = std::abs(m(0, 0));
    for (std::size_t i = 0; i < Rows; ++i) {
        for (std::size_t j = 0; j < Cols; ++j) {
            max_val = std::max(max_val, std::abs(m(i, j)));
        }
    }
    return max_val;
}

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
T matrix_1_norm(const Matrix<T, Rows, Cols>& m) {
    T max_col_sum = T{0};
    for (std::size_t j = 0; j < Cols; ++j) {
        T col_sum = T{0};
        for (std::size_t i = 0; i < Rows; ++i) {
            col_sum += std::abs(m(i, j));
        }
        max_col_sum = std::max(max_col_sum, col_sum);
    }
    return max_col_sum;
}

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
T matrix_inf_norm(const Matrix<T, Rows, Cols>& m) {
    T max_row_sum = T{0};
    for (std::size_t i = 0; i < Rows; ++i) {
        T row_sum = T{0};
        for (std::size_t j = 0; j < Cols; ++j) {
            row_sum += std::abs(m(i, j));
        }
        max_row_sum = std::max(max_row_sum, row_sum);
    }
    return max_row_sum;
}

}

#endif
