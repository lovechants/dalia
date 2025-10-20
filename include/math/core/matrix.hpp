#ifndef MATH_CORE_MATRIX_HPP
#define MATH_CORE_MATRIX_HPP

#include "concepts/arithmetic.hpp"
#include "concepts/linalg.hpp"
#include "vector.hpp"
#include <array>
#include <iostream>

namespace math {

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols> 
    requires concepts::ValidMatrixDims<Rows, Cols>
class Matrix {
    std::array<T, Rows * Cols> data_;

public:
    using value_type = T;
    static constexpr std::size_t rows_value = Rows;
    static constexpr std::size_t cols_value = Cols;

    constexpr Matrix() : data_{} {}

    constexpr std::size_t rows() const { return Rows; }
    constexpr std::size_t cols() const { return Cols; }

    constexpr T& operator()(std::size_t i, std::size_t j) {
        return data_[i * Cols + j];
    }

    constexpr const T& operator()(std::size_t i, std::size_t j) const {
        return data_[i * Cols + j];
    }

    static constexpr Matrix identity() requires (Rows == Cols) {
        Matrix result;
        for (std::size_t i = 0; i < Rows; ++i) {
            result(i, i) = T{1};
        }
        return result;
    }

    static constexpr Matrix zeros() {
        return Matrix{};
    }

    static constexpr Matrix ones() {
        Matrix result;
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < Cols; ++j) {
                result(i, j) = T{1};
            }
        }
        return result;
    }

    constexpr Matrix operator+(const Matrix& other) const {
        Matrix result;
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < Cols; ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }

    constexpr Matrix operator-(const Matrix& other) const {
        Matrix result;
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < Cols; ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

    constexpr Matrix operator*(T scalar) const {
        Matrix result;
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < Cols; ++j) {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }
        return result;
    }

    template<std::size_t OtherCols>
    constexpr Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols>& other) const {
        Matrix<T, Rows, OtherCols> result;
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < OtherCols; ++j) {
                T sum = T{0};
                for (std::size_t k = 0; k < Cols; ++k) {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    constexpr Vector<T, Rows> operator*(const Vector<T, Cols>& v) const {
        Vector<T, Rows> result;
        for (std::size_t i = 0; i < Rows; ++i) {
            T sum = T{0};
            for (std::size_t j = 0; j < Cols; ++j) {
                sum += (*this)(i, j) * v[j];
            }
            result[i] = sum;
        }
        return result;
    }

    constexpr Matrix<T, Cols, Rows> transpose() const {
        Matrix<T, Cols, Rows> result;
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < Cols; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    constexpr T trace() const requires (Rows == Cols) {
        T sum = T{0};
        for (std::size_t i = 0; i < Rows; ++i) {
            sum += (*this)(i, i);
        }
        return sum;
    }

    constexpr bool operator==(const Matrix& other) const {
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < Cols; ++j) {
                if ((*this)(i, j) != other(i, j)) return false;
            }
        }
        return true;
    }
};

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
constexpr Matrix<T, Rows, Cols> operator*(T scalar, const Matrix<T, Rows, Cols>& m) {
    return m * scalar;
}

template<concepts::Arithmetic T, std::size_t Rows, std::size_t Cols>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Rows, Cols>& m) {
    os << "[";
    for (std::size_t i = 0; i < Rows; ++i) {
        if (i > 0) os << " ";
        os << "[";
        for (std::size_t j = 0; j < Cols; ++j) {
            os << m(i, j);
            if (j < Cols - 1) os << ", ";
        }
        os << "]";
        if (i < Rows - 1) os << "\n";
    }
    os << "]";
    return os;
}

template<concepts::Arithmetic T, std::size_t N>
using SquareMatrix = Matrix<T, N, N>;

}

#endif
