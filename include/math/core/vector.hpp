#ifndef MATH_CORE_VECTOR_HPP
#define MATH_CORE_VECTOR_HPP

#include "concepts/arithmetic.hpp"
#include <array>
#include <cmath>
#include <iostream>

namespace math {

template<concepts::Arithmetic T, std::size_t N>
    requires concepts::NonZero<N>
class Vector {
    std::array<T, N> data_;

public:
    using value_type = T;
    static constexpr std::size_t size_value = N;

    constexpr Vector() : data_{} {}

    constexpr Vector(std::array<T, N> data) : data_(data) {}

    template<typename... Args> requires (sizeof...(Args) == N)
    constexpr Vector(Args... args) : data_{static_cast<T>(args)...} {}

    constexpr std::size_t size() const { return N; }

    constexpr T& operator[](std::size_t i) { return data_[i]; }
    constexpr const T& operator[](std::size_t i) const { return data_[i]; }

    constexpr T& operator()(std::size_t i) { return data_[i]; }
    constexpr const T& operator()(std::size_t i) const { return data_[i]; }

    constexpr auto begin() { return data_.begin(); }
    constexpr auto end() { return data_.end(); }
    constexpr auto begin() const { return data_.begin(); }
    constexpr auto end() const { return data_.end(); }

    constexpr Vector operator+(const Vector& other) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = data_[i] + other[i];
        }
        return result;
    }

    constexpr Vector operator-(const Vector& other) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = data_[i] - other[i];
        }
        return result;
    }

    constexpr Vector operator*(T scalar) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = data_[i] * scalar;
        }
        return result;
    }

    constexpr Vector operator/(T scalar) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i) {
            result[i] = data_[i] / scalar;
        }
        return result;
    }

    constexpr Vector& operator+=(const Vector& other) {
        for (std::size_t i = 0; i < N; ++i) {
            data_[i] += other[i];
        }
        return *this;
    }

    constexpr Vector& operator-=(const Vector& other) {
        for (std::size_t i = 0; i < N; ++i) {
            data_[i] -= other[i];
        }
        return *this;
    }

    constexpr Vector& operator*=(T scalar) {
        for (std::size_t i = 0; i < N; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    constexpr Vector& operator/=(T scalar) {
        for (std::size_t i = 0; i < N; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    constexpr bool operator==(const Vector& other) const {
        for (std::size_t i = 0; i < N; ++i) {
            if (data_[i] != other[i]) return false;
        }
        return true;
    }

    constexpr T dot(const Vector& other) const {
        T sum = T{0};
        for (std::size_t i = 0; i < N; ++i) {
            sum += data_[i] * other[i];
        }
        return sum;
    }

    T norm() const {
        return std::sqrt(dot(*this));
    }

    Vector normalized() const {
        T n = norm();
        return *this / n;
    }
};

template<concepts::Arithmetic T, std::size_t N>
constexpr Vector<T, N> operator*(T scalar, const Vector<T, N>& v) {
    return v * scalar;
}

template<concepts::Arithmetic T, std::size_t N>
constexpr T dot(const Vector<T, N>& a, const Vector<T, N>& b) {
    return a.dot(b);
}

template<concepts::Arithmetic T, std::size_t N>
T norm(const Vector<T, N>& v) {
    return v.norm();
}

template<concepts::Arithmetic T, std::size_t N>
Vector<T, N> normalize(const Vector<T, N>& v) {
    return v.normalized();
}

template<concepts::Arithmetic T>
constexpr Vector<T, 3> cross(const Vector<T, 3>& a, const Vector<T, 3>& b) {
    return Vector<T, 3>(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    );
}

template<concepts::Arithmetic T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T, N>& v) {
    os << "[";
    for (std::size_t i = 0; i < N; ++i) {
        os << v[i];
        if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
}

template<concepts::Arithmetic T>
using Vec2 = Vector<T, 2>;

template<concepts::Arithmetic T>
using Vec3 = Vector<T, 3>;

template<concepts::Arithmetic T>
using Vec4 = Vector<T, 4>;

}

#endif
