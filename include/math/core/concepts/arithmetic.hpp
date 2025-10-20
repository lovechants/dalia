#ifndef MATH_CORE_CONCEPTS_ARITHMETIC_HPP
#define MATH_CORE_CONCEPTS_ARITHMETIC_HPP

#include <concepts>

namespace math::concepts {

template<typename T>
concept Arithmetic = std::floating_point<T> || std::integral<T>;

template<typename T>
concept FloatingPoint = std::floating_point<T>;

template<typename T>
concept Integral = std::integral<T>;

template<typename T>
concept Numeric = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
};

template<std::size_t N>
concept NonZero = N > 0;

}

#endif
