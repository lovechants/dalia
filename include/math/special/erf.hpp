#ifndef MATH_SPECIAL_ERF_HPP
#define MATH_SPECIAL_ERF_HPP

#include "../core/concepts/arithmetic.hpp"
#include <cmath>
#include <numbers>
#include <limits>

namespace math::special {

template<concepts::FloatingPoint T>
T erf(T x) {
    constexpr T a1 = T{0.254829592};
    constexpr T a2 = T{-0.284496736};
    constexpr T a3 = T{1.421413741};
    constexpr T a4 = T{-1.453152027};
    constexpr T a5 = T{1.061405429};
    constexpr T p = T{0.3275911};
    
    int sign = (x < T{0}) ? -1 : 1;
    x = std::abs(x);
    
    T t = T{1} / (T{1} + p * x);
    T y = T{1} - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);
    
    return sign * y;
}

template<concepts::FloatingPoint T>
T erfc(T x) {
    return T{1} - erf(x);
}

template<concepts::FloatingPoint T>
T erf_inv(T x) {
    // Handle edge cases
    if (x <= T{-1}) {
        return x == T{-1} ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::quiet_NaN();
    }
    if (x >= T{1}) {
        return x == T{1} ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::quiet_NaN();
    }
    
    // Handle very small values directly to avoid precision issues
    if (std::abs(x) < T{1e-10}) return x; // erf_inv(x) approx x for small x
    
    // Save the sign of x
    int sign = (x < T{0}) ? -1 : 1;
    x = std::abs(x);
    
    // Use a rational approximation for better accuracy
    // This is from Winitzki's approximation
    T a = T{0.147};
    T log_1mx2 = std::log(T{1} - x * x);
    T term1 = T{2} / (std::numbers::pi_v<T> * a) + log_1mx2 / T{2};
    T term2 = log_1mx2 / a;
    T result = std::sqrt(-term1 + std::sqrt(term1 * term1 - term2));
    
    // Apply the saved sign to the result
    return sign * result;
}

template<concepts::FloatingPoint T>
T normal_cdf(T x, T mu = T{0}, T sigma = T{1}) {
    return T{0.5} * (T{1} + erf((x - mu) / (sigma * std::sqrt(T{2}))));
}

template<concepts::FloatingPoint T>
T normal_pdf(T x, T mu = T{0}, T sigma = T{1}) {
    T z = (x - mu) / sigma;
    return std::exp(-T{0.5} * z * z) / (sigma * std::sqrt(T{2} * std::numbers::pi_v<T>));
}

}

#endif
