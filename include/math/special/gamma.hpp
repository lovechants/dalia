#ifndef MATH_SPECIAL_GAMMA_HPP
#define MATH_SPECIAL_GAMMA_HPP

#include "../core/concepts/arithmetic.hpp"
#include <cmath>
#include <numbers>
#include <limits>

namespace math::special {

template<concepts::FloatingPoint T>
T lanczos_gamma(T z) {
    constexpr T g = T{7};
    constexpr T coef[] = {
        T{0.99999999999980993},
        T{676.5203681218851},
        T{-1259.1392167224028},
        T{771.32342877765313},
        T{-176.61502916214059},
        T{12.507343278686905},
        T{-0.13857109526572012},
        T{9.9843695780195716e-6},
        T{1.5056327351493116e-7}
    };
    
    if (z < T{0.5}) {
        return std::numbers::pi_v<T> / (std::sin(std::numbers::pi_v<T> * z) * lanczos_gamma(T{1} - z));
    }
    
    z -= T{1};
    T x = coef[0];
    for (std::size_t i = 1; i < 9; ++i) {
        x += coef[i] / (z + static_cast<T>(i));
    }
    
    T t = z + g + T{0.5};
    return std::sqrt(T{2} * std::numbers::pi_v<T>) * std::pow(t, z + T{0.5}) * std::exp(-t) * x;
}

template<concepts::FloatingPoint T>
T gamma(T x) {
    return lanczos_gamma(x);
}

template<concepts::FloatingPoint T>
T log_gamma(T x) {
    if (x <= T{0}) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    
    constexpr T g = T{7};
    constexpr T coef[] = {
        T{0.99999999999980993},
        T{676.5203681218851},
        T{-1259.1392167224028},
        T{771.32342877765313},
        T{-176.61502916214059},
        T{12.507343278686905},
        T{-0.13857109526572012},
        T{9.9843695780195716e-6},
        T{1.5056327351493116e-7}
    };
    
    if (x < T{0.5}) {
        return std::log(std::numbers::pi_v<T>) - std::log(std::sin(std::numbers::pi_v<T> * x)) - log_gamma(T{1} - x);
    }
    
    x -= T{1};
    T sum = coef[0];
    for (std::size_t i = 1; i < 9; ++i) {
        sum += coef[i] / (x + static_cast<T>(i));
    }
    
    T t = x + g + T{0.5};
    return T{0.5} * std::log(T{2} * std::numbers::pi_v<T>) + (x + T{0.5}) * std::log(t) - t + std::log(sum);
}

template<concepts::FloatingPoint T>
T digamma(T x) {
    constexpr T euler_mascheroni = T{0.5772156649015329};
    
    if (x <= T{0}) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    
    T result = -euler_mascheroni;
    T z = x;
    
    while (z < T{10}) {
        result -= T{1} / z;
        z += T{1};
    }
    
    T z2 = z * z;
    result += std::log(z) - T{1} / (T{2} * z) - T{1} / (T{12} * z2);
    result += T{1} / (T{120} * z2 * z2);
    result -= T{1} / (T{252} * z2 * z2 * z2);
    
    return result;
}

template<concepts::FloatingPoint T>
T factorial(int n) {
    if (n < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (n == 0 || n == 1) {
        return T{1};
    }
    return gamma(static_cast<T>(n + 1));
}

template<concepts::FloatingPoint T>
T binomial_coefficient(int n, int k) {
    if (k < 0 || k > n) {
        return T{0};
    }
    if (k == 0 || k == n) {
        return T{1};
    }
    return std::exp(log_gamma(static_cast<T>(n + 1)) - 
                    log_gamma(static_cast<T>(k + 1)) - 
                    log_gamma(static_cast<T>(n - k + 1)));
}

}

#endif
