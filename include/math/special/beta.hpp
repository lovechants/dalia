#ifndef MATH_SPECIAL_BETA_HPP
#define MATH_SPECIAL_BETA_HPP

#include "../core/concepts/arithmetic.hpp"
#include "gamma.hpp"
#include <cmath>

namespace math::special {

template<concepts::FloatingPoint T>
T beta(T a, T b) {
    return std::exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b));
}

template<concepts::FloatingPoint T>
T log_beta(T a, T b) {
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b);
}


// TODO -> These betas return a seg fault "passing test" but seg fault so tests are commented out
template<concepts::FloatingPoint T>
T incomplete_beta(T x, T a, T b) {
    if (x < T{0} || x > T{1}) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (x == T{0}) {
        return T{0};
    }
    if (x == T{1}) {
        return T{1};
    }
    
    T bt = std::exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) +
                    a * std::log(x) + b * std::log(T{1} - x));
    
    if (x < (a + T{1}) / (a + b + T{2})) {
        T result = T{0};
        T term = T{1};
        
        for (int n = 0; n < 100; ++n) {
            T num = T{1};
            T den = T{1};
            
            for (int i = 0; i < n; ++i) {
                num *= (a + static_cast<T>(i)) * x;
                den *= (a + b + static_cast<T>(i)) * (static_cast<T>(i) + T{1});
            }
            
            term = num / den;
            result += term;
            
            if (std::abs(term) < T{1e-10}) {
                break;
            }
        }
        
        return bt * result / a;
    } else {
        return T{1} - incomplete_beta(T{1} - x, b, a);
    }
}

template<concepts::FloatingPoint T>
T regularized_incomplete_beta(T x, T a, T b) {
    return incomplete_beta(x, a, b) / beta(a, b);
}

}

#endif
