#ifndef MATH_STATS_DESCRIPTIVE_DISPERSION_HPP
#define MATH_STATS_DESCRIPTIVE_DISPERSION_HPP

#include "../../core/concepts/arithmetic.hpp"
#include "../../core/vector.hpp"
#include "central.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace math::stats::descriptive {

template<concepts::Arithmetic T>
T variance(const std::vector<T>& data, bool sample = true) {
    if (data.empty() || (sample && data.size() == 1)) {
        return T{0};
    }
    
    T mu = mean(data);
    T sum_sq = T{0};
    
    for (const auto& x : data) {
        T diff = x - mu;
        sum_sq += diff * diff;
    }
    
    std::size_t n = sample ? data.size() - 1 : data.size();
    return sum_sq / static_cast<T>(n);
}

template<concepts::Arithmetic T, std::size_t N>
T variance(const Vector<T, N>& data, bool sample = true) {
    if (sample && N == 1) {
        return T{0};
    }
    
    T mu = mean(data);
    T sum_sq = T{0};
    
    for (std::size_t i = 0; i < N; ++i) {
        T diff = data[i] - mu;
        sum_sq += diff * diff;
    }
    
    std::size_t n = sample ? N - 1 : N;
    return sum_sq / static_cast<T>(n);
}

template<concepts::Arithmetic T>
T std_dev(const std::vector<T>& data, bool sample = true) {
    return std::sqrt(variance(data, sample));
}

template<concepts::Arithmetic T, std::size_t N>
T std_dev(const Vector<T, N>& data, bool sample = true) {
    return std::sqrt(variance(data, sample));
}

template<concepts::Arithmetic T>
T range(const std::vector<T>& data) {
    if (data.empty()) {
        return T{0};
    }
    
    auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
    return *max_it - *min_it;
}

template<concepts::Arithmetic T>
T iqr(std::vector<T> data) {
    if (data.size() < 2) {
        return T{0};
    }
    
    T q1 = quantile(data, 0.25);
    T q3 = quantile(data, 0.75);
    return q3 - q1;
}

template<concepts::Arithmetic T>
T mad(const std::vector<T>& data) {
    if (data.empty()) {
        return T{0};
    }
    
    T mu = mean(data);
    std::vector<T> deviations;
    deviations.reserve(data.size());
    
    for (const auto& x : data) {
        deviations.push_back(std::abs(x - mu));
    }
    
    return mean(deviations);
}

}

#endif
