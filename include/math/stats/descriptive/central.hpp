#ifndef MATH_STATS_DESCRIPTIVE_CENTRAL_HPP
#define MATH_STATS_DESCRIPTIVE_CENTRAL_HPP

#include "../../core/concepts/arithmetic.hpp"
#include "../../core/vector.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

namespace math::stats::descriptive {

template<concepts::Arithmetic T>
T mean(const std::vector<T>& data) {
    if (data.empty()) {
        return T{0};
    }
    
    T sum = T{0};
    for (const auto& x : data) {
        sum += x;
    }
    return sum / static_cast<T>(data.size());
}

template<concepts::Arithmetic T, std::size_t N>
T mean(const Vector<T, N>& data) {
    T sum = T{0};
    for (std::size_t i = 0; i < N; ++i) {
        sum += data[i];
    }
    return sum / static_cast<T>(N);
}

template<concepts::Arithmetic T>
T median(std::vector<T> data) {
    if (data.empty()) {
        return T{0};
    }
    
    std::sort(data.begin(), data.end());
    std::size_t n = data.size();
    
    if (n % 2 == 0) {
        return (data[n/2 - 1] + data[n/2]) / T{2};
    } else {
        return data[n/2];
    }
}

template<concepts::Arithmetic T>
std::vector<T> mode(const std::vector<T>& data) {
    if (data.empty()) {
        return {};
    }
    
    std::vector<T> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    std::vector<T> modes;
    T current = sorted[0];
    std::size_t current_count = 1;
    std::size_t max_count = 1;
    
    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i] == current) {
            current_count++;
        } else {
            if (current_count > max_count) {
                max_count = current_count;
                modes.clear();
                modes.push_back(current);
            } else if (current_count == max_count) {
                modes.push_back(current);
            }
            current = sorted[i];
            current_count = 1;
        }
    }
    
    if (current_count > max_count) {
        modes.clear();
        modes.push_back(current);
    } else if (current_count == max_count) {
        modes.push_back(current);
    }
    
    return modes;
}

template<concepts::Arithmetic T>
T quantile(std::vector<T> data, double q) {
    if (data.empty() || q < 0.0 || q > 1.0) {
        return T{0};
    }
    
    std::sort(data.begin(), data.end());
    
    if (q == 0.0) return data.front();
    if (q == 1.0) return data.back();
    
    double pos = q * (data.size() - 1);
    std::size_t lower = static_cast<std::size_t>(std::floor(pos));
    std::size_t upper = static_cast<std::size_t>(std::ceil(pos));
    
    if (lower == upper) {
        return data[lower];
    }
    
    double weight = pos - lower;
    return data[lower] * (1.0 - weight) + data[upper] * weight;
}

}

#endif
