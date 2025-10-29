#ifndef MATH_STATS_DESCRIPTIVE_CORRELATION_HPP
#define MATH_STATS_DESCRIPTIVE_CORRELATION_HPP

#include "../../core/concepts/arithmetic.hpp"
#include "../../core/vector.hpp"
#include "central.hpp"
#include "dispersion.hpp"
#include <vector>
#include <cmath>

namespace math::stats::descriptive {

template<concepts::Arithmetic T>
T covariance(const std::vector<T>& x, const std::vector<T>& y, bool sample = true) {
    if (x.size() != y.size() || x.empty() || (sample && x.size() == 1)) {
        return T{0};
    }
    
    T mean_x = mean(x);
    T mean_y = mean(y);
    
    T sum = T{0};
    for (std::size_t i = 0; i < x.size(); ++i) {
        sum += (x[i] - mean_x) * (y[i] - mean_y);
    }
    
    std::size_t n = sample ? x.size() - 1 : x.size();
    return sum / static_cast<T>(n);
}

template<concepts::Arithmetic T>
T correlation(const std::vector<T>& x, const std::vector<T>& y) {
    if (x.size() != y.size() || x.size() < 2) {
        return T{0};
    }
    
    T cov = covariance(x, y, true);
    T std_x = std_dev(x, true);
    T std_y = std_dev(y, true);
    
    if (std_x == T{0} || std_y == T{0}) {
        return T{0};
    }
    
    return cov / (std_x * std_y);
}

template<concepts::Arithmetic T>
T spearman_correlation(std::vector<T> x, std::vector<T> y) {
    if (x.size() != y.size() || x.size() < 2) {
        return T{0};
    }
    
    std::vector<std::size_t> rank_x(x.size());
    std::vector<std::size_t> rank_y(y.size());
    
    std::vector<std::size_t> indices_x(x.size());
    std::vector<std::size_t> indices_y(y.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        indices_x[i] = i;
        indices_y[i] = i;
    }
    
    std::sort(indices_x.begin(), indices_x.end(),
              [&x](std::size_t i, std::size_t j) { return x[i] < x[j]; });
    std::sort(indices_y.begin(), indices_y.end(),
              [&y](std::size_t i, std::size_t j) { return y[i] < y[j]; });
    
    for (std::size_t i = 0; i < x.size(); ++i) {
        rank_x[indices_x[i]] = i + 1;
        rank_y[indices_y[i]] = i + 1;
    }
    
    std::vector<T> rank_x_t(x.size());
    std::vector<T> rank_y_t(y.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        rank_x_t[i] = static_cast<T>(rank_x[i]);
        rank_y_t[i] = static_cast<T>(rank_y[i]);
    }
    
    return correlation(rank_x_t, rank_y_t);
}

}

#endif
