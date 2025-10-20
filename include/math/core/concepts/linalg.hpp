#ifndef MATH_CORE_CONCEPTS_LINALG_HPP
#define MATH_CORE_CONCEPTS_LINALG_HPP

#include <cstddef>

namespace math::concepts {

template<std::size_t Rows, std::size_t Cols>
concept ValidMatrixDims = (Rows > 0) && (Cols > 0);

template<typename M>
concept MatrixLike = requires(M m) {
    { m.rows() } -> std::convertible_to<std::size_t>;
    { m.cols() } -> std::convertible_to<std::size_t>;
    { m(std::size_t{}, std::size_t{}) };
};

template<typename V>
concept VectorLike = requires(V v) {
    { v.size() } -> std::convertible_to<std::size_t>;
    { v[std::size_t{}] };
};

}

#endif
