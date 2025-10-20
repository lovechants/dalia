#include <math/core/concepts.hpp>
#include "test_framework.hpp"

using namespace math::test;

TEST(arithmetic_concept) {
    assert_true(math::Arithmetic<double>);
    assert_true(math::Arithmetic<int>);
    assert_true(math::FloatingPoint<double>);
}

RUN_ALL_TESTS()
