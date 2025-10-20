#include <math/core/io.hpp>
#include "test_framework.hpp"


TEST(basic_print) {
    using namespace math::io;
    println("test");
    math::test::assert_true(true);
}

RUN_ALL_TESTS()
