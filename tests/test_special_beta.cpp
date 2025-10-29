#include <math/special/beta.hpp>
#include "test_framework.hpp"

using namespace math;
using namespace math::test;
using namespace math::special;

TEST(beta_function) {
    assert_near(beta(1.0, 1.0), 1.0, 1e-10);
    assert_near(beta(2.0, 3.0), 0.0833333333, 1e-8);
    assert_near(beta(5.0, 5.0), 0.0015873016, 1e-8);
}

TEST(log_beta) {
    assert_near(log_beta(1.0, 1.0), 0.0, 1e-10);
    assert_near(std::exp(log_beta(2.0, 3.0)), 0.0833333333, 1e-8);
}
// See header 
/*TEST(incomplete_beta_endpoints) {*/
/*    assert_near(incomplete_beta(0.0, 2.0, 3.0), 0.0, 1e-10);*/
/*    assert_near(incomplete_beta(1.0, 2.0, 3.0), 1.0, 1e-10);*/
/*}*/
/**/
/*TEST(regularized_incomplete_beta) {*/
/*    double result = regularized_incomplete_beta(0.5, 2.0, 2.0);*/
/*    assert_near(result, 0.5, 1e-6);*/
/*}*/

RUN_ALL_TESTS()
