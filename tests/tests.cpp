#include <gtest/gtest.h>
#include "roots_test.hpp"
#include "torus_test.hpp"
#include "helicoid_test.hpp"
#include "vec_test.hpp"
#include "matrix_test.hpp"
#include "sphere_test.hpp"
#include "ellipsoid_test.hpp"
#include "newton_test.hpp"
#include "halley_test.hpp"
#include "tensor_test.hpp"
#include "ko_method_test.hpp"
#include "interval_test.hpp"

TEST(NanTest, Nans) {
   EXPECT_TRUE(std::isnan(std::numeric_limits<double>::quiet_NaN())) << " You must compile with NaNs enabled.";

   double x = 1.0/0.0;
   EXPECT_TRUE(std::isnan(x)) << " You must compile with NaNs enabled.";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
