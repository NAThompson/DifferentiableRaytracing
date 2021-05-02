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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
