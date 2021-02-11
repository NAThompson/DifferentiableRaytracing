#include <gtest/gtest.h>
#include "roots_test.hpp"
#include "torus_test.hpp"
#include "helicoid_test.hpp"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
