#include <iostream>
#include <cmath>
#include <cstring>
#include <random>
#include <gtest/gtest.h>
#include "roots_test.hpp"
#include "torus_test.hpp"
using namespace drt;


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
