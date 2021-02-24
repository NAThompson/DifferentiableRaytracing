#ifndef HALLEY_TEST_HPP
#define HALLEY_TEST_HPP
#include <drt/halley.hpp>

using namespace drt;
using Real = double;

TEST(HalleyTest, Division) {
    Real a = 1;
    Real b = 2;
    auto f = [&](Real x) {
        Real y = a*x - b;
        Real dydx = a;
        Real d2ydx2 = 0;
        return std::make_tuple(y, dydx, d2ydx2);
    };

    Real bdiva = halley<Real>(f, -2.0, 2.0);
    EXPECT_FLOAT_EQ(bdiva, 2.0);

    bdiva = halley<Real>(f, -2.0, 1.0);
    EXPECT_TRUE(std::isnan(bdiva));
}


TEST(HalleyTest, SquareRoot) {
    auto f = [&](Real x) {
        Real y = x*x - 2;
        Real dydx = 2*x;
        return std::make_tuple(y, dydx, 2);
    };

    // It should get the minimal t root:
    Real neg_rt_two = halley<Real>(f, -2.0, 2.0);
    EXPECT_FLOAT_EQ(neg_rt_two, -sqrt(2.0));

    Real rt_two = halley<Real>(f, 1, 2.0);
    EXPECT_FLOAT_EQ(rt_two, sqrt(2.0));

    // This one starts at a place where f'(0) = 0.
    rt_two = halley<Real>(f, 0.0, 2.0);
    EXPECT_FLOAT_EQ(rt_two, sqrt(2.0));
}


TEST(HalleyTest, NoSolution) {
    auto f = [&](Real x) {
        Real y = x*x + 2;
        Real dydx = 2*x;
        return std::make_tuple(y, dydx, 2);
    };

    // No solution!
    Real rt = halley<Real>(f, -4.0, 4.0);
    EXPECT_TRUE(std::isnan(rt));
}


TEST(HalleyTest, CubeRoot) {
    auto f = [&](Real t) {
        Real y = std::cbrt(t);
        Real dydx;
        if (y == 0) {
            dydx = std::numeric_limits<Real>::infinity();
        } else {
            dydx = -1/(3*y*y);
        }
        Real d2ydx2 = -2*std::pow(y, -5)/9;
        return std::make_tuple(y, dydx, d2ydx2);
    };

    
    Real rt = halley<Real>(f, -4.0, 4.0);
    EXPECT_TRUE(abs(rt) <= sqrt(std::numeric_limits<Real>::epsilon()));
}


#endif
