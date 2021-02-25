#ifndef NEWTON_TEST_HPP
#define NEWTON_TEST_HPP
#include <drt/newton.hpp>

using namespace drt;
using Real = double;

TEST(NewtonTest, Division) {
    Real a = 1;
    Real b = 2;
    auto f = [&](Real x) {
        Real y = a*x - b;
        Real dydx = a;
        return std::make_pair(y, dydx);
    };

    Real bdiva = newton<Real>(f, -2.0, 2.0);
    EXPECT_FLOAT_EQ(bdiva, 2.0);

    bdiva = newton<Real>(f, -2.0, 1.0);
    EXPECT_TRUE(std::isnan(bdiva));
}

TEST(NewtonTest, SquareRoot) {
    auto f = [&](Real x) {
        Real y = x*x - 2;
        Real dydx = 2*x;
        return std::make_pair(y, dydx);
    };

    // It should get the minimal t root:
    Real neg_rt_two = newton<Real>(f, -2.0, 2.0);
    EXPECT_FLOAT_EQ(neg_rt_two, -sqrt(2.0));

    Real rt_two = newton<Real>(f, 1, 2.0);
    EXPECT_FLOAT_EQ(rt_two, sqrt(2.0));

    // This one starts at a place where f'(0) = 0.
    rt_two = newton<Real>(f, 0.0, 2.0);
    EXPECT_FLOAT_EQ(rt_two, sqrt(2.0));
}

TEST(NewtonTest, NoSolution) {
    auto f = [&](Real x) {
        Real y = x*x + 2;
        Real dydx = 2*x;
        return std::make_pair(y, dydx);
    };

    // No solution!
    Real rt = newton<Real>(f, -4.0, 4.0);
    EXPECT_TRUE(std::isnan(rt));
}

TEST(NewtonTest, CubeRoot) {
    auto f = [&](Real t) {
        Real y = std::cbrt(t);
        Real dydx;
        if (y == 0) {
            dydx = std::numeric_limits<Real>::infinity();
        } else {
            dydx = -1/(3*y*y);
        }
        return std::make_pair(y, dydx);
    };

    // Famously Newton's method doesn't converge here.
    Real rt = newton<Real>(f, -4.0, 4.0);
    EXPECT_TRUE(abs(rt) <= sqrt(std::numeric_limits<Real>::epsilon()));
}


TEST(NewtonTest, 2D) {
    // From Annie Cuyt's paper
    // "Computational Implementation of the Multivariate Halley Method for Solving Nonlinear Systems of Equations."
    auto f = [](Real t, Real u) {
        vec<Real, 2> v;
        Real expmp = exp(-t+u);
        Real expmm = exp(-t-u);
        v[0] = expmp - 0.1;
        v[1] = expmm - 0.1;

        mat<Real, 2,2> M;
        M(0,0) = -expmp;
        M(0,1) = expmp;
        M(1,0) = -expmm;
        M(1,1) = -expmm;
        return std::make_pair(v, M);
    };

    auto [t, u] = newton<Real>(f, 0.0, 5.0, -1.0, 3.0);
    EXPECT_FLOAT_EQ(t, -log(0.1));
    EXPECT_LE(abs(u), std::numeric_limits<Real>::epsilon());
}

#endif
