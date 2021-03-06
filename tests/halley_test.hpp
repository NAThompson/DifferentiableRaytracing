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

TEST(HalleyTest, CuytExpProblem)
{
    // Cuyt, "Computational Implementation of the Multivariate Halley Method for Solving Nonlinear Systems of Equations",
    // Equation 4.8
    // exp(-x₀ + x₁) = d₀
    // exp(-x₀ - x₁) = d₁
    // Has solution x₀ = -ln(d₁d₀)/2, x₁ = ln(d₀/d₁)/2.
    Real d0 = 0.1;
    Real d1 = 0.2;
    auto f = [&d0, &d1](vec<Real, 2> const & x) {
        vec<Real,2> y;
        Real expmp = exp(-x[0] + x[1]);
        Real expmm = exp(-x[0] - x[1]);
        y[0] = expmp - d0;
        y[1] = expmm - d1;

        matrix<Real,2,2> J;
        J(0,0) = -expmp;
        J(0,1) = expmp;
        J(1,0) = -expmm;
        J(1,1) = -expmm;

        tensor<Real, 2,2,2> H;
        H(0,0,0) = expmp;
        H(0,0,1) = -expmp;
        H(0,1,0) = -expmp;
        H(0,1,1) = expmp;
        H(1,0,0) = expmm;
        H(1,0,1) = expmm;
        H(1,1,0) = expmm;
        H(1,1,1) = expmm;
        return std::make_tuple(y, J, H);
    };

    bounds<Real, 2> b({-5,5}, {-5,5});

    auto sol = halley<Real,2>(f, b, b.center());
    vec<Real,2> expected_sol(-std::log(d0*d1)/2, std::log(d0/d1)/2);
    EXPECT_FLOAT_EQ(sol[0], expected_sol[0]);
    EXPECT_FLOAT_EQ(sol[1], expected_sol[1]);
}


#endif
