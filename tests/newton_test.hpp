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


TEST(NewtonTest, 2DCuyt) {
    // From Annie Cuyt's paper
    // "Computational Implementation of the Multivariate Halley Method for Solving Nonlinear Systems of Equations."
    auto f = [](vec<Real, 2> const & w) {
        Real t = w[0];
        Real u = w[1];
        vec<Real, 2> v;
        Real expmp = exp(-t+u);
        Real expmm = exp(-t-u);
        v[0] = expmp - 0.1;
        v[1] = expmm - 0.1;

        matrix<Real, 2,2> M;
        M(0,0) = -expmp;
        M(0,1) = expmp;
        M(1,0) = -expmm;
        M(1,1) = -expmm;
        return std::make_pair(v, M);
    };

    bounds<Real, 2> b({0.0, 5.0}, {-1.0, 3.0});
    vec<Real,2> v = newton<Real, 2>(f, b, vec<Real,2>(0.0, 2.0));
    EXPECT_FLOAT_EQ(v[0], -log(0.1));
    EXPECT_LE(abs(v[1]), std::numeric_limits<Real>::epsilon());
}

TEST(NewtonTest, Corless) {
    // From Corless, "A Graduate Introduction to Numerical Methods", Example 3.14.
    auto f = [](vec<Real, 2> w) {
        Real x = w[0];
        Real y = w[1];
        vec<Real, 2> v;
        v[0] = x*x + y*y - 1;
        v[1] = 25*x*y - 12;

        matrix<Real, 2,2> M;
        M(0,0) = 2*x;
        M(0,1) = 2*y;
        M(1,0) = 25*y;
        M(1,1) = 25*x;
        return std::make_pair(v, M);
    };

    bounds<Real, 2> b({0.0, 1.0}, {0.0, 1.0});
    auto v = newton<Real, 2>(f, b, vec<Real,2>(0, 0.5));
    EXPECT_FLOAT_EQ(v[0], 3.0/5);
    EXPECT_FLOAT_EQ(v[1], 4.0/5);
}


TEST(NewtonTest, 2DBroyden) {
    // From Broyden, "A Class of  Methods for  Solving Nonlinear Simultaneous Equations" Case 9:
    // (t,u) = (-0.330435, -0.869239), as well as many other initial conditions, induces backtracking.
    auto f = [](vec<Real,2> w) {
        Real x0 = w[0];
        Real x1 = w[1];
        vec<Real, 2> v;
        v[0] = 10*(x1 - x0*x0);
        v[1] = 1 - x0;

        matrix<Real, 2,2> J;
        J(0,0) = -20*x0;
        J(0,1) = 10;
        J(1,0) = -1;
        J(1,1) = 0;
        return std::make_pair(v, J);
    };
    bounds<Real,2> b({-2.0,2.0}, {-10.0,10.0});

    auto x = newton<Real,2>(f, b, vec<Real,2>(-2,0));
    auto [v, M] = f(x);
    EXPECT_LE(squared_norm(v), std::numeric_limits<Real>::epsilon());
    EXPECT_EQ(x[0], 1);
    EXPECT_EQ(x[1], 1);
}


TEST(NewtonTest, Helicoid) {
        vec<Real, 3> o(-4, 0,1);
        vec<Real, 3> d(1, 0, -1.0/3.0);
        Real speed_ = 1;
        Real radius_ = 1;
        auto f = [&o, &d, &speed_, &radius_](vec<Real,2> w) {
            Real t = w[0];
            Real v = w[1];
            vec<Real, 2> g;
            Real k = 2*M_PI/speed_;
            Real x = o[0] + t*d[0];
            Real y = o[1] + t*d[1];
            Real z = o[2] + t*d[2];
            Real rckz = radius_*cos(k*z);
            Real rskz = radius_*sin(k*z);
            Real rvckz = v*rckz;
            Real rvskz = v*rskz;
            g[0] = x + rvckz;
            g[1] = y + rvskz;
            matrix<Real, 2,2> J;
            J(0,0) = d[0] - k*d[2]*rvskz;
            J(1,0) = d[1] + k*d[2]*rvckz;
            J(0,1) = rckz;
            J(1,1) = rskz;
            return std::make_pair(g, J);
        };

        bounds<Real,2> b({3,5}, {0,1});
        auto w = newton<Real,2>(f, b, vec<Real,2>(3, 0.5));
        Real t = w[0];
        Real v = w[1];
        EXPECT_LE(t, b[0].second);
        EXPECT_GE(t, b[0].first);
        EXPECT_LE(v, b[1].second);
        EXPECT_GE(v, b[1].first);
        Real u = (o[2] + t*d[2])/speed_ + 0.5;

        // v = 1, u = 1/2.
        EXPECT_FLOAT_EQ(t, 3.0);
        EXPECT_FLOAT_EQ(v, 1.0);
        EXPECT_FLOAT_EQ(u, 0.5);
        Real x1 = o[0] + t*d[0];
        Real x2 = radius_*v*cos(2*M_PI*u);
        EXPECT_FLOAT_EQ(x1, x2);

        Real y1 = o[1] + t*d[1];
        Real y2 = radius_*v*sin(2*M_PI*u);
        EXPECT_LE(abs(y1), 10*std::numeric_limits<Real>::epsilon());
        EXPECT_LE(abs(y2), 10*std::numeric_limits<Real>::epsilon());

        Real z1 = o[2] + t*d[2];
        Real z2 = speed_*(u - 1.0/2.0);
        EXPECT_FLOAT_EQ(z1, z2);

}

TEST(NewtonTest, 3DCuyt) {
    // From Cuyt, "Computational implementation of the Multivariate Halley Method"
    auto f = [](vec<Real,3> w) {
        Real x0 = w[0];
        Real x1 = w[1];
        Real x2 = w[2];
        vec<Real, 3> v;
        Real x0_sq = x0*x0;
        Real x0_4 = x0_sq*x0_sq;
        Real x1_sq = x1*x1;
        Real x1_4 = x1_sq*x1_sq;
        Real x2_sq = x2*x2;
        Real x2_4 = x2_sq*x2_sq;

        v[0] = 16*x0_4 + 16*x1_4 + x2_4 - 16;
        v[1] = x0_sq + x1_sq + x2_sq - 3;
        v[2] = x0_sq*x0 - x1;

        matrix<Real,3,3> J;
        J(0,0) = 16*4*x0_sq*x0;
        J(0,1) = 16*4*x1_sq*x1;
        J(0,2) = 4*x2_sq*x2;
        J(1,0) = 2*x0;
        J(1,1) = 2*x1;
        J(1,2) = 2*x2;
        J(2,0) = 3*x0_sq;
        J(2,1) = -1;
        J(2,2) = 0;
        return std::make_pair(v, J);
    };

    // Cuyt uses initial guesses x0 = x1 = x2 = 1.0.
    // She lists the solution as

    vec<Real,3> expected(0.877965760274, 0.676756970518, 1.33085541162);
    auto [v, M] = f(expected);
    // Her solution is indeed very accurate:
    EXPECT_LE(squared_norm(v), std::numeric_limits<Real>::epsilon());

    bounds<Real,3> b({0,2},{0,2},{0,2});
    auto sol = newton<Real,3>(f, b, b.center());
    EXPECT_FLOAT_EQ(sol[0], expected[0]);
    EXPECT_FLOAT_EQ(sol[1], expected[1]);
    EXPECT_FLOAT_EQ(sol[2], expected[2]);
}

#endif
