#ifndef ROOTS_TEST_HPP
#define ROOTS_TEST_HPP
#include <array>
#include <random>
#include <drt/roots.hpp>

TEST(RootsTest, Quadratic) {
    // (x - 1)x = x^2 - x
    using Real = double;
    auto roots = drt::quadratic_roots(Real(1), Real(-1), Real(0));
    EXPECT_EQ(roots.size(), size_t(2));
    EXPECT_FLOAT_EQ(roots[0], 0);
    EXPECT_FLOAT_EQ(roots[1], 1);

    // (x-3)(x+5) = x^2 + 2x - 15
    roots = drt::quadratic_roots(Real(1), Real(2), Real(-15));
    EXPECT_EQ(roots.size(), size_t(2));
    EXPECT_FLOAT_EQ(roots[0], -5);
    EXPECT_FLOAT_EQ(roots[1], 3);

    // x^2 
    roots = drt::quadratic_roots(Real(1), Real(0), Real(0));
    EXPECT_EQ(roots.size(), size_t(2));
    EXPECT_FLOAT_EQ(roots[0], 0);
    EXPECT_FLOAT_EQ(roots[1], 0);
}

TEST(RootsTest, Cubic) {
    // x^3
    using Real = double;
    auto roots = drt::cubic_roots(Real(1), Real(0), Real(0), Real(0));
    EXPECT_EQ(roots.size(), size_t(3));
    EXPECT_FLOAT_EQ(roots[0], 0);
    EXPECT_FLOAT_EQ(roots[1], 0);
    EXPECT_FLOAT_EQ(roots[2], 0);

    // x^3 + 1. Single real root at x = -1.
    roots = drt::cubic_roots(Real(1), Real(0), Real(0), Real(1));
    EXPECT_EQ(roots.size(), size_t(1));
    EXPECT_FLOAT_EQ(roots[0], -1);

    // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6 
    roots = drt::cubic_roots(Real(1), Real(-6), Real(11), Real(-6));
    EXPECT_EQ(roots.size(), size_t(3));
    EXPECT_FLOAT_EQ(roots[0], 1);
    EXPECT_FLOAT_EQ(roots[1], 2);
    EXPECT_FLOAT_EQ(roots[2], 3);

    std::uniform_real_distribution<Real> dis(-2,2);
    std::mt19937 gen(12345);
    // Expected roots
    std::array<Real, 3> r;
    int i = 0;
    do {
        // Mathematica:
        // Expand[(x - r0)*(x - r1)*(x - r2)]
        // r0 r1 r2 r3 - (r0 r1 r2 + r0 r1 r3 + r0 r2 r3 + r1 r2 r3) x
        // + (r0 r1 + r0 r2 + r1 r2 + r0 r3 + r1 r3 + r2 r3)x^2 - (r0 + r1 + r2 + r3) x^3 + x^4
        for (auto & root : r) {
            root = dis(gen);
        }
        std::sort(r.begin(), r.end());
        Real a = 1;
        Real b = -(r[0] + r[1] + r[2]);
        Real c = r[0]*r[1] + r[0]*r[2] + r[1]*r[2];
        Real d = -r[0]*r[1]*r[2];

        auto roots = drt::cubic_roots(a, b, c, d);
        EXPECT_EQ(roots.size(), size_t(3));
        EXPECT_FLOAT_EQ(roots[0], r[0]);
        EXPECT_FLOAT_EQ(roots[1], r[1]);
        EXPECT_FLOAT_EQ(roots[2], r[2]);
    } while (i++ < 500);

}

TEST(RootsTest, Quartic) {
    // x^4
    using Real = double;
    auto roots = drt::quartic_roots(Real(1), Real(0), Real(0), Real(0), Real(0));
    EXPECT_EQ(roots.size(), size_t(4));
    EXPECT_FLOAT_EQ(roots[0], 0);
    EXPECT_FLOAT_EQ(roots[1], 0);
    EXPECT_FLOAT_EQ(roots[2], 0);
    EXPECT_FLOAT_EQ(roots[3], 0);

    // x^4 + 1. No real roots.
    roots = drt::quartic_roots(Real(1), Real(0), Real(0), Real(0), Real(1));
    EXPECT_EQ(roots.size(), size_t(0));

    // x^4 - 1. Two real roots at -1,1. 
    roots = drt::quartic_roots(Real(1), Real(0), Real(0), Real(0), Real(-1));
    EXPECT_EQ(roots.size(), size_t(2));
    EXPECT_FLOAT_EQ(roots[0], -1);
    EXPECT_FLOAT_EQ(roots[1], 1);

    // (x-1)(x-2)(x-3)(x-4) = 24 - 50 x + 35 x^2 - 10 x^3 + x^4
    roots = drt::quartic_roots(Real(1), Real(-10), Real(35), Real(-50), Real(24));
    EXPECT_EQ(roots.size(), size_t(4));
    EXPECT_FLOAT_EQ(roots[0], 1);
    EXPECT_FLOAT_EQ(roots[1], 2);
    EXPECT_FLOAT_EQ(roots[2], 3);
    EXPECT_FLOAT_EQ(roots[3], 4);

    // (x-1.5)(x+0.5)(x+1.5)(x+2) = -2.25 - 5.625 x - 1.25 x^2 + 2.5 x^3 + x^4
    roots = drt::quartic_roots(Real(1), 2.5, -1.25, -5.625, -2.25);
    EXPECT_EQ(roots.size(), size_t(4));
    EXPECT_FLOAT_EQ(roots[0], -2);
    EXPECT_FLOAT_EQ(roots[1], -1.5);
    EXPECT_FLOAT_EQ(roots[2], -0.5);
    EXPECT_FLOAT_EQ(roots[3], 1.5);

    // Test the biquadratic branch:
    // Solve[1 - 6 z^2 + z^4 == 0, z]
    // {{z -> -1 - Sqrt[2]}, {z -> 1 - Sqrt[2]}, {z -> -1 + Sqrt[2]}, {z -> 1 + Sqrt[2]}}
    roots = drt::quartic_roots(Real(1), Real(0), Real(-6), Real(0), Real(1));
    EXPECT_EQ(roots.size(), size_t(4));
    Real rt2 = std::sqrt(2.0);
    EXPECT_FLOAT_EQ(roots[0], -1 - rt2);
    EXPECT_FLOAT_EQ(roots[1], 1 - rt2);
    EXPECT_FLOAT_EQ(roots[2], -1 + rt2);
    EXPECT_FLOAT_EQ(roots[3], 1 + rt2);


    std::uniform_real_distribution<Real> dis(-2,2);
    std::mt19937 gen(12345);
    // Expected roots
    std::array<Real, 4> r;
    int i = 0;
    do {
        // Mathematica:
        // Expand[(x - r0)*(x - r1)*(x - r2)*(x - r3)]
        // r0 r1 r2 r3 - (r0 r1 r2 + r0 r1 r3 + r0 r2 r3 + r1 r2 r3) x
        // + (r0 r1 + r0 r2 + r1 r2 + r0 r3 + r1 r3 + r2 r3)x^2 - (r0 + r1 + r2 + r3) x^3 + x^4
        for (auto & root : r) {
            root = dis(gen);
        }
        std::sort(r.begin(), r.end());
        Real a = 1;
        Real b = -(r[0] + r[1] + r[2] + r[3]);
        Real c = r[0]*r[1] + r[0]*r[2] + r[1]*r[2] + r[0]*r[3] + r[1]*r[3] + r[2]*r[3];
        Real d = - (r[0]*r[1]*r[2] + r[0]*r[1]*r[3] + r[0]*r[2]*r[3] + r[1]*r[2]*r[3]);
        Real e = r[0]*r[1]*r[2]*r[3];

        auto roots = drt::quartic_roots(a, b, c, d, e);
        EXPECT_EQ(roots.size(), size_t(4));
        EXPECT_FLOAT_EQ(roots[0], r[0]);
        EXPECT_FLOAT_EQ(roots[1], r[1]);
        EXPECT_FLOAT_EQ(roots[2], r[2]);
        EXPECT_FLOAT_EQ(roots[3], r[3]);
    } while (i++ < 500);
}

#endif