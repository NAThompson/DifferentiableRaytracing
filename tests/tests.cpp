#include <iostream>
#include <cmath>
#include <cstring>
#include <gtest/gtest.h>
#include <drt/roots.hpp>
#include <drt/torus.hpp>
#include <drt/ray.hpp>
#include <drt/lambertian.hpp>

using namespace drt;

// See: https://randomascii.wordpress.com/2012/01/23/stupid-float-tricks-2/ for why this works.
template<typename Real>
inline uint64_t float_distance(Real x, Real y)
{
    static_assert(std::numeric_limits<Real>::has_denorm == std::denorm_present,
                "float_distance presumes the floating-point type has subnormal numbers.");

    if (!std::isfinite(x) || !std::isfinite(y)) {
        return 0xFFFFFFFFFFFFFFFFL;
    }

    // Signed zero is the sworn enemy of this process.
    if (y == 0) {
        y = std::abs(y);
    }
    if (x == 0) {
        x = std::abs(x);
    }

    if ( (x < 0 && y >= 0) || (x >= 0 && y < 0) )
    {
        uint64_t dx, dy;
        if (x < 0)
        {
            dy = float_distance(Real(0), y);
            dx = float_distance(Real(0), -x);
        }
        else
        {
            dy = float_distance(Real(0), -y);
            dx = float_distance(Real(0), x);
        }
        return dx + dy;
    }

    if (x < 0 && y < 0)
    {
        return float_distance(-x, -y);
    }

    // Note that:
    // int64_t xi = *reinterpret_cast<int64_t*>(&x);
    // int64_t yi = *reinterpret_cast<int64_t*>(&y);
    // also works, but generates warnings.
    // Good option to have if we get compile errors off memcpy or don't want to #include <cstring> though.
    // At least on gcc, both versions generate the same assembly.
    uint64_t xi, yi;
    if constexpr (std::is_same_v<Real, double>)
    {
        memcpy(&xi, &x, sizeof(uint64_t));
        memcpy(&yi, &y, sizeof(uint64_t));
    }
    else if constexpr (std::is_same_v<Real, float>)
    {
        uint32_t xi_32, yi_32;
        memcpy(&xi_32, &x, sizeof(uint32_t));
        memcpy(&yi_32, &y, sizeof(uint32_t));
        xi = xi_32;
        yi = yi_32;
    }
    else
    {
        std::cerr << "Type is neither float nor double; this is unsupported.\n";
        return 0xFFFFFFFFFFFFFFFFL;
    }
    if (yi > xi) {
        return yi - xi;
    }
    return xi - yi;
}


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
}

TEST(Torus, Intersection)
{
    using Real = double;
    auto mat = std::make_shared<lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> origin(0,0,-5);
    vec<Real> direction(0,0,1);
    vec<Real> center(0,0,0);
    auto r = ray<Real>(origin, direction);
    auto tor = torus<Real>(center, 3.0, 1.0, mat);

    hit_record<Real> rec;
    bool hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_FALSE(hits);

    // We should get two solutions here:
    origin[0] = 3;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(rec.p[0], 3);
    EXPECT_FLOAT_EQ(rec.p[1], 0);
    EXPECT_FLOAT_EQ(rec.p[2], -1);
    EXPECT_FLOAT_EQ(rec.t, 4);

    EXPECT_FLOAT_EQ(rec.normal[0], 0);
    EXPECT_FLOAT_EQ(rec.normal[1], 0);
    EXPECT_FLOAT_EQ(rec.normal[2], -1);

    // We should get no solutions here:
    origin[0] = 5;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_FALSE(hits);

    origin[0] = -5;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_FALSE(hits);

}

TEST(Torus, AABB)
{
    using Real = double;
    auto mat = std::make_shared<lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> center(0,0,0);
    auto tor = torus<Real>(center, 3.0, 1.0, mat);
    aabb<Real> box;
    EXPECT_TRUE(tor.bounding_box(box));
    EXPECT_EQ(box.min_[0], -4);
    EXPECT_EQ(box.min_[1], -4);
    EXPECT_EQ(box.min_[2], -1);

    EXPECT_EQ(box.max_[0], 4);
    EXPECT_EQ(box.max_[1], 4);
    EXPECT_EQ(box.max_[2], 1);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
