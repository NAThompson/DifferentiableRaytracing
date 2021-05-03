#ifndef INTERVAL_TEST_HPP
#define INTERVAL_TEST_HPP
#include <cmath>
#include <drt/interval.hpp>

using namespace drt;

TEST(IntervalTest, Multiplication) {
    using Real = double;
    auto x = interval<Real>(2.2, 2.4);
    auto y = interval<Real>(2.4, 2.6);

    EXPECT_GE(y, x);
    // Reflexivity fails for non-tight intervals: y >= y and x >= x fail,
    // But for a tight interval, it should succeed:
    auto z = interval<Real>(2.3, 2.3);
    EXPECT_GE(z,z);

    x = interval<double>(2, 3);
    y = interval<double>(4, 5);

    z = x + y;
    EXPECT_EQ(z.lower(), 6);
    EXPECT_EQ(z.upper(), 8);

    // [2,3]x[4,5] = [8,15]
    z = x*y;
    EXPECT_EQ(z.lower(), 8);
    EXPECT_EQ(z.upper(), 15);

    auto z1 = x/y;
    auto z2 = x*interval<Real>(Real(1)/5, Real(1)/4);
    EXPECT_EQ(z1.lower(), z2.lower());
    EXPECT_EQ(z1.upper(), z2.upper());

    x = interval<Real>(1,2);
    y = interval<Real>(3*std::numeric_limits<Real>::epsilon()/4, 7*std::numeric_limits<Real>::epsilon()/8);
    z = x + y;
    // We need to round down for the lower:
    EXPECT_EQ(z.lower(), 1);
    // We need to round up for the upper:
    EXPECT_EQ(z.upper(), 2*(1 + std::numeric_limits<Real>::epsilon()));

    // Set intersection:
    x = interval<Real>(1,2);
    y = interval<Real>(3,4);
    z = intersection(x,y);
    EXPECT_TRUE(z.empty());

    x = interval<Real>(1,3);
    y = interval<Real>(2,4);
    z = intersection(x,y);
    EXPECT_EQ(z.lower(), 2);
    EXPECT_EQ(z.upper(), 3);

    z = intersection(y,x);
    EXPECT_EQ(z.lower(), 2);
    EXPECT_EQ(z.upper(), 3);

    x = interval<Real>(-1,1);
    y = interval<Real>(0,2);
    z = intersection(x,y);
    EXPECT_EQ(z.lower(), 0);
    EXPECT_EQ(z.upper(), 1);

    z = intersection(y,x);
    EXPECT_EQ(z.lower(), 0);
    EXPECT_EQ(z.upper(), 1);

    // Subtraction:
    x = interval<Real>(1,2);
    y = interval<Real>(0,2);
    z = x - y;
    EXPECT_EQ(z.lower(), -1);
    EXPECT_EQ(z.upper(), 2);

}

TEST(IntervalTest, BrouwerFixPoint) {

    auto f = [](interval<Real> x) {
        if (x.lower() < 1 || x.upper() > 2) {
            std::cerr << "c'mon.\n";
        }
        const int originalRounding = fegetround();
        std::fesetround(FE_UPWARD);
        Real f1 = x.lower()/2 + 1/x.lower();
        Real f2 = x.upper()/2 + 1/x.upper();
        std::fesetround(originalRounding);
        Real f = std::min(f1,f2);
        interval<Real> y = x/interval<Real>(2,2) + interval<Real>(1,1)/x;
        if (y.upper() > f) {
            y.b_ = f;
        }
        return y;
    };

    interval<Real> x0(1,2);

    while (x0.width() > std::numeric_limits<Real>::epsilon()) {
        interval<Real> x1 = f(x0);
        x0 = x1;
    }
    // The upper bound is 1 bit too low!!
    //std::cout << std::hexfloat;
    //std::cout << x0 << "\n";
    //std::cout << "[" << sqrt(2) << ", " << sqrt(2) << "]\n";
    EXPECT_LE(x0.lower(), sqrt(2));
}

#endif
