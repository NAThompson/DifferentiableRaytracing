#ifndef INTERVAL_TEST_HPP
#define INTERVAL_TEST_HPP
#include <drt/interval.hpp>

using namespace drt;

TEST(IntervalTest, Multiplication) {
    using Real = double;
    auto x = interval<Real>(2.2, 2.4);
    std::cout << "x = " << x << "\n";

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

    z = x*y;
    EXPECT_EQ(z.lower(), 8);
    EXPECT_EQ(z.upper(), 15);
}

#endif