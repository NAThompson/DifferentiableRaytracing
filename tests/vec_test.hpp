#ifndef VEC_TEST_HPP
#define VEC_TEST_HPP
#include <drt/vec.hpp>

using namespace drt;

TEST(VecTest, Cross) {
    using Real = double;
    vec<Real> v1(1,0,0);
    vec<Real> v2(1,0,0);
    vec<Real> v3 = cross(v1, v2);

    Real length = norm(v3);
    EXPECT_FLOAT_EQ(length, 0);

    v2[0] = 0;
    v2[1] = 1;
    v3 = cross(v1, v2);
    EXPECT_FLOAT_EQ(v3[0], 0);
    EXPECT_FLOAT_EQ(v3[1], 0);
    EXPECT_FLOAT_EQ(v3[2], 1);

    v2[0] = 0;
    v2[1] = 0;
    v2[2] = 1;
    v3 = cross(v1, v2);
    EXPECT_FLOAT_EQ(v3[0], 0);
    EXPECT_FLOAT_EQ(v3[1], -1);
    EXPECT_FLOAT_EQ(v3[2], 0);
}

#endif