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

    std::uniform_real_distribution<Real> dis(-1,1);
    std::mt19937_64 gen;
    int i = 0;
    while (i++ < 500) {
        v1[0] = dis(gen);
        v1[1] = dis(gen);
        v1[2] = dis(gen);

        v2[0] = dis(gen);
        v2[1] = dis(gen);
        v2[2] = dis(gen);

        auto w1 = cross(v1, v2);
        auto w2 = cross(v2, v1);

        EXPECT_FLOAT_EQ(w1[0], -w2[0]);
        EXPECT_FLOAT_EQ(w1[1], -w2[1]);
        EXPECT_FLOAT_EQ(w1[2], -w2[2]);

        w1 = cross(v1, v1);
        EXPECT_LE(abs(w1[0]), std::numeric_limits<Real>::epsilon());
        EXPECT_LE(abs(w1[1]), std::numeric_limits<Real>::epsilon());
        EXPECT_LE(abs(w1[2]), std::numeric_limits<Real>::epsilon());

        vec<Real> v3(dis(gen), dis(gen), dis(gen));
        w1 = cross(v1, v2 + v3);
        w2 = cross(v1, v2) + cross(v1, v3);
        EXPECT_FLOAT_EQ(w1[0], w2[0]);
        EXPECT_FLOAT_EQ(w1[1], w2[1]);
        EXPECT_FLOAT_EQ(w1[2], w2[2]);
    }
}

#endif