#ifndef MAT_TEST_HPP
#define MAT_TEST_HPP
#include <drt/mat.hpp>

using namespace drt;

TEST(MatTest, TwoByTwo) {
    using Real = double;
    mat<Real,2,2> M;
    M(0,0) = 1;
    M(0,1) = -1;
    M(1,0) = 2;
    M(1,1) = 3;
    EXPECT_FLOAT_EQ(determinant(M), 5);

    vec<Real, 2> v;
    v[0] = 3;
    v[1] = 8;
    auto w = M*v;
    EXPECT_FLOAT_EQ(w[0], -5);
    EXPECT_FLOAT_EQ(w[1], 30);

    auto M_inv = inverse(M);
    EXPECT_FLOAT_EQ(M_inv(0,0), Real(3)/5);
    EXPECT_FLOAT_EQ(M_inv(0,1), Real(1)/5);
    EXPECT_FLOAT_EQ(M_inv(1,0), -Real(2)/5);
    EXPECT_FLOAT_EQ(M_inv(1,1), Real(1)/5);

    v = M_inv*w;
    EXPECT_FLOAT_EQ(v[0], 3);
    EXPECT_FLOAT_EQ(v[1], 8);

    v = M.solve(w);
    EXPECT_FLOAT_EQ(v[0], 3);
    EXPECT_FLOAT_EQ(v[1], 8);

    std::mt19937_64 gen(12345);
    std::uniform_real_distribution<Real> dis;

    int i = 0;
    while (i++ < 512) {
        vec<Real, 2> b;
        b[0] = dis(gen);
        b[1] = dis(gen);
        M(0,0) = dis(gen);
        M(0,1) = dis(gen);
        M(1,0) = dis(gen);
        M(1,1) = dis(gen);
        // Mx = b
        auto x1 = M.solve(b);
        auto M_inv = inverse(M);
        auto x2 = M_inv*b;
        EXPECT_FLOAT_EQ(x1[0], x2[0]);
        EXPECT_FLOAT_EQ(x1[1], x2[1]);

        M(0,0) = 0;
        M(0,1) = dis(gen);
        M(1,0) = dis(gen);
        M(1,1) = dis(gen);
        // Mx = b
        x1 = M.solve(b);
        M_inv = inverse(M);
        x2 = M_inv*b;
        EXPECT_FLOAT_EQ(x1[0], x2[0]);
        EXPECT_FLOAT_EQ(x1[1], x2[1]);

        M(0,0) = dis(gen);
        M(0,1) = dis(gen);
        M(1,0) = dis(gen);
        M(1,1) = 0;
        // Mx = b
        x1 = M.solve(b);
        M_inv = inverse(M);
        x2 = M_inv*b;
        EXPECT_FLOAT_EQ(x1[0], x2[0]);
        EXPECT_FLOAT_EQ(x1[1], x2[1]);

        M(0,0) = dis(gen);
        M(0,1) = 0;
        M(1,0) = dis(gen);
        M(1,1) = dis(gen);
        // Mx = b
        x1 = M.solve(b);
        M_inv = inverse(M);
        x2 = M_inv*b;
        EXPECT_FLOAT_EQ(x1[0], x2[0]);
        EXPECT_FLOAT_EQ(x1[1], x2[1]);

        M(0,0) = dis(gen);
        M(0,1) = dis(gen);
        M(1,0) = dis(gen);
        M(1,1) = 0;
        // Mx = b
        x1 = M.solve(b);
        M_inv = inverse(M);
        x2 = M_inv*b;
        EXPECT_FLOAT_EQ(x1[0], x2[0]);
        EXPECT_FLOAT_EQ(x1[1], x2[1]);

        M(0,0) = dis(gen);
        M(0,1) = 0;
        M(1,0) = 0;
        M(1,1) = dis(gen);
        // Mx = b
        x1 = M.solve(b);
        M_inv = inverse(M);
        x2 = M_inv*b;
        EXPECT_FLOAT_EQ(x1[0], x2[0]);
        EXPECT_FLOAT_EQ(x1[1], x2[1]);

    }
}


TEST(MatTest, ThreeByThree) {
    using Real = double;
    // Inverse[{{1, -1, 1}, {1, 2, 1}, {2, 1, 1}}]
    // {{-(1/3), -(2/3), 1}, {-(1/3), 1/3, 0}, {1, 1, -1}}
    mat<Real, 3, 3> M;
    M(0,0) = 1;
    M(0,1) = -1;
    M(0,2) = 1;
    M(1,0) = 1;
    M(1,1) = 2;
    M(1,2) = 1;
    M(2,0) = 2;
    M(2,1) = 1;
    M(2,2) = 1;
    EXPECT_FLOAT_EQ(determinant(M), -3);

    vec<Real, 3> v(1,2,3);
    auto w = M*v;
    EXPECT_FLOAT_EQ(w[0], 2);
    EXPECT_FLOAT_EQ(w[1], 8);
    EXPECT_FLOAT_EQ(w[2], 7);
    auto M_inv = inverse(M);
    Real t = 1.0/3.0;
    EXPECT_FLOAT_EQ(M_inv(0,0), -t);
    EXPECT_FLOAT_EQ(M_inv(0,1), -2*t);
    EXPECT_FLOAT_EQ(M_inv(0,2), 1);

    EXPECT_FLOAT_EQ(M_inv(1,0), -t);
    EXPECT_FLOAT_EQ(M_inv(1,1), t);
    EXPECT_FLOAT_EQ(M_inv(1,2), 0);

    EXPECT_FLOAT_EQ(M_inv(2,0), 1);
    EXPECT_FLOAT_EQ(M_inv(2,1), 1);
    EXPECT_FLOAT_EQ(M_inv(2,2), -1);
    v = M_inv*w;
    EXPECT_FLOAT_EQ(v[0], 1);
    EXPECT_FLOAT_EQ(v[1], 2);
    EXPECT_FLOAT_EQ(v[2], 3);
}
#endif