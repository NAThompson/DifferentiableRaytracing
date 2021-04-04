#ifndef MATRIX_TEST_HPP
#define MATRIX_TEST_HPP
#include <random>
#include <drt/matrix.hpp>

using namespace drt;

TEST(MatrixTest, TwoByTwo) {
    using Real = double;
    matrix<Real,2,2> M;
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

    // A depressing example found in the Ko method test:
    M(0,0) = 39.4784176043574;
    M(0,1) = 8.39603594799628e-33;
    M(1,0) = M(0,1);
    M(1,1) = 9.8696044010893;
    w[0] = -6.28318530717959;
    w[1] = 3.91572165271905e-15;

    v = M.solve(w);
    auto v2 = inverse(M)*w;
    EXPECT_FLOAT_EQ(v[0], v2[0]);
    EXPECT_FLOAT_EQ(v[1], v2[1]);

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


TEST(MatrixTest, ThreeByThree) {
    using Real = double;
    // Inverse[{{1, -1, 1}, {1, 2, 1}, {2, 1, 1}}]
    // {{-(1/3), -(2/3), 1}, {-(1/3), 1/3, 0}, {1, 1, -1}}
    matrix<Real, 3, 3> M;
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

TEST(MatrixTest, 3x3Solve)
{
    using Real = double;
    matrix<Real, 3, 3> M;
    M(0,0) = 1; M(0,1) = 2; M(0,2) = 3;
    M(1,0) = 4; M(1,1) = 5; M(1,2) = 6;
    M(2,0) = 7; M(2,1) = 8; M(2,2) = 9;

    vec<Real, 3> b(2,2,2);

    vec<Real, 3> x(-2,2,0);

    auto Mx = M*x;
    EXPECT_FLOAT_EQ(Mx[0], b[0]);
    EXPECT_FLOAT_EQ(Mx[1], b[1]);
    EXPECT_FLOAT_EQ(Mx[2], b[2]);

    auto x1 = M.solve(b);

    EXPECT_FLOAT_EQ(x1[0], x[0]);
    EXPECT_FLOAT_EQ(x1[1], x[1]);
    EXPECT_FLOAT_EQ(x1[2], x[2]);


    std::uniform_real_distribution<Real> dis(-1,1);
    std::mt19937_64 gen;
    int64_t n = 0;
    while (n++ < 512) {
        for (int64_t i = 0; i < 2; ++i) {
            for (int64_t j = 0; j < 2; ++j) {
                M(i,j) = dis(gen);
            }
        }
        vec<Real, 3> b(dis(gen), dis(gen), dis(gen));

        auto x = M.solve(b);
        auto Mx = M*x;
        EXPECT_FLOAT_EQ(Mx[0], b[0]);
        EXPECT_FLOAT_EQ(Mx[1], b[1]);
        EXPECT_FLOAT_EQ(Mx[2], b[2]);

        M(0,0) = 0;
        x = M.solve(b);
        Mx = M*x;
        EXPECT_FLOAT_EQ(Mx[0], b[0]);
        EXPECT_FLOAT_EQ(Mx[1], b[1]);
        EXPECT_FLOAT_EQ(Mx[2], b[2]);

    }
}
#endif
