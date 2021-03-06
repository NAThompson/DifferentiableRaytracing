#ifndef TENSOR_TEST_HPP
#define TENSOR_TEST_HPP
#include <random>
#include <drt/tensor.hpp>

using namespace drt;

TEST(TensorTest, 2x2x2) {
    using Real = double;
    constexpr const int64_t d = 2;
    tensor<Real,d,d,d> T;
    for (int64_t i = 0; i < d; ++i) {
        for (int64_t j = 0; j < d; ++j) {
            for (int64_t k = 0; k < d; ++k) {
                if (i == j && j == k) {
                    T(i,j,k) = 1;
                } else {
                    T(i,j,k) = 0;
                }
            }
        }
    }

    std::mt19937_64 gen(12345);
    std::uniform_real_distribution<Real> dis(-1,1);
    vec<Real, d> x;
    vec<Real, d> y;
    for (int64_t i = 0; i < d; ++i) {
        x[i] = dis(gen);
        y[i] = dis(gen);
    }

    auto z = T(x,y);

    EXPECT_FLOAT_EQ(z[0], x[0]*y[0]);
    EXPECT_FLOAT_EQ(z[1], x[1]*y[1]);

    for (int64_t i = 0; i < d; ++i) {
        for (int64_t j = 0; j < d; ++j) {
            for (int64_t k = 0; k < d; ++k) {
                T(i,j,k) = dis(gen);
            }
        }
    }
    vec<Real, d> zero;
    zero[0] = 0;
    zero[1] = 0;
    z = T(x, zero);

    EXPECT_FLOAT_EQ(z[0], 0);
    EXPECT_FLOAT_EQ(z[1], 0);

    // Bilinearity:
    Real k1 = dis(gen);
    Real k2 = dis(gen);
    auto z1 = T(k1*x, k2*y);
    z = k1*k2*T(x,y);
    EXPECT_FLOAT_EQ(z[0], z1[0]);
    EXPECT_FLOAT_EQ(z[1], z1[1]);

    vec<Real,d> w;
    w[0] = dis(gen);
    w[1] = dis(gen);
    z = T(x, k1*y + k2*w);
    z1 = k1*T(x, y) + k2*T(x,w);

    EXPECT_FLOAT_EQ(z[0], z1[0]);
    EXPECT_FLOAT_EQ(z[1], z1[1]);
}


TEST(TensorTest, 3x3x3) {
    using Real = double;
    constexpr const int64_t d = 3;
    tensor<Real,d,d,d> T;
    for (int64_t i = 0; i < d; ++i) {
        for (int64_t j = 0; j < d; ++j) {
            for (int64_t k = 0; k < d; ++k) {
                if (i == j && j == k) {
                    T(i,j,k) = 1;
                } else {
                    T(i,j,k) = 0;
                }
            }
        }
    }

    std::mt19937_64 gen(12345);
    std::uniform_real_distribution<Real> dis(-1,1);
    vec<Real, d> x;
    vec<Real, d> y;
    for (int64_t i = 0; i < d; ++i) {
        x[i] = dis(gen);
        y[i] = dis(gen);
    }

    auto z = T(x,y);

    EXPECT_FLOAT_EQ(z[0], x[0]*y[0]);
    EXPECT_FLOAT_EQ(z[1], x[1]*y[1]);
    EXPECT_FLOAT_EQ(z[2], x[2]*y[2]);

    for (int64_t i = 0; i < d; ++i) {
        for (int64_t j = 0; j < d; ++j) {
            for (int64_t k = 0; k < d; ++k) {
                T(i,j,k) = dis(gen);
            }
        }
    }
    vec<Real, d> zero;
    zero[0] = 0;
    zero[1] = 0;
    zero[2] = 0;
    z = T(x, zero);

    EXPECT_FLOAT_EQ(z[0], 0);
    EXPECT_FLOAT_EQ(z[1], 0);
    EXPECT_FLOAT_EQ(z[2], 0);

    // Bilinearity:
    Real k1 = dis(gen);
    Real k2 = dis(gen);
    auto z1 = T(k1*x, k2*y);
    z = k1*k2*T(x,y);
    EXPECT_FLOAT_EQ(z[0], z1[0]);
    EXPECT_FLOAT_EQ(z[1], z1[1]);
    EXPECT_FLOAT_EQ(z[2], z1[2]);

    vec<Real,d> w;
    w[0] = dis(gen);
    w[1] = dis(gen);
    w[2] = dis(gen);
    z = T(x, k1*y + k2*w);
    z1 = k1*T(x, y) + k2*T(x,w);

    EXPECT_FLOAT_EQ(z[0], z1[0]);
    EXPECT_FLOAT_EQ(z[1], z1[1]);
    EXPECT_FLOAT_EQ(z[2], z1[2]);
}

// Useful for parametric surface intersection as  σ:ℝ² -> ℝ³ generates a Hessian which is a 3x2x2 tensor.
// For ray intersections, the Hessian does not have to be square since the linear part gets blasted away.
TEST(TensorTest, 3x2x2) {
    using Real = double;
    tensor<Real,3,2,2> T;
    for (int64_t i = 0; i < 3; ++i) {
        for (int64_t j = 0; j < 2; ++j) {
            for (int64_t k = 0; k < 2; ++k) {
                if (i == j && j == k) {
                    T(i,j,k) = 1;
                } else {
                    T(i,j,k) = 0;
                }
            }
        }
    }

    std::mt19937_64 gen(12345);
    std::uniform_real_distribution<Real> dis(-1,1);
    vec<Real, 2> x;
    vec<Real, 2> y;
    for (int64_t i = 0; i < 2; ++i) {
        x[i] = dis(gen);
        y[i] = dis(gen);
    }

    auto z = T(x,y);

    EXPECT_FLOAT_EQ(z[0], x[0]*y[0]);
    EXPECT_FLOAT_EQ(z[1], x[1]*y[1]);
    EXPECT_FLOAT_EQ(z[2], 0);

    for (int64_t i = 0; i < 3; ++i) {
        for (int64_t j = 0; j < 2; ++j) {
            for (int64_t k = 0; k < 2; ++k) {
                T(i,j,k) = dis(gen);
            }
        }
    }
    vec<Real, 2> zero;
    zero[0] = 0;
    zero[1] = 0;
    z = T(x, zero);

    EXPECT_FLOAT_EQ(z[0], 0);
    EXPECT_FLOAT_EQ(z[1], 0);
    EXPECT_FLOAT_EQ(z[2], 0);

    // Bilinearity:
    Real k1 = dis(gen);
    Real k2 = dis(gen);
    auto z1 = T(k1*x, k2*y);
    z = k1*k2*T(x,y);
    EXPECT_FLOAT_EQ(z[0], z1[0]);
    EXPECT_FLOAT_EQ(z[1], z1[1]);
    EXPECT_FLOAT_EQ(z[2], z1[2]);

    vec<Real,2> w;
    w[0] = dis(gen);
    w[1] = dis(gen);
    z = T(x, k1*y + k2*w);
    z1 = k1*T(x, y) + k2*T(x,w);

    EXPECT_FLOAT_EQ(z[0], z1[0]);
    EXPECT_FLOAT_EQ(z[1], z1[1]);
    EXPECT_FLOAT_EQ(z[2], z1[2]);
}
#endif
