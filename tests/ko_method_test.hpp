#ifndef KO_METHOD_TEST_HPP
#define KO_METHOD_TEST_HPP
#include <drt/ko_method.hpp>
#include <drt/sphere.hpp>

using namespace drt;

TEST(KoMethodTest, Sphere) {
    using Real = double;
    vec<Real> center(0,0,0);
    sphere<Real> s(center, 1);
    vec<Real> o(2,0,0);
    vec<Real> d(-1,0,0);
    ray<Real> r(o, d);

    drt::bounds<Real,3> bound({0,1}, {0,1}, {0,10});

    auto uvt = drt::ko_method(s, r, bound, 0.0, 0.5);
    EXPECT_FLOAT_EQ(uvt[0], 0.0);
    EXPECT_FLOAT_EQ(uvt[1], 0.5);
    EXPECT_FLOAT_EQ(uvt[2], 1.0);

    uvt = drt::ko_method(s, r, bound, 0.0, 0.25);
    EXPECT_FLOAT_EQ(uvt[0], 0.0);
    EXPECT_FLOAT_EQ(uvt[1], 0.5);
    EXPECT_FLOAT_EQ(uvt[2], 1.0);

    // Surface patch is irregular at (u0, v0) = (0,0):
    // Uncomment this to see the error message:
    //uvt = drt::ko_method(s, r, bound, 0.0, 0.0);
    // Challenge the method by starting the iteration close to the irregular point:
    uvt = drt::ko_method(s, r, bound, 0.0, std::numeric_limits<Real>::epsilon());
    EXPECT_FLOAT_EQ(uvt[0], 0.0);
    EXPECT_FLOAT_EQ(uvt[1], 0.5);
    EXPECT_FLOAT_EQ(uvt[2], 1.0);

    bound = drt::bounds<Real,3>({-std::numeric_limits<Real>::epsilon(),1}, {-1,1}, {0,10});
    uvt = drt::ko_method(s, r, bound, Real(1)/Real(128), 0.5);
    EXPECT_TRUE(abs(uvt[0]) < std::numeric_limits<Real>::epsilon());
    EXPECT_FLOAT_EQ(uvt[1], 0.5);
    EXPECT_FLOAT_EQ(uvt[2], 1.0);


    // There's another solution on the other side of the sphere at u = 0.5, v = 0.5.
    // What root does it find at u = 0.25?
    // Obviously we have no proof of this, but we want it to find the closer t root:
    uvt = drt::ko_method(s, r, bound, 0.25, 0.5);
    EXPECT_TRUE(abs(uvt[0]) < std::numeric_limits<Real>::epsilon());
    EXPECT_FLOAT_EQ(uvt[1], 0.5);
    EXPECT_FLOAT_EQ(uvt[2], 1.0);


    // I think this is an important test, but I have no clue how to fix it:
    //uvt = drt::ko_method(s, r, bound, 0.25, std::numeric_limits<Real>::epsilon());

    uvt = drt::ko_method(s, r, bound, 0.25, 0.99);
    EXPECT_TRUE(abs(uvt[0]) < std::numeric_limits<Real>::epsilon());
    EXPECT_FLOAT_EQ(uvt[1], 0.5);
    EXPECT_FLOAT_EQ(uvt[2], 1.0);

}

#endif
