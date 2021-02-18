#ifndef ELLIPSOID_TEST_HPP
#define ELLIPSOID_TEST_HPP
#include <drt/ray.hpp>
#include <drt/ellipsoid.hpp>

using namespace drt;

TEST(EllipsoidTest, Intersection) {
    using Real = double;
    vec<Real> center(0,0,0);
    vec<Real> scales(1,1,1);
    ellipsoid<Real> s(center, scales);
    vec<Real> o(2,0,0);
    vec<Real> d(-1,0,0);
    ray<Real> r(o, d);
    drt::hit_record<Real> hr;
    bool hits = s.hit(r, 0.0, 10.0, hr);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(hr.p[0], 1);
    EXPECT_FLOAT_EQ(hr.p[1], 0);
    EXPECT_FLOAT_EQ(hr.p[2], 0);
    EXPECT_FLOAT_EQ(hr.t, 1);
    EXPECT_FLOAT_EQ(hr.normal[0], 1);
    EXPECT_FLOAT_EQ(hr.normal[1], 0);
    EXPECT_FLOAT_EQ(hr.normal[2], 0);
    EXPECT_FLOAT_EQ(hr.u, 0.0);
    EXPECT_FLOAT_EQ(hr.v, 0.5);
    EXPECT_FLOAT_EQ(hr.F, 0);
    EXPECT_FLOAT_EQ(hr.mean_curvature(), 1);
    EXPECT_FLOAT_EQ(hr.gaussian_curvature(), 1);
    auto [k1, k2] = hr.principal_curvatures();
    EXPECT_FLOAT_EQ(k1, 1);
    EXPECT_FLOAT_EQ(k2, 1);


    scales = vec<Real>(2,2,2);
    s = ellipsoid<Real>(center, scales);
    o = vec<Real>(4,0,0);
    r = ray<Real>(o, d);
    hits = s.hit(r, 0.0, 10.0, hr);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(hr.p[0], 2);
    EXPECT_FLOAT_EQ(hr.p[1], 0);
    EXPECT_FLOAT_EQ(hr.p[2], 0);
    EXPECT_FLOAT_EQ(hr.t, 2);
    EXPECT_FLOAT_EQ(hr.normal[0], 1);
    EXPECT_FLOAT_EQ(hr.normal[1], 0);
    EXPECT_FLOAT_EQ(hr.normal[2], 0);
    EXPECT_FLOAT_EQ(hr.u, 0.0);
    EXPECT_FLOAT_EQ(hr.v, 0.5);
    EXPECT_FLOAT_EQ(hr.F, 0);
    // H = 1/R for a sphere
    EXPECT_FLOAT_EQ(hr.mean_curvature(), 0.5);
    // K = 1/R^2 for a sphere.
    EXPECT_FLOAT_EQ(hr.gaussian_curvature(), 0.25);
    // k1 = k2 = H on a sphere.
    std::tie(k1, k2) = hr.principal_curvatures();
    EXPECT_FLOAT_EQ(k1, 0.5);
    EXPECT_FLOAT_EQ(k2, 0.5);
}

#endif
