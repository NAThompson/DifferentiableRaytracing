#ifndef SPHERE_TEST_HPP
#define SPHERE_TEST_HPP
#include <drt/ray.hpp>
#include <drt/sphere.hpp>

using namespace drt;

TEST(SphereTest, Intersection) {
    using Real = double;
    vec<Real> center(0,0,0);
    sphere<Real> s(center, 1, nullptr);
    vec<Real> o(2,0,0);
    vec<Real> d(-1,0,0);
    ray<Real> r(o, d);
    hit_record<Real> hr;
    bool hits = s.hit(r, 0.0, 10.0, hr);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(hr.p[0], 1);
    EXPECT_FLOAT_EQ(hr.p[1], 0);
    EXPECT_FLOAT_EQ(hr.p[2], 0);
    /*
    EXPECT_FLOAT_EQ(hr.u, 0.5);
    EXPECT_FLOAT_EQ(hr.v, 0);
    EXPECT_FLOAT_EQ(hr.t, 2);
    EXPECT_FLOAT_EQ(hr.E, 4*M_PI*M_PI*hr.v*hr.v + 1);
    EXPECT_FLOAT_EQ(hr.F, 0);
    EXPECT_FLOAT_EQ(hr.G, 1);
    EXPECT_FLOAT_EQ(hr.mean_curvature(), 0);
    auto [k1, k2] = hr.principle_curvatures();
    */
}

#endif