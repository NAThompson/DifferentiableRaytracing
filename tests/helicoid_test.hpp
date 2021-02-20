#ifndef HELICOID_TEST_HPP
#define HELICOID_TEST_HPP
#include <drt/ray.hpp>
#include <drt/helicoid.hpp>

using namespace drt;
using std::sqrt;
using std::sin;
using std::cos;

TEST(HelicoidTest, Intersection) {
    using Real = double;
    helicoid<Real> h(1, 1);
    vec<Real> o(2,0,0);
    vec<Real> d(-1,0,0);
    ray<Real> r(o, d);
    drt::hit_record<Real> hr;
    // This is a somewhat pathological intersection, since the line intersects at uncountably many points.
    bool hits = h.hit(r, 0.0, 10.0, hr);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(hr.p[0], 0);
    EXPECT_FLOAT_EQ(hr.p[1], 0);
    EXPECT_FLOAT_EQ(hr.p[2], 0);
    EXPECT_FLOAT_EQ(hr.u, 0.5);
    EXPECT_FLOAT_EQ(hr.v, 0);
    EXPECT_FLOAT_EQ(hr.t, 2);
    // E = 4π^2r^2v^2 + λ^2
    EXPECT_FLOAT_EQ(hr.E, 4*M_PI*M_PI*hr.v*hr.v + 1);
    EXPECT_FLOAT_EQ(hr.F, 0);
    EXPECT_FLOAT_EQ(hr.G, 4*M_PI*M_PI);
    EXPECT_FLOAT_EQ(hr.mean_curvature(), 0);
    auto [k1, k2] = hr.principal_curvatures();
    EXPECT_FLOAT_EQ(k1, -k2);

    o[0] = -2;
    d[0] = 1;
    r = ray<Real>(o, d);
    // Again a pathological intersection.
    hits = h.hit(r, 0.0, 10.0, hr);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(hr.p[0], -1);
    EXPECT_FLOAT_EQ(hr.p[1], 0);
    EXPECT_FLOAT_EQ(hr.p[2], 0);
    EXPECT_FLOAT_EQ(hr.u, 0.5);
    EXPECT_FLOAT_EQ(hr.v, 1);
    EXPECT_FLOAT_EQ(hr.t, 1);
    EXPECT_FLOAT_EQ(hr.E, 4*M_PI*M_PI*hr.v*hr.v + 1);
    EXPECT_FLOAT_EQ(hr.F, 0);
    EXPECT_FLOAT_EQ(hr.G, 4*M_PI*M_PI);
    EXPECT_FLOAT_EQ(hr.mean_curvature(), 0);
    std::tie(k1, k2) = hr.principal_curvatures();
    EXPECT_FLOAT_EQ(k1, -k2);

    // Looking down from +z.
    // v = 0, u = 1. Pathological since only constraint is 3/2 = u+t; minimal t occurs for maximal u.
    // Of course there are lots of pathological cases since the helicoid is a ruled surface.
    d[0] = 0;
    d[1] = 0;
    d[2] = -1;
    o[0] = 0;
    o[1] = 0;
    o[2] = 1;
    r = ray<Real>(o, d);
    hits = h.hit(r, 0.0, 10.0, hr);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(hr.p[0], 0);
    EXPECT_FLOAT_EQ(hr.p[1], 0);
    EXPECT_FLOAT_EQ(hr.p[2], 0.5);
    EXPECT_FLOAT_EQ(hr.u, 1);
    EXPECT_FLOAT_EQ(hr.v, 0);
    EXPECT_FLOAT_EQ(hr.t, 0.5);
    EXPECT_FLOAT_EQ(hr.E, 4*M_PI*M_PI*hr.v*hr.v + 1);
    EXPECT_FLOAT_EQ(hr.F, 0);
    EXPECT_FLOAT_EQ(hr.G, 4*M_PI*M_PI);
    EXPECT_FLOAT_EQ(hr.mean_curvature(), 0);
    std::tie(k1, k2) = hr.principal_curvatures();
    EXPECT_FLOAT_EQ(k1, -k2);


    std::uniform_real_distribution<Real> dist(0, 1);
    std::mt19937_64 gen;
    int i = 0;
    while (i++ < 512) {
        d[0] = 0;
        d[1] = 0;
        d[2] = -1;
        Real u = dist(gen);
        Real v = dist(gen);
        o[0] = v*std::cos(2*M_PI*u);
        o[1] = v*std::sin(2*M_PI*u);
        o[2] = 1;
        r = ray<Real>(o, d);
        hits = h.hit(r, 0.0, 10.0, hr);
        EXPECT_TRUE(hits);
        EXPECT_FLOAT_EQ(hr.u, u);
        EXPECT_FLOAT_EQ(hr.v, v);
        EXPECT_FLOAT_EQ(hr.p[0], o[0]);
        EXPECT_FLOAT_EQ(hr.p[1], o[1]);
        EXPECT_FLOAT_EQ(hr.p[2], u - 0.5);
        EXPECT_FLOAT_EQ(hr.t, 1.5 - u);
        EXPECT_FLOAT_EQ(squared_norm(hr.normal), 1);
        Real rdenom = 1.0/sqrt(1 + 4*M_PI*M_PI*v*v);
        EXPECT_FLOAT_EQ(hr.normal[0], sin(2*M_PI*u)*rdenom);
        EXPECT_FLOAT_EQ(hr.normal[1], -cos(2*M_PI*u)*rdenom);
        EXPECT_FLOAT_EQ(hr.normal[2], 2*M_PI*v*rdenom);

        EXPECT_FLOAT_EQ(hr.E, 4*M_PI*M_PI*v*v + 1);
        EXPECT_LE(abs(hr.F), 10*std::numeric_limits<Real>::epsilon());
        EXPECT_FLOAT_EQ(hr.G, 4*M_PI*M_PI);
        EXPECT_FLOAT_EQ(hr.e, 0);
        // This differs in sign from Mathworld; no problem since the normal is defined only up to sign.
        EXPECT_FLOAT_EQ(hr.f, 2*M_PI/sqrt(1 + 4*M_PI*M_PI*v*v));
        EXPECT_FLOAT_EQ(hr.g, 0);

        EXPECT_FLOAT_EQ(hr.mean_curvature(), 0);
        std::tie(k1, k2) = hr.principal_curvatures();
        EXPECT_FLOAT_EQ(k1, -k2);

        Real gc_d = 1 + 4*M_PI*M_PI*v*v;
        EXPECT_FLOAT_EQ(hr.gaussian_curvature(), -1.0/(gc_d*gc_d));
        EXPECT_FLOAT_EQ(k1*k2, -1.0/(gc_d*gc_d));
    }

}

#endif
