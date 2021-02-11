#ifndef HELICOID_TEST_HPP
#define HELICOID_TEST_HPP
#include <drt/ray.hpp>
#include <drt/helicoid.hpp>

using namespace drt;

TEST(HelicoidTest, Intersection) {
    using Real = double;
    auto mat = std::make_shared<drt::lambertian<Real>>(vec<Real>(1,0,0));
    helicoid<Real> h(1, 1, mat);
    vec<Real> o(2,0,0);
    vec<Real> d(-1,0,0);
    ray<Real> r(o, d);
    hit_record<Real> hr;
    bool hits = h.hit(r, 0.0, 10.0, hr);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(hr.p[0], 1);
    EXPECT_FLOAT_EQ(hr.p[1], 0);
    EXPECT_FLOAT_EQ(hr.p[2], 0);
    EXPECT_FLOAT_EQ(hr.u, 0.5);
    EXPECT_FLOAT_EQ(hr.v, 1);
    EXPECT_FLOAT_EQ(hr.t, 1);

    std::cerr << hr << "\n";
}

#endif