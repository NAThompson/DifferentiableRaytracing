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

    auto uvt = drt::ko_method(s, r, 0.5, 0.5);
}

#endif
