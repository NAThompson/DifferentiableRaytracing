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

    auto uvt = drt::ko_method(s, r, 0.0, 0.5, bound);
    EXPECT_FLOAT_EQ(uvt[0], 0.0);
    EXPECT_FLOAT_EQ(uvt[1], 0.5);
    EXPECT_FLOAT_EQ(uvt[2], 1.0);

    //uvt = drt::ko_method(s, r, 0.0, 0.25, bound);
    //uvt = drt::ko_method(s, r, 0.0, 0.483019009238001, bound);
    //uvt = drt::ko_method(s, r, 0.0, 0.499991963816114, bound);
    uvt = drt::ko_method(s, r, 0.0, 0.499999999999999, bound);

}

#endif
