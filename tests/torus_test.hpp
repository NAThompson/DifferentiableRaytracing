#ifndef TORUS_TEST_HPP
#define TORUS_TEST_HPP

#include <drt/torus.hpp>
#include <drt/ray.hpp>
#include <drt/lambertian.hpp>

using drt::vec;
using drt::ray;
using drt::torus;
using drt::lambertian;
using drt::aabb;

TEST(Torus, Intersection)
{
    using Real = double;
    auto mat = std::make_shared<drt::lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> origin(0,0,-5);
    vec<Real> direction(0,0,1);
    vec<Real> center(0,0,0);
    auto r = ray<Real>(origin, direction);
    auto tor = torus<Real>(center, 3.0, 1.0, mat);

    drt::hit_record<Real> rec;
    bool hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_FALSE(hits);

    // We should get two solutions here:
    origin[0] = 3;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(rec.p[0], 3);
    EXPECT_FLOAT_EQ(rec.p[1], 0);
    EXPECT_FLOAT_EQ(rec.p[2], -1);
    EXPECT_FLOAT_EQ(rec.t, 4);
    Real res = tor.residual(rec.p);
    Real ex_res = tor.expected_residual(rec.p);
    EXPECT_LE(abs(res), ex_res);
    EXPECT_FLOAT_EQ(rec.t, 4);
    EXPECT_FLOAT_EQ(rec.normal[0], 0);
    EXPECT_FLOAT_EQ(rec.normal[1], 0);
    EXPECT_FLOAT_EQ(rec.normal[2], -1);
    EXPECT_FLOAT_EQ(rec.v, 0.75);
    EXPECT_FLOAT_EQ(rec.u, 0.0);
    EXPECT_LE(abs(tor.gaussian_curvature(rec.p)), std::numeric_limits<Real>::epsilon());
    EXPECT_LE(abs(rec.gaussian_curvature()), std::numeric_limits<Real>::epsilon());
    EXPECT_FLOAT_EQ(rec.mean_curvature(), -0.5);

    origin[0] = -3;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(rec.p[0], -3);
    EXPECT_FLOAT_EQ(rec.p[1], 0);
    EXPECT_FLOAT_EQ(rec.p[2], -1);
    res = tor.residual(rec.p);
    ex_res = tor.expected_residual(rec.p);
    EXPECT_LE(abs(res), ex_res);
    EXPECT_FLOAT_EQ(rec.t, 4);
    EXPECT_FLOAT_EQ(rec.normal[0], 0);
    EXPECT_FLOAT_EQ(rec.normal[1], 0);
    EXPECT_FLOAT_EQ(rec.normal[2], -1);
    EXPECT_LE(abs(tor.gaussian_curvature(rec.p)), std::numeric_limits<Real>::epsilon());
    EXPECT_LE(abs(rec.gaussian_curvature()), std::numeric_limits<Real>::epsilon());

    origin[0] = 0;
    origin[1] = 3;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(rec.p[0], 0);
    EXPECT_FLOAT_EQ(rec.p[1], 3);
    EXPECT_FLOAT_EQ(rec.p[2], -1);
    res = tor.residual(rec.p);
    ex_res = tor.expected_residual(rec.p);
    EXPECT_LE(abs(res), ex_res);
    EXPECT_FLOAT_EQ(rec.t, 4);
    EXPECT_FLOAT_EQ(rec.normal[0], 0);
    EXPECT_FLOAT_EQ(rec.normal[1], 0);
    EXPECT_FLOAT_EQ(rec.normal[2], -1);

    origin[0] = 0;
    origin[1] = -3;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(rec.p[0], 0);
    EXPECT_FLOAT_EQ(rec.p[1], -3);
    EXPECT_FLOAT_EQ(rec.p[2], -1);
    res = tor.residual(rec.p);
    ex_res = tor.expected_residual(rec.p);
    EXPECT_LE(abs(res), ex_res);
    EXPECT_FLOAT_EQ(rec.t, 4);
    EXPECT_FLOAT_EQ(rec.normal[0], 0);
    EXPECT_FLOAT_EQ(rec.normal[1], 0);
    EXPECT_FLOAT_EQ(rec.normal[2], -1);

    // We should get no solutions here:
    origin[0] = 5;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_FALSE(hits);

    origin[0] = -5;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_FALSE(hits);


    origin[0] = 0;
    origin[1] = -5;
    origin[2] = 0;
    direction[0] = 0;
    direction[1] = 1;
    direction[2] = 0;
    r = ray<Real>(origin, direction);
    hits = tor.hit(r, std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::infinity(), rec);
    EXPECT_TRUE(hits);
    EXPECT_FLOAT_EQ(rec.p[0], 0);
    EXPECT_FLOAT_EQ(rec.p[1], -4);
    EXPECT_FLOAT_EQ(rec.p[2], 0);
    res = tor.residual(rec.p);
    ex_res = tor.expected_residual(rec.p);
    EXPECT_LE(abs(res), ex_res);
    EXPECT_FLOAT_EQ(rec.t, 1);
    EXPECT_FLOAT_EQ(rec.normal[0], 0);
    EXPECT_FLOAT_EQ(rec.normal[1], -1);
    EXPECT_FLOAT_EQ(rec.normal[2], 0);
    EXPECT_FLOAT_EQ(rec.v, 0);
    EXPECT_FLOAT_EQ(rec.u, 0.75);

}

TEST(Torus, AABB)
{
    using Real = double;
    auto mat = std::make_shared<lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> center(0,0,0);
    auto tor = torus<Real>(center, 3.0, 1.0, mat);
    aabb<Real> box;
    EXPECT_TRUE(tor.bounding_box(box));
    EXPECT_EQ(box.min_[0], -4);
    EXPECT_EQ(box.min_[1], -4);
    EXPECT_EQ(box.min_[2], -1);

    EXPECT_EQ(box.max_[0], 4);
    EXPECT_EQ(box.max_[1], 4);
    EXPECT_EQ(box.max_[2], 1);
}

TEST(Torus, Curvature)
{
    using Real = double;
    auto mat = std::make_shared<lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> center(0,0,0);
    auto tor = torus<Real>(center, 3.0, 1.0, mat);
    auto [kappa_min, kappa_max] = tor.gaussian_curvature_bounds();
    EXPECT_FLOAT_EQ(kappa_min, -0.5);
    EXPECT_FLOAT_EQ(kappa_max, 0.25);
}

#endif