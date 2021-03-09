#ifndef DRT_AARECT_HPP
#define DRT_AARECT_HPP
#include <cmath>
#include <limits>
#include <drt/aabb.hpp>
#include <drt/ray.hpp>
#include <drt/hittable.hpp>

namespace drt {
template<typename Real>
class xy_rect : public hittable<Real> {
public:
    xy_rect() {}

    xy_rect(Real x0, Real x1, Real y0, Real y1, Real k)
        : x0_(x0), x1_(x1), y0_(y0), y1_(y1), k_(k) {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override {
        // The bounding box must have non-zero width in each dimension, so pad the Z
        // dimension a small amount.
        using std::sqrt;
        Real rt_eps = sqrt(std::numeric_limits<Real>::epsilon())/100;
        output_box = aabb(vec<Real>(x0_, y0_, k_ - rt_eps),
                          vec<Real>(x1_, y1_, k_ + rt_eps));
        return true;
    }

    virtual ~xy_rect() = default;

public:
    Real x0_, x1_, y0_, y1_, k_;
};

template<typename Real>
class xz_rect : public hittable<Real> {
public:
    xz_rect() {}

    xz_rect(Real x0, Real x1, Real z0, Real z1, Real k)
        : x0_(x0), x1_(x1), z0_(z0), z1_(z1), k_(k) {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override {
        using std::sqrt;
        Real rt_eps = sqrt(std::numeric_limits<Real>::epsilon())/100;
        output_box = aabb(vec<Real>(x0_, k_ - rt_eps, z0_),
                          vec<Real>(x1_, k_ + rt_eps, z1_));

        return true;
    }

    virtual ~xz_rect() = default;

public:
    Real x0_, x1_, z0_, z1_, k_;
};

template<typename Real>
class yz_rect : public hittable<Real> {
public:

    yz_rect(Real y0, Real y1, Real z0, Real z1, Real k)
        : y0_(y0), y1_(y1), z0_(z0), z1_(z1), k_(k) {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override {
        using std::sqrt;
        Real rt_eps = sqrt(std::numeric_limits<Real>::epsilon())/100;
        output_box = aabb(vec<Real>(k_ - rt_eps, y0_, z0_),
                          vec<Real>(k_ + rt_eps, y1_, z1_));

        return true;
    }

    virtual ~yz_rect() = default;

public:
    Real y0_, y1_, z0_, z1_, k_;
};


template<typename Real>
bool xy_rect<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    Real t = (k_ - r.origin()[2]) / r.direction()[2];
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin()[0] + t*r.direction()[0];
    auto y = r.origin()[1] + t*r.direction()[1];
    if (x < x0_ || x > x1_ || y < y0_ || y > y1_) {
        return false;
    }
    rec.u = (x-x0_)/(x1_ - x0_);
    rec.v = (y-y0_)/(y1_ - y0_);
    rec.t = t;
    auto outward_normal = vec<Real>(0, 0, 1);
    rec.set_face_normal(r, outward_normal);
    rec.p = r(t);
    return true;
}

template<typename Real>
bool xz_rect<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    auto t = (k_ - r.origin()[1]) / r.direction()[1];
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin()[0] + t*r.direction()[0];
    auto z = r.origin()[2] + t*r.direction()[2];
    if (x < x0_ || x > x1_ || z < z0_ || z > z1_)
        return false;
    rec.u = (x-x0_)/(x1_ - x0_);
    rec.v = (z-z0_)/(z1_ - z0_);
    rec.t = t;
    auto outward_normal = vec<Real>(0, 1, 0);
    rec.set_face_normal(r, outward_normal);
    rec.p = r(t);
    return true;
}

template<typename Real>
bool yz_rect<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    auto t = (k_ - r.origin()[0]) / r.direction()[0];
    if (t < t_min || t > t_max)
        return false;
    auto y = r.origin()[1] + t*r.direction()[1];
    auto z = r.origin()[2] + t*r.direction()[2];
    if (y < y0_ || y > y1_ || z < z0_ || z > z1_)
        return false;
    rec.u = (y - y0_)/(y1_ - y0_);
    rec.v = (z - z0_)/(z1_ - z0_);
    rec.t = t;
    auto outward_normal = vec<Real>(1, 0, 0);
    rec.set_face_normal(r, outward_normal);
    rec.p = r(t);
    return true;
}


}
#endif
