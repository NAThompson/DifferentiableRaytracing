#ifndef DRT_AARECT_HPP
#define DRT_AARECT_HPP
#include <cmath>
#include <limits>
#include <drt/aabb.hpp>
#include <drt/ray.hpp>
#include <drt/hittable.hpp>
#include <drt/material.hpp>

namespace drt {
template<typename Real>
class xy_rect : public hittable<Real> {
public:
    xy_rect() {}

    xy_rect(Real x0, Real x1, Real y0, Real y1, Real k,
           std::shared_ptr<material<Real>> mat)
        : x0_(x0), x1_(x1), y0_(y0), y1_(y1), k_(k), mp_(mat) {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override {
        // The bounding box must have non-zero width in each dimension, so pad the Z
        // dimension a small amount.
        using std::sqrt;
        Real rt_eps = sqrt(std::numeric_limits<Real>::epsilon())/100;
        output_box = aabb(vec<Real>(x0_, y0_, k_*(1-rt_eps)),
                          vec<Real>(x1_, y1_, k_*(1+rt_eps)));
        return true;
    }

    virtual ~xy_rect() = default;

public:
    Real x0_, x1_, y0_, y1_, k_;
    std::shared_ptr<material<Real>> mp_;
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
    rec.mat_ptr = mp_;
    rec.p = r(t);
    return true;
}

}
#endif