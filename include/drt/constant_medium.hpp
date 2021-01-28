#ifndef DRT_CONSTANT_MEDIUM_HPP
#define DRT_CONSTANT_MEDIUM_HPP
#include <memory>
#include <cmath>
#include <drt/hittable.hpp>
#include <drt/material.hpp>
#include <drt/texture.hpp>
#include <drt/isotropic.hpp>

namespace drt {

template<typename Real>
class constant_medium : public hittable<Real> {
public:
    constant_medium(std::shared_ptr<hittable<Real>> b, Real d, std::shared_ptr<texture<Real>> a)
        : boundary_(b), neg_inv_density_(-1/d), phase_function_(std::make_shared<isotropic<Real>>(a))
    {
        dis_ = std::uniform_real_distribution<Real>(0,1);
    }

    constant_medium(std::shared_ptr<hittable<Real>> b, Real d, vec<Real> c)
        : boundary_(b), neg_inv_density_(-1/d), phase_function_(std::make_shared<isotropic<Real>>(c))
    {
        dis_ = std::uniform_real_distribution<Real>(0,1);
    }

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override {
        return boundary_->bounding_box(output_box);
    }

    virtual ~constant_medium() = default;

public:
    std::shared_ptr<hittable<Real>> boundary_;
    Real neg_inv_density_;
    std::shared_ptr<material<Real>> phase_function_;
    mutable std::uniform_real_distribution<Real> dis_;
    mutable std::mt19937_64 gen_;
};

template<typename Real>
bool constant_medium<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const
{
    hit_record<Real> rec1, rec2;

    if (!boundary_->hit(r, -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity(), rec1))
        return false;

    if (!boundary_->hit(r, rec1.t + 100*std::numeric_limits<Real>::epsilon(), std::numeric_limits<Real>::infinity(), rec2))
        return false;


    if (rec1.t < t_min) rec1.t = t_min;
    if (rec2.t > t_max) rec2.t = t_max;

    if (rec1.t >= rec2.t)
        return false;

    if (rec1.t < 0)
        rec1.t = 0;

    const Real ray_length = drt::norm(r.direction());
    const Real distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
    const Real hit_distance = neg_inv_density_ * std::log(dis_(gen_));

    if (hit_distance > distance_inside_boundary)
        return false;

    rec.t = rec1.t + hit_distance / ray_length;
    rec.p = r(rec.t);


    rec.normal = vec<Real>(1,0,0);  // arbitrary
    rec.front_face = true;     // also arbitrary
    rec.mat_ptr = phase_function_;

    return true;
}

}
#endif