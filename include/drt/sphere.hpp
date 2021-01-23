#ifndef DRT_SPHERE_H
#define DRT_SPHERE_H

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>

namespace drt {
template<typename Real>
class sphere : public hittable<Real> {
public:
    sphere() {}
    sphere(vec<Real, 3> const & center, Real radius, std::shared_ptr<material<Real>> mat_ptr)
       : center_(center), radius_(radius), mat_ptr_(mat_ptr)
    {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual ~sphere() = default;
public:
    vec<Real, 3> center_;
    Real radius_;
    std::shared_ptr<material<Real>> mat_ptr_;
};

template<typename Real>
bool sphere<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> oc = r.origin() - center_;
    Real a = drt::squared_norm(r.direction());
    Real half_b = drt::dot(oc, r.direction());
    auto c = drt::squared_norm(oc) - radius_*radius_;

    auto discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return false;
    Real sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r(root);
    auto outward_normal = (rec.p - center_) / radius_;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr_;
    return true;
}

}
#endif