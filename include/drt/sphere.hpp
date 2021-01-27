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

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~sphere() = default;
public:

    static void get_sphere_uv(const vec<Real>& p, Real& u, Real& v) {
        auto theta = std::acos(-p[1]);
        auto phi = std::atan2(-p[2], p[0]) + M_PI;

        u = phi / (2*M_PI);
        v = theta / M_PI;
    }

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
    get_sphere_uv(outward_normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr_;
    return true;
}

template<typename Real>
bool sphere<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(
        center_ - vec<Real>(radius_, radius_, radius_),
        center_ + vec<Real>(radius_, radius_, radius_));
    return true;
}

}
#endif