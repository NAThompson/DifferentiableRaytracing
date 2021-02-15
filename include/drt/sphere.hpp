#ifndef DRT_SPHERE_H
#define DRT_SPHERE_H

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>
#include <drt/roots.hpp>

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

    // The parametrization is:
    // σ(u,v) = (cos(2πu)sin(2πv), sin(2πu)sin(2πv), cos(2πv)), u,v \in [0,1].
    void get_sphere_uv(const vec<Real>& p, Real& u, Real& v) const {
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        Real z = p[2] - center_[2];

        v = std::acos(z)/(2*M_PI);
        if (v < 0) {
            v += 1;
        }
        u = std::atan2(y, x)/(2*M_PI);
        if (u < 0) {
            u += 1;
        }
    }

private:
    // Private since rec.u, rec.v must be set before this can be called:
    void set_fundamental_forms(hit_record<Real>& rec) const {
        // First fundamental form:
        Real dsigmadu = 2*M_PI*std::sin(2*M_PI*rec.v);
        rec.E = dsigmadu*dsigmadu;
        rec.F = 0;
        rec.G = 4*M_PI*M_PI;

        // The second fundamental form of a sphere is the same as the first,
        // see Pressley, Elementary Differential Geometry, Section 8.1.
        rec.L = rec.E;
        rec.M = rec.F;
        rec.N = rec.G;
    }

    vec<Real, 3> center_;
    Real radius_;
    std::shared_ptr<material<Real>> mat_ptr_;
};

template<typename Real>
bool sphere<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> oc = r.origin() - center_;
    Real a = drt::squared_norm(r.direction());
    Real b = 2*drt::dot(oc, r.direction());
    Real c = drt::squared_norm(oc) - radius_*radius_;
    auto roots = quadratic_roots(a, b, c);
    if (roots.size() != 2) {
        return false;
    }
    Real root = std::numeric_limits<Real>::quiet_NaN();
    if (roots[0] < t_min || roots[0] > t_max) {
        if (roots[1] >= t_min && roots[1] <= t_max) {
            root = roots[1];
        }
    }
    else {
        root = roots[0];
    }

    if (std::isnan(root)) {
        return false;
    }
    rec.t = root;
    rec.p = r(root);
    auto outward_normal = (rec.p - center_) / radius_;
    rec.set_face_normal(r, outward_normal);
    get_sphere_uv(outward_normal, rec.u, rec.v);
    set_fundamental_forms(rec);
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