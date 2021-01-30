#ifndef DRT_ELLIPSOID_HPP
#define DRT_ELLIPSOID_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>

namespace drt {
template<typename Real>
class ellipsoid : public hittable<Real> {
public:
    ellipsoid() {}
    ellipsoid(vec<Real, 3> const & center, Real a, Real b, Real c, std::shared_ptr<material<Real>> mat_ptr)
       : center_(center), a_(a), b_(b), c_(c), mat_ptr_(mat_ptr)
    {
        if (a_ <= 0) {
            std::cerr << "a <= 0 is not allowed for an ellipsoid.\n";
        }
        if (b_ <= 0) {
            std::cerr << "b <= 0 is not allowed for an ellipsoid.\n";
        }
        if (c_ <= 0) {
            std::cerr << "c <= 0 is not allowed for an ellipsoid.\n";
        }
    };

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;


    Real gaussian_curvature(const vec<Real>& p) {
        Real asq = a_*a_;
        Real bsq = b_*b_;
        Real csq = c_*c_;
        // https://mathworld.wolfram.com/Ellipsoid.html, equation 14:s
        Real numerator = asq*bsq*bsq*bsq*csq*csq*csq;
        Real sqrt_denom = csq*csq*bsq*bsq + csq*csq*(asq - bsq)*(p[1]-center_[1])*(p[1] -center_[1]) + bsq*bsq*(asq-csq)*(p[2] -center_[2])*(p[2] -center_[2]);
        return numerator/(sqrt_denom*sqrt_denom);
    }

    virtual ~ellipsoid() = default;
public:

    static void get_sphere_uv(const vec<Real>& p, Real& u, Real& v) {
        auto theta = std::acos(-p[1]);
        auto phi = std::atan2(-p[2], p[0]) + M_PI;

        u = phi / (2*M_PI);
        v = theta / M_PI;
    }

    vec<Real, 3> center_;
    Real a_, b_, c_;
    std::shared_ptr<material<Real>> mat_ptr_;
};

template<typename Real>
bool ellipsoid<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> oc = r.origin() - center_;
    oc[0] /= a_;
    oc[1] /= b_;
    oc[2] /= c_;
    vec<Real> sd = r.direction();
    sd[0] /= a_;
    sd[1] /= b_;
    sd[2] /= c_;

    Real a = drt::squared_norm(sd);
    Real half_b = drt::dot(oc, sd);
    auto c = drt::squared_norm(oc) - 1;

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
    vec<Real> outward_normal = (rec.p - center_);
    outward_normal[0] *= 2/(a_*a_);
    outward_normal[1] *= 2/(b_*b_);
    outward_normal[2] *= 2/(c_*c_);
    normalize(outward_normal);

    rec.set_face_normal(r, outward_normal);
    get_sphere_uv(outward_normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr_;
    return true;
}

template<typename Real>
bool ellipsoid<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(
        center_ - vec<Real>(a_, b_, c_),
        center_ + vec<Real>(a_, b_, c_));
    return true;
}

}
#endif