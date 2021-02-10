#ifndef DRT_CYLINDER_HPP
#define DRT_CYLINDER_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>
#include <drt/roots.hpp>


namespace drt {

namespace {
    static int64_t error_count = 0;
}

template<typename Real>
class cylinder : public hittable<Real> {
public:

    cylinder(Real radius, Real z_min, Real z_max, std::shared_ptr<material<Real>> mat_ptr)
       : radius_(radius), z_min_(z_min), z_max_(z_max), mat_ptr_(mat_ptr)
    {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~cylinder() = default;

    Real residual(vec<Real> const & p) const {
        return p[0]*p[0] + p[1]*p[1] - radius_*radius_;
    }

    Real expected_residual([[maybe_unused]] vec<Real> const & p) const {
        return 2*std::numeric_limits<Real>::epsilon()*radius_*radius_;
    }

private:

    // σ(u,v) = (rcos(2πu), rsin(2πu), z_min + v(z_max - z_min))
    void set_cylinder_uv(vec<Real> const & p, Real& u, Real& v) const {
        v = (p[2] - z_min_)/(z_max_ - z_min_);
        v = std::clamp(v, Real(0), Real(1));

        u = std::atan2(p[1], p[0])/(2*M_PI);
        if (u < 0) {
            u += 1;
        }
    }

    void set_fundamental_forms(hit_record<Real> & rec) const {
        rec.E = 4*M_PI*M_PI*radius_*radius_;
        rec.F = 0;
        rec.G = (z_max_ - z_min_)*(z_max_ - z_min_);

        rec.L = 4*M_PI*M_PI*radius_;
        rec.M = 0;
        rec.N = 0;
    }

    Real radius_;
    Real z_min_;
    Real z_max_;
    std::shared_ptr<material<Real>> mat_ptr_;
};

template<typename Real>
bool cylinder<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> o = r.origin();
    auto dir = r.direction();
    Real a = dir[0]*dir[0] + dir[1]*dir[1];
    Real b = 2*(o[0]*dir[0] + o[1]*dir[1]);
    Real c = o[0]*o[0] + o[1]*o[1] - radius_*radius_;

    std::vector<Real> roots = quadratic_roots(a, b, c);
    if (roots.size() == 0) {
        return false;
    }
    Real root = roots[0];
    if (root < t_min || t_max < root) {
        if (roots.size() == 1) {
            return false;
        }
        root = roots[1];
        if (root < t_min || t_max < root) {
            return false;
        }
    }
    rec.t = root;
    rec.p = r(root);
    Real z = rec.p[2];
    if (z > z_max_ || z < z_min_) {
        return false;
    }

    vec<Real> gradient(2*rec.p[0], 2*rec.p[1], 0);
    rec.gradient_magnitude = norm(gradient);
    vec<Real> outward_normal = gradient/rec.gradient_magnitude;
    rec.set_face_normal(r, outward_normal);
    set_cylinder_uv(outward_normal, rec.u, rec.v);
    set_fundamental_forms(rec);
    rec.mat_ptr = mat_ptr_;

    Real residual = this->residual(rec.p);
    Real expected_residual = this->expected_residual(rec.p);
    if (std::abs(residual) > expected_residual) {
#ifdef DEBUG
        std::cerr << "Residual is " << residual << ", but expected residual is " << expected_residual << "\n";
        std::cerr << "Ray: " << r << ", [t_min, t_max] = ["  << t_min << ", " << t_max << "]\n";
        std::cerr << rec << "\n";
        std::cerr << "Error count: " << ++error_count << "\n";
        std::cerr << "\n";
#endif
        return false;
    }

    return true;
}

template<typename Real>
bool cylinder<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(
        vec<Real>(-radius_, -radius_, z_min_),
        vec<Real>(radius_, radius_, z_max_));
    return true;
}

}
#endif