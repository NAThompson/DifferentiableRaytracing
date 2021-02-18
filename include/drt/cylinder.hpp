#ifndef DRT_CYLINDER_HPP
#define DRT_CYLINDER_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>
#include <drt/roots.hpp>


namespace drt {

namespace {
    static int64_t cylinder_error_count = 0;
}

template<typename Real>
class cylinder : public hittable<Real> {
public:

    cylinder(Real radius, Real z_min, Real z_max)
       : radius_(radius), z_min_(z_min), z_max_(z_max)
    {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~cylinder() = default;

    Real residual(vec<Real> const & p) const {
        return p[0]*p[0] + p[1]*p[1] - radius_*radius_;
    }

    Real expected_residual([[maybe_unused]] vec<Real> const & p) const {
        return 20*std::numeric_limits<Real>::epsilon()*radius_*radius_;
    }

    // σ(u,v) = (rcos(2πu), rsin(2πu), z_min + v(z_max - z_min))
    std::pair<Real, Real> get_uv(vec<Real> const & p) const {
        Real v = (p[2] - z_min_)/(z_max_ - z_min_);
        v = std::clamp(v, Real(0), Real(1));

        Real u = std::atan2(p[1], p[0])/(2*M_PI);
        if (u < 0) {
            u += 1;
        }
        return std::make_pair(u,v);
    }

private:

    void set_fundamental_forms(hit_record<Real> & rec) const {
        rec.E = 4*M_PI*M_PI*radius_*radius_;
        rec.F = 0;
        rec.G = (z_max_ - z_min_)*(z_max_ - z_min_);

        rec.e = 4*M_PI*M_PI*radius_;
        rec.f = 0;
        rec.g = 0;
    }

    Real radius_;
    Real z_min_;
    Real z_max_;
};

template<typename Real>
bool cylinder<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> o = r.origin();
    auto dir = r.direction();
    Real a = dir[0]*dir[0] + dir[1]*dir[1];
    Real b = 2*(o[0]*dir[0] + o[1]*dir[1]);
    Real c = o[0]*o[0] + o[1]*o[1] - radius_*radius_;
    // This might be wrong if the first root is in the t-range, but gives a ray
    // which is not in the z-range.
    auto opt_root = first_quadratic_root_in_range(a, b, c, t_min, t_max);
    if (!opt_root) {
        return false;
    }
    rec.t = *opt_root;
    rec.p = r(rec.t);
    Real z = rec.p[2];
    if (z > z_max_ || z < z_min_) {
        return false;
    }

    vec<Real> gradient(2*rec.p[0], 2*rec.p[1], 0);
    rec.gradient_magnitude = norm(gradient);
    vec<Real> outward_normal = gradient/rec.gradient_magnitude;
    rec.set_face_normal(r, outward_normal);
    std::tie(rec.u, rec.v) = get_uv(rec.p);
    set_fundamental_forms(rec);

    Real residual = this->residual(rec.p);
    Real expected_residual = this->expected_residual(rec.p);
    if (std::abs(residual) > expected_residual) {
#ifdef DEBUG
        std::cerr << "Residual is " << residual << ", but expected residual is " << expected_residual << "\n";
        std::cerr << "Ray: " << r << ", [t_min, t_max] = ["  << t_min << ", " << t_max << "]\n";
        std::cerr << rec << "\n";
        std::cerr << "Error count: " << ++cylinder_error_count << "\n";
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