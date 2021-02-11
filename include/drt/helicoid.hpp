#ifndef DRT_HELICOID_HPP
#define DRT_HELICOID_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>
#include <drt/roots.hpp>


namespace drt {

namespace {
    static int64_t helicoid_error_count = 0;
}

template<typename Real>
class helicoid : public hittable<Real> {
public:

    helicoid(Real radius, Real speed, std::shared_ptr<material<Real>> mat_ptr)
       : radius_(radius), speed_(speed), mat_ptr_(mat_ptr)
    {
        if (radius_ <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Radius > 0 is required for a helicoid.\n";
        }
        if (speed_ <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Speed > 0 is required for a helicoid.\n";
        }
    };

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~helicoid() = default;

    Real residual([[maybe_unused]] vec<Real> const & p) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    Real expected_residual([[maybe_unused]] vec<Real> const & p) const {
        return std::numeric_limits<Real>::quiet_NaN();
    }

private:

    // σ(u,v) = (rvcos(2πu), rvsin(2πu), λ(u-1/2))
    void set_helicoid_uv(vec<Real> const & p, Real& u, Real& v) const {
        v = std::sqrt(p[0]*p[0] + p[1]*p[1])/radius_;
        v = std::clamp(v, Real(0), Real(1));

        u = p[2]/speed_ + Real(1)/2;
    }

    void set_fundamental_forms(hit_record<Real> & rec) const {
        rec.E = 4*M_PI*M_PI*radius_*radius_;
        rec.F = 0;
        rec.G = 0;

        rec.L = 4*M_PI*M_PI*radius_;
        rec.M = 0;
        rec.N = 0;
    }

    std::tuple<Real, Real, Real> f(Real t, ray<Real> const & r) const {
        // Mathematica:
        // f[t_] := (Subscript[\[Sigma], x] + t*Subscript[d, x])*Tan[2*Pi*(Subscript[\[Sigma], z] + 
        // t*Subscript[d, z])/\[Lambda]] - (Subscript[\[Sigma], y] + t*Subscript[d, y])
        using std::tan;
        using std::cos;
        vec<Real> o = r.origin();
        vec<Real> d = r.direction();
        Real k = 2*M_PI/speed_;
        Real arg = k*(o[2] + t*d[2]);
        Real tn = tan(arg);
        Real sc = 1/cos(arg);
        Real x = o[0] + t*d[0];
        Real y = x*tn - (o[1] + t*d[1]);
        Real dydt = -d[1] + d[0]*tn + k*sc*sc*d[2]*x;
        Real d2ydt2 = 2*k*sc*sc*d[0]*d[2]+ 2*k*k*sc*sc*d[2]*d[2]*x*tn;
        return std::make_tuple(y, dydt, d2ydt2);
    }

    Real radius_;
    Real speed_;
    std::shared_ptr<material<Real>> mat_ptr_;
};

template<typename Real>
bool helicoid<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    using std::abs;
    vec<Real, 3> o = r.origin();
    auto dir = r.direction();
    Real a = dir[0]*dir[0] + dir[1]*dir[1];
    Real b = 2*(o[0]*dir[0] + o[1]*dir[1]);
    Real c = o[0]*o[0] + o[1]*o[1] - radius_*radius_;

    std::vector<Real> roots = quadratic_roots(a, b, c);
    if (roots.size() != 2) {
        return false;
    }
    Real helicoid_tmin = std::max(t_min, roots[0]);
    Real helicoid_tmax = std::min(t_max, roots[1]);

    Real helicoid_z_max = speed_/2;
    Real helicoid_z_min = -speed_/2;
    auto cylinder_p1 = r(helicoid_tmax);
    auto cylinder_p2 = r(helicoid_tmin);
    Real cylinder_intersection_z_max = std::max(cylinder_p1[2], cylinder_p2[2]);
    Real cylinder_intersection_z_min = std::min(cylinder_p1[2], cylinder_p2[2]);
    if (cylinder_intersection_z_min > helicoid_z_max) {
        return false;
    }
    if (cylinder_intersection_z_max < helicoid_z_min) {
        return false;
    }

    Real t = helicoid_tmin;
    auto [y, dydt, d2ydt2] = f(t, r);
    int i = 0;
    Real fudge_scale = 100;
    while (abs(y) > 10*std::sqrt(std::numeric_limits<Real>::epsilon())) {
        //t -= y/dydt;
        t -= 2*y*dydt/(2*dydt*dydt - y*d2ydt2);
        if (t > helicoid_tmax) {
            std::cerr << "Newton's method went too far forward!\n";
            std::cerr << "updated t = " << t << "\n";
            std::cerr << "error count " << ++helicoid_error_count << "\n";
            return false;
        }
        if (t < helicoid_tmin) {
            std::cerr << "Newton's method went too far backward!\n";
            std::cerr << "error count " << ++helicoid_error_count << "\n";
            return false;
        }
        std::tie(y, dydt, d2ydt2) = f(t, r);
        if (i++ > 500) {
            std::cerr << "Bad news: Hit max iterations!\n";
            std::cerr << "Residual for Newton's method: " << abs(y) << "\n";
            std::cerr << "Acceptable residual = " << fudge_scale*std::numeric_limits<Real>::epsilon()*abs(t*dydt) << "\n";
            std::cerr << "dydt = " << dydt << "\n";
            std::cerr << "t = " << t << ", should be in [tmin, tmax] = [" << helicoid_tmin << ", " << helicoid_tmax << "]\n";
            std::cerr << "error count " << ++helicoid_error_count << "\n";
            break;
        }
    }
    rec.p = r(t);
    rec.t = t;
    set_helicoid_uv(rec.p, rec.u, rec.v);

    vec<Real> dsigmadu(-2*M_PI*radius_*rec.v*sin(2*M_PI*rec.u), 2*M_PI*radius_*rec.v*cos(2*M_PI*rec.u), speed_);
    vec<Real> dsigmadv(radius_*cos(2*M_PI*rec.u), radius_*sin(2*M_PI*rec.u), 0);
 
    vec<Real> gradient = drt::cross(dsigmadu, dsigmadv);
    rec.gradient_magnitude = norm(gradient);
    vec<Real> outward_normal = gradient/rec.gradient_magnitude;
    rec.set_face_normal(r, outward_normal);
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
bool helicoid<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(
        vec<Real>(-radius_, -radius_, speed_/2),
        vec<Real>(radius_, radius_, speed_/2));
    return true;
}

}
#endif