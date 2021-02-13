#ifndef DRT_HELICOID_HPP
#define DRT_HELICOID_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>
#include <drt/roots.hpp>


namespace drt {

namespace {
    static int64_t helicoid_error_count = 0;
    static int64_t helicoid_hits = 0;
    static int64_t helicoid_misses = 0;
}

// Parametrized via σ(u,v) = (rvcos(2πu), rvsin(2πu), λ(u-1/2))
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

    virtual ~helicoid() {
        std::cerr << "Hits = " << helicoid_hits << ", misses = " << helicoid_misses << "\n";
    };

    vec<Real> sigma(Real u, Real v) const {
        Real x = radius_*v*std::cos(2*M_PI*u);
        Real y = radius_*v*std::sin(2*M_PI*u);
        Real z = speed_*(u-0.5);
        return vec<Real>(x,y,z);
    }

private:

    std::tuple<Real, Real, Real> f(Real t, vec<Real> const & o, vec<Real> const & d) const {
        // Mathematica:
        // f[t_] := (Subscript[\[Sigma], x] + t*Subscript[d, x])*Tan[2*Pi*(Subscript[\[Sigma], z] + 
        // t*Subscript[d, z])/\[Lambda]] - (Subscript[\[Sigma], y] + t*Subscript[d, y])
        using std::tan;
        using std::cos;
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

    std::tuple<Real, Real, Real> compute_tuv(vec<Real> const & o, vec<Real> const & d, Real t_min, Real t_max) const
    {
        using std::abs;
        using std::tan;
        using std::sqrt;
        std::tuple<Real, Real, Real> w(std::numeric_limits<Real>::quiet_NaN(),
                                       std::numeric_limits<Real>::quiet_NaN(),
                                       std::numeric_limits<Real>::quiet_NaN());
        Real u = std::numeric_limits<Real>::quiet_NaN();
        Real t = std::numeric_limits<Real>::quiet_NaN();
        Real v = std::numeric_limits<Real>::quiet_NaN();
        if (abs(d[2]) <= std::numeric_limits<Real>::min()) {
            u = o[2]/speed_ + Real(1)/Real(2);
            if (u < 0 || u > 1) {
                return w;
            }
            Real theta = 2*M_PI*o[2]/speed_;
            Real tant = tan(theta);
            if (abs(tant) <= std::numeric_limits<Real>::min()) {
                if (abs(d[1]) <= std::numeric_limits<Real>::min()) {
                    if (abs(o[1]) > std::numeric_limits<Real>::min()) {
                        // no solution
                        return w;
                    }
                    else {
                        // There are infinitely many solutions along the line o_x + td_x = -rv
                        // Extremal values of t occur when v = 0 or v = 1,
                        // unless o_x in one the line!
                        // Ignore this case for now.
                        Real t0 = -o[0]/d[0];
                        Real t1 = -(o[0] + 1)/d[0];
                        if (t0 < t1 && t0 >= t_min) {
                            t = t0;
                            v = 0;
                        }
                        else {
                            t = t1;
                            v = 1;
                        }
                    }
                }
                else {
                    t = -o[1]/d[1];
                    v = (o[0] + t*d[0])/(radius_*cos(2*M_PI*u));
                }
            }
            else {
                t = (o[0]*tant - o[1])/(d[1] - d[0]*tant);
                Real x = o[0] + t*d[0];
                Real y = o[1] + t*d[1];
                v = sqrt(x*x + y*y)/radius_;
            }

            if (t >= t_min && t <= t_max) {
                std::get<0>(w) = t;
                std::get<1>(w) = std::clamp(u, Real(0), Real(1));
                std::get<2>(w) = std::clamp(v, Real(0), Real(1));
            }
            return w;
        }
        // It's easier to recover u first, then t.
        Real umin = (o[2] + t_min*d[2])/speed_ + Real(1)/2;
        Real umax = (o[2] + t_max*d[2])/speed_ + Real(1)/2;
        if (umin > umax) {
            std::swap(umin, umax);
        }
        if (umin > 1) {
            return w;
        }
        umin = std::max(Real(0), umin);
        if (umax < 0) {
            return w;
        }
        umax = std::min(Real(1), umax);
        std::function<Real(Real)> g2 = [&o, &d, this](Real u) -> Real {
            Real t = (speed_*(u - 0.5) - o[2])/d[2];
            Real angle = std::atan2(o[1] + t*d[1], o[0] + t*d[0])/(2*M_PI);
            if (angle < 0) {
                angle += 1;
            }
            return angle - u;
        };
        std::tie(umin, umax) = drt::bisect(g2, umin, umax);
        if (std::isnan(umin)) {
            return w;
        }
        u = umin;
        t = (speed_*(u-Real(1)/Real(2)) - o[2])/d[2];

        Real x_r = (o[0] + t*d[0]);
        Real y_r = (o[1] + t*d[1]);
        Real z_r = (o[2] + t*d[2]);
        v = std::sqrt(x_r*x_r + y_r*y_r)/radius_;
        bool xbad = abs(x_r - radius_*v*cos(2*M_PI*u)) > 0.1;
        bool ybad = abs(y_r - radius_*v*sin(2*M_PI*u)) > 0.1;
        bool zbad = abs(z_r - speed_*(u-0.5)) > 0.1;
        if (xbad || ybad || zbad) {
            std::function<Real(Real)> g1 = [&o, &d, this](Real t) -> Real {
                Real angle = std::atan2(o[1] + t*d[1], o[0] + t*d[0]);
                if (angle < 0) {
                    angle += 2*M_PI;
                }
                return angle - 2*M_PI*(o[2] + t*d[2])/speed_ - M_PI;
            };

            std::tie(t_min, t_max) = drt::bisect(g1, t_min, t_max);
            if (std::isnan(t_min)) {
                return w;
            }
            t = t_min;
            x_r = (o[0] + t*d[0]);
            y_r = (o[1] + t*d[1]);
            z_r = (o[2] + t*d[2]);
            u = z_r/speed_ + 0.5;
            v = std::sqrt(x_r*x_r + y_r*y_r)/radius_;
            xbad = abs(x_r - radius_*v*cos(2*M_PI*u)) > 0.1;
            ybad = abs(y_r - radius_*v*sin(2*M_PI*u)) > 0.1;
            zbad = abs(z_r - speed_*(u-0.5)) > 0.1;
            if (xbad || ybad || zbad) {
                return w;
            }
        }
        std::get<0>(w) = t;
        std::get<1>(w) = std::clamp(u, Real(0), Real(1));
        std::get<2>(w) = std::clamp(v, Real(0), Real(1));
        return w;
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
    if (roots[1] < t_min) {
        return false;
    }
    t_min = std::max(t_min, roots[0]);
    t_max = std::min(t_max, roots[1]);

    auto cylinder_p1 = r(t_max);
    auto cylinder_p2 = r(t_min);
    Real cylinder_intersection_z_max = std::max(cylinder_p1[2], cylinder_p2[2]);
    Real cylinder_intersection_z_min = std::min(cylinder_p1[2], cylinder_p2[2]);
    if (cylinder_intersection_z_min > speed_/2) {
        return false;
    }
    if (cylinder_intersection_z_max <  -speed_/2) {
        return false;
    }

    std::tie(rec.t, rec.u, rec.v) = compute_tuv(r.origin(), r.direction(), t_min, t_max);
    if (std::isnan(rec.t)) {
        ++helicoid_misses;
        return false;
    }
    rec.p = r(rec.t);

    vec<Real> dsigmadu(-2*M_PI*radius_*rec.v*sin(2*M_PI*rec.u), 2*M_PI*radius_*rec.v*cos(2*M_PI*rec.u), speed_);
    vec<Real> dsigmadv(radius_*cos(2*M_PI*rec.u), radius_*sin(2*M_PI*rec.u), 0);
    vec<Real> d2sigmadu2(-4*M_PI*M_PI*radius_*rec.v*cos(2*M_PI*rec.u), -4*M_PI*M_PI*radius_*rec.v*sin(2*M_PI*rec.u), 0);
    vec<Real> d2sigmadudv(-2*M_PI*radius_*sin(2*M_PI*rec.u), 2*M_PI*radius_*cos(2*M_PI*rec.u), 0);
 
    vec<Real> gradient = drt::cross(dsigmadu, dsigmadv);
    rec.gradient_magnitude = norm(gradient);
    vec<Real> outward_normal = gradient/rec.gradient_magnitude;
    // Could verify that dot(d2sigmadu2, outward_normal) == 0.
    rec.L = 0;
    rec.M = dot(d2sigmadudv, outward_normal);
    rec.N = 0;
    rec.E = dot(dsigmadu, dsigmadu);
    rec.F = 0;
    rec.G = radius_*radius_;
    rec.mat_ptr = mat_ptr_;
    rec.set_face_normal(r, outward_normal);
    Real residual = norm(rec.p - this->sigma(rec.u, rec.v));
    //Real expected_residual = this->expected_residual(rec.p);
    if (std::abs(residual) > 1) {
        ++error_count;
#ifdef DEBUG
        std::cerr << "Residual for a helicoid with r = " << radius_ << " and λ = " << speed_ << " is " << residual << ".\n";
        std::cerr << "r(" << rec.t << ") = " << r(rec.t) << ", but σ(" << rec.u << ", " << rec.v << ") = " << this->sigma(rec.u, rec.v) << "\n";
        std::cerr << "Ray: " << r << ", [t_min, t_max] = ["  << t_min << ", " << t_max << "]\n";
        std::cerr << rec << "\n";
        std::cerr << "Error count: " << error_count << "\n";
        std::cerr << "Hits = " << helicoid_hits << ", misses = " << helicoid_misses << "\n";
        std::cerr << "\n";
#endif
        ++helicoid_misses;
        return false;
    }

    ++helicoid_hits;
    return true;
}

template<typename Real>
bool helicoid<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(
        vec<Real>(-radius_, -radius_, -speed_/2),
        vec<Real>(radius_, radius_, speed_/2));
    return true;
}

}
#endif