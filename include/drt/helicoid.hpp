#ifndef DRT_HELICOID_HPP
#define DRT_HELICOID_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/matrix.hpp>
#include <drt/roots.hpp>
#include <drt/newton.hpp>
#include <drt/halley.hpp>
#include <drt/tensor.hpp>
#include <drt/bounds.hpp>


namespace drt {

namespace {
    static int64_t helicoid_error_count = 0;
    static int64_t helicoid_hits = 0;
    static int64_t helicoid_misses = 0;
}

// Parametrized via σ(u,v) = (rvcos(2πu), rvsin(2πu), λ(u-1/2))
// Mathematica: ParametricPlot3D[{v*Cos[2*\[Pi]*u], v*Sin[2*\[Pi]*u], (u - 1/2)}, {u, 0, 1}, {v, 0, 1}]
template<typename Real>
class helicoid : public hittable<Real> {
public:

    helicoid(Real radius, Real speed)
       : radius_(radius), speed_(speed)
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
        #ifdef DEBUG
        if (helicoid_error_count > 0) {
            std::cerr << "Helicoid error count = " << helicoid_error_count << "\n";
        }
        #endif
    };

    vec<Real,3> operator()(Real u, Real v) const override {
        vec<Real> w;
        w[0] = radius_*v*std::cos(2*M_PI*u);
        w[1] = radius_*v*std::sin(2*M_PI*u);
        w[2] = speed_*(u-0.5);
        return w;
    }

    vec<Real,3> operator()(Real u, Real v, matrix<Real,3,3>& J, tensor<Real, 3,2,2>& H) const override {
        vec<Real> w;
        w[0] = radius_*v*std::cos(2*M_PI*u);
        w[1] = radius_*v*std::sin(2*M_PI*u);
        w[2] = speed_*(u-0.5);

        J(0,2) = std::numeric_limits<Real>::quiet_NaN();
        J(1,2) = std::numeric_limits<Real>::quiet_NaN();
        J(2,2) = std::numeric_limits<Real>::quiet_NaN();

        J(0,0) = -2*M_PI*w[1];
        J(1,0) = 2*M_PI*w[0];
        J(2,0) = speed_;

        J(0,1) = radius_*std::cos(2*M_PI*u);
        J(1,1) = radius_*std::sin(2*M_PI*u);
        J(2,1) = 0;

        H(0,0,0) = -4*M_PI*M_PI*w[0];
        H(0,0,1) = -2*M_PI*sin(2*M_PI*u);

        H(0,1,0) = H(0,0,1);
        H(0,1,1) = 0;

        H(1,0,0) = -4*M_PI*M_PI*w[1];
        H(1,0,1) = 2*M_PI*cos(2*M_PI*u);

        H(1,1,0) = H(1,0,1);
        H(1,1,1) = 0;

        for (int64_t j = 0; j < 2; ++j) {
            for (int64_t k = 0; k < 2; ++k) {
                H(2,j,k) = 0;
            }
        }

        return w;
    }

    template<int64_t p = 0>
    auto derivatives(Real u, Real v) const {
        static_assert(p >= 0 && p <= 2, "Not implemented past two derivatives.");
        vec<Real> w;
        w[0] = radius_*v*std::cos(2*M_PI*u);
        w[1] = radius_*v*std::sin(2*M_PI*u);
        w[2] = speed_*(u-0.5);
        if constexpr (p >= 1) {
            // Note that a function σ:ℝ² -> ℝ³, the Jacobian is a 3x2 matrix.
            // However, that's pretty much useless for the purpose people actually use Jacobians for.
            // Hence, this returns a 3x3 matrix. The final column is filled with nans,
            // and must be populated by the user. The most common filling would be with -d for a ray
            // with r(t) = o + td.
            matrix<Real,3,3> J;
            J(0,2) = std::numeric_limits<Real>::quiet_NaN();
            J(1,2) = std::numeric_limits<Real>::quiet_NaN();
            J(2,2) = std::numeric_limits<Real>::quiet_NaN();

            J(0,0) = -2*M_PI*w[1];
            J(1,0) = 2*M_PI*w[0];
            J(2,0) = speed_;

            J(0,1) = radius_*std::cos(2*M_PI*u);
            J(1,1) = radius_*std::sin(2*M_PI*u);
            J(2,1) = 0;

            // H(i,j,k) := ∂ⱼ∂ₖσᵢ = H(i,k,j).
            if constexpr (p == 2) {
                tensor<Real, 3, 2, 2> H;
                H(0,0,0) = -4*M_PI*M_PI*w[0];
                H(0,0,1) = -2*M_PI*sin(2*M_PI*u);

                H(0,1,0) = H(0,0,1);
                H(0,1,1) = 0;

                H(1,0,0) = -4*M_PI*M_PI*w[1];
                H(1,0,1) = 2*M_PI*cos(2*M_PI*u);

                H(1,1,0) = H(1,0,1);
                H(1,1,1) = 0;

                for (int64_t j = 0; j < 2; ++j) {
                    for (int64_t k = 0; k < 2; ++k) {
                        H(2,j,k) = 0;
                    }
                }
                return std::make_tuple(w, J, H);
            }
            else {
                return std::make_pair(w, J);
            }
        }
        else {
            return w;
        }
    }

    std::pair<Real, Real> gaussian_curvature_bounds() const {
        Real Kmin = -1/(speed_*speed_);
        Real rt_denom = speed_*speed_ + 4*M_PI*M_PI*radius_*radius_;
        Real Kmax = -speed_*speed_/(rt_denom*rt_denom);
        return std::make_pair(Kmin, Kmax);
    }

private:

    // For a ray that intersects a cylinder at {tmin, tmax}, what is the corresponding [umin, umax]?
    std::pair<Real, Real> ubounds(vec<Real> const & o, vec<Real> const & d, Real t_min, Real t_max) const
    {
        // No constraints if we're already in the cylinder:
        if (o[0]*o[0] + o[1]*o[1] < radius_*radius_ && o[2] < speed_/2 && o[2] > -speed_/2) {
            return std::make_pair<Real,Real>(0, 1);
        }

        Real umin = (o[2] + t_min*d[2])/speed_ + Real(1)/2;
        Real umax = (o[2] + t_max*d[2])/speed_ + Real(1)/2;
        if (umin > umax) {
            std::swap(umin, umax);
        }
        umin = std::max(Real(0), umin);
        umax = std::min(Real(1), umax);
        if (umin > 1 || umax < 0) {
            return std::make_pair(std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN());
        }
        return std::make_pair(umin, umax);
    }

    // For a ray that intersects a cylinder at {tmin, tmax}, what is the corresponding {vmin, vmax}?
    std::pair<Real, Real> vbounds(vec<Real> const & o, vec<Real> const & d, Real tmin, Real tmax) const
    {
        using std::sqrt;
        // If we are already in the cylinder, then we have essentially no constrains on v:
        if (o[0]*o[0] + o[1]*o[1] < radius_*radius_ && o[2] < speed_/2 && o[2] > -speed_/2) {
            return std::make_pair<Real,Real>(0, 1);
        }
        // v²r² = (oₓ+tdₓ)² + (oy+tdy)² which is minimized at
        // t_c = -(oₓdₓ + oydz)/(dₓ² + dy²).
        Real tc = -(o[0]*d[0] + o[1]*d[1])/(d[0]*d[0] + d[1]*d[1]);
        if (tc < tmin || tc > tmax)
        {
            // Then (ox+tdx)² + (oy+tdy)² in monotonic on [t_min, t_max].
            Real x = o[0] + tmin*d[0];
            Real y = o[1] + tmin*d[1];
            Real vmin = sqrt(x*x + y*y)/radius_;
            x = o[0] + tmax*d[0];
            y = o[1] + tmax*d[1];
            Real vmax = sqrt(x*x + y*y)/radius_;

            if (vmin > vmax) {
                std::swap(vmin, vmax);
            }
            if (vmin < 0) {
                vmin = 0;
            }
            if (vmax > 1) {
                vmax = 1;
            }
            if (vmin > 1 || vmax < 0) {
                return std::make_pair(std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN());
            }
        }

        // The critical point is in the interval [tmin, tmax].
        Real x = o[0] + tmin*d[0];
        Real y = o[1] + tmin*d[1];
        Real vmax1 = sqrt(x*x + y*y)/radius_;
        x = o[0] + tmax*d[0];
        y = o[1] + tmax*d[1];
        Real vmax2 = sqrt(x*x + y*y)/radius_;
        x = o[0] + tc*d[0];
        y = o[1] + tc*d[1];
        Real vmin = sqrt(x*x + y*y)/radius_;
        if (vmax1 < vmax2) {
            return std::make_pair(vmin, vmax1);
        }
        return std::make_pair(vmin, vmax2);
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

        // Yes, we do have to take care of this case separately. Sad!
        if (d[0]*d[0] + d[1]*d[1] <= std::numeric_limits<Real>::min())
        {
            v = sqrt(o[0]*o[0] + o[1]*o[1])/radius_;
            // We're shooting a ray down, but the origin is outside of r:
            if (v > 1) {
                return w;
            }
            if (v < std::numeric_limits<Real>::min()) {
                // Infinitely many solutions on
                // oz + tdz = λ(u-1/2)
                if (d[2] > 0) {
                    u = 0;
                    t = (-speed_/2 - o[2])/d[2];
                }
                else
                {
                    u = 1;
                    t = (speed_/2 - o[2])/d[2];
                }
            }
            else {
                Real rv = radius_*v;
                u = atan2(o[1]/rv, o[0]/rv)/(2*M_PI);
                if (u < 0) {
                    u += 1;
                }
                t = (speed_*(u-0.5) - o[2])/d[2];
            }
            std::get<0>(w) = t;
            std::get<1>(w) = std::clamp(u, Real(0), Real(1));
            std::get<2>(w) = std::clamp(v, Real(0), Real(1));
            return w;
        }

        auto [vmin, vmax] = this->vbounds(o, d, t_min, t_max);
        drt::bounds<Real, 2> b({t_min, t_max}, {vmin, vmax});
        auto f = [&o, &d, this](vec<Real,2> w) {
            Real t = w[0];
            Real v = w[1];
            vec<Real, 2> g;
            Real k = 2*M_PI/speed_;
            Real x = o[0] + t*d[0];
            Real y = o[1] + t*d[1];
            Real z = o[2] + t*d[2];
            Real rckz = radius_*cos(k*z);
            Real rskz = radius_*sin(k*z);
            Real rvckz = v*rckz;
            Real rvskz = v*rskz;
            g[0] = x + rvckz;
            g[1] = y + rvskz;
            matrix<Real, 2,2> J;
            J(0,0) = d[0] - k*d[2]*rvskz;
            J(1,0) = d[1] + k*d[2]*rvckz;
            J(0,1) = rckz;
            J(1,1) = rskz;
            return std::make_pair(g, J);
        };

        vec<Real, 2> sol = newton<Real,2>(f, b, vec<Real,2>(t_min, vmax*(1- 1000*std::numeric_limits<Real>::epsilon())));
        /*auto f = [&o, &d, this](vec<Real,2> w) {
            Real t = w[0];
            Real v = w[1];
            vec<Real, 2> g;
            Real k = 2*M_PI/speed_;
            Real x = o[0] + t*d[0];
            Real y = o[1] + t*d[1];
            Real z = o[2] + t*d[2];
            Real rckz = radius_*cos(k*z);
            Real rskz = radius_*sin(k*z);
            Real rvckz = v*rckz;
            Real rvskz = v*rskz;
            g[0] = x + rvckz;
            g[1] = y + rvskz;
            matrix<Real, 2,2> J;
            Real kdz = k*d[2];
            J(0,0) = d[0] - kdz*rvskz;
            J(1,0) = d[1] + kdz*rvckz;
            J(0,1) = rckz;
            J(1,1) = rskz;
            tensor<Real, 2,2,2> H;
            H(0,0,0) = -kdz*kdz*rvckz;
            H(0,0,1) = -kdz*rskz;
            H(0,1,0) = H(0,0,1);
            H(0,1,1) = 0;
            H(1,0,0) = -kdz*kdz*rvskz;
            H(1,0,1) = kdz*rckz;
            H(1,1,0) = H(1,0,1);
            H(1,1,1) = 0;
            return std::make_tuple(g, J, H);
        };
        vec<Real, 2> sol = halley<Real,2>(f, b, vec<Real,2>(t_min, vmax*(1- 1000*std::numeric_limits<Real>::epsilon())));*/
        t = sol[0];
        v = sol[1];
        u = (o[2] + t*d[2])/speed_ + 0.5;

        std::get<0>(w) = t;
        std::get<1>(w) = std::clamp(u, Real(0), Real(1));
        std::get<2>(w) = std::clamp(v, Real(0), Real(1));
        return w;
    }

    Real radius_;
    Real speed_;
};

template<typename Real>
bool helicoid<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    using std::abs;
    vec<Real, 3> o = r.origin();
    auto dir = r.direction();
    Real a = dir[0]*dir[0] + dir[1]*dir[1];
    // Rays vertically down do hit the helicoid, without hitting the walls of the cylinder.
    if (a > std::numeric_limits<Real>::min()) {
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
        if (cylinder_intersection_z_max < -speed_/2) {
            return false;
        }
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
    // Simplified:
    //vec<Real> gradient2(-radius_*speed_*sin(2*M_PI*rec.u), radius_*speed_*cos(2*M_PI*rec.u), -2*M_PI*radius_*radius_*rec.v);
    rec.gradient_magnitude = norm(gradient);
    vec<Real> outward_normal = gradient/rec.gradient_magnitude;
    // Could verify that dot(d2sigmadu2, outward_normal) == 0.
    rec.e = 0;
    rec.f = dot(d2sigmadudv, outward_normal);
    rec.g = 0;
    rec.E = dot(dsigmadu, dsigmadu);
    rec.F = 0;
    rec.G = 4*M_PI*M_PI*radius_*radius_;
    rec.set_face_normal(r, outward_normal);
    auto [p, J] = this->derivatives<1>(rec.u, rec.v);
    J(0,2) = -dir[0];
    J(1,2) = -dir[1];
    J(2,2) = -dir[2];
    Real residual = max_norm(rec.p - p);
    Real r0 = abs(rec.t*J(0,2)) + abs(J(0,0)*rec.u) + abs(J(0,1)*rec.v);
    Real r1 = abs(rec.t*J(1,2)) + abs(J(1,0)*rec.u) + abs(J(1,1)*rec.v);
    Real r2 = abs(rec.t*J(2,2)) + abs(J(2,0)*rec.u) + abs(J(2,1)*rec.v);

    auto M = inverse(J);
    rec.condition_number = M.max_norm()/std::max({rec.u, rec.v,rec.t});

    Real expected_residual = std::max({r0,r1,r2});
    expected_residual *= std::numeric_limits<Real>::epsilon();
    if (residual > 10*expected_residual) {
        ++helicoid_error_count;
#ifdef DEBUG
        std::cerr << "Residual for a helicoid with r = " << radius_ << " and λ = " << speed_ << " is " << residual;
        std::cerr << ", but expected residual is " << expected_residual << ".\n";
        std::cerr << "r(" << rec.t << ") = " << r(rec.t) << ", but σ(" << rec.u << ", " << rec.v << ") = " << (*this)(rec.u, rec.v) << "\n";
        std::cerr << "Ray: " << r << ", [t_min, t_max] = ["  << t_min << ", " << t_max << "]\n";
        std::cerr << rec << "\n";
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
