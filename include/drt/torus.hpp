#ifndef DRT_TORUS_HPP
#define DRT_TORUS_HPP
#include <utility>
#include <algorithm>
#include <cmath>
#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>
#include <drt/roots.hpp>

static int64_t error_count = 0;
namespace drt {
template<typename Real>
class torus : public hittable<Real> {
public:
    torus(vec<Real, 3> const & center, Real major_radius, Real minor_radius, std::shared_ptr<material<Real>> mat_ptr)
       : center_(center), R_(major_radius), r_(minor_radius), mat_ptr_(mat_ptr)
    {
        if (R_ < 0) {
            std::cerr << "Major radius must be >= 0\n";
        }
        if (r_ < 0) {
            std::cerr << "Minor radius must be >= 0\n";
        }
    };

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    Real area() const {
        return 4*M_PI*M_PI*R_*r_;
    }

    Real volume() const {
        return 2*M_PI*M_PI*R_*r_*r_;
    }

    vec<Real,3> normal(vec<Real, 3> const & p) const {
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        Real z = p[2] - center_[2];
        Real s = x*x + y*y + z*z - r_*r_;
        vec<Real> n(4*x*(s- R_*R_), 4*y*(s - R_*R_), 4*z*(s + R_*R_));
        normalize(n);
        return n;
    }

    Real gaussian_curvature(vec<Real, 3> const & p) const {
        using std::sqrt;
        Real zr = (p[2] - center_[2])/r_;
        std::clamp(zr, -Real(1), Real(1));
        // cos(v) = cos(arcsin(z/r)) = ±√(1-z²/r²).
        Real cosv = sqrt(1-zr*zr);
        // We need additional information to resolve the sign of the square root.
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        // x² + y² = (R+rcos(v))² = R² + r²cos(v)² + 2rRcos(v)
        if (x*x + y*y - (R_*R_ + r_*r_*cosv*cosv) < 0) {
            cosv = -cosv;
        }
        return cosv/(r_*(R_ + r_*cosv));
    }

    std::pair<Real, Real> gaussian_curvature_bounds() const {
        return std::make_pair<Real,Real>(-1/(r_*(R_ - r_)), 1/(r_*(R_ + r_)));
    }

    // For a point p approximately on the torus,
    // what is f(p)? In exact arithmetic, f(p) = 0 identically,
    // but due to floating point error, f(p) = residual.
    Real residual(vec<Real, 3> const & p) const {
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        Real z = p[2] - center_[2];
        Real tmp = x*x + y*y + z*z + R_*R_ - r_*r_;
        Real residual = tmp*tmp - 4*R_*R_*(x*x + y*y);
        return residual;
    }

    // Under the floating point rounding model x = x(1+u),
    // what is the residual we should expect?
    Real expected_residual(vec<Real, 3> const & p) const {
        using std::abs;
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        Real z = p[2] - center_[2];
        Real length_sq = x*x + y*y + z*z;

        Real dfdx = 4*x*(length_sq - R_*R_ - r_*r_);
        Real dfdy = 4*y*(length_sq - R_*R_ - r_*r_);
        Real dfdz = 4*z*(length_sq + R_*R_ - r_*r_);

        Real residual = abs(x*dfdx) + abs(y*dfdy) + abs(z*dfdz);
        // Factor of 100 just to save some annoyance:
        residual *= 100*std::sqrt(std::numeric_limits<Real>::epsilon());
        return residual;
    }

    // x = (R + r cos(v)) cos(u)
    // y = (R + r cos(v)) sin(u)
    // z = r sin(v)
    std::pair<Real, Real> vec_to_uv(vec<Real> const & p)
    {
        using std::asin;
        using std::atan2;
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        Real z = p[2] - center_[2];
        // Bunch of clamps to stop nans:
        std::clamp(x, -(R_+r_), (R_+r_));
        std::clamp(y, -(R_+r_), (R_+r_));
        std::clamp(z, -r_, r_);
        Real u = atan2(x,y);
        if(u < 0) {
            u += 2*M_PI;
        }
        Real v = asin(z/r_);
        if (v < 0) {
            v += 2*M_PI;
        }
        return std::pair<Real, Real>(u,v);
    }

    virtual ~torus() = default;

private:
    vec<Real, 3> center_;
    Real R_; // major radius
    Real r_; // minor radius
    std::shared_ptr<material<Real>> mat_ptr_;
};

template<typename Real>
bool torus<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> oc = r.origin() - center_;
    vec<Real, 5> c;
    // http://blog.marcinchwedczuk.pl/ray-tracing-torus is informative,
    // but gives a torus radially symmetric about the y-axis.
    // Wikipedia gives a torus radially symmetric about the z-axis.
    // The Mathematica command is:
    // FullSimplify[Expand[((ocx + t*dx)^2 + (ocy + t*dy)^2 + (ocz + t*dz)^2 + R^2 - r^2)^2 - 4*R^2*((ocx + t*dx)^2 + (ocy + t*dy)^2)]]
    Real nsq = oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2];
    Real tmp = (nsq - (r_*r_ + R_*R_));
    auto d = r.direction();
    Real dsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    Real ocd = dot(oc, d);
    c[0] = tmp*tmp - 4*R_*R_*(r_*r_ - oc[2]*oc[2]);
    c[1] = 4*tmp*ocd + 8*R_*R_*oc[2]*d[2];
    c[2] = 2*dsq*tmp + 4*ocd*ocd + 4*R_*R_*d[2]*d[2];
    c[3] = 4*dsq*ocd;
    c[4] = dsq*dsq;

    auto roots = quartic_roots(c[4], c[3], c[2], c[1], c[0]);
    if (roots.size() == 0) {
        return false;
    }

    Real t = std::numeric_limits<Real>::quiet_NaN();
    // Roots are sorted so this gets the minimal root:
    for (auto root : roots) {
        if (root >= t_min && root <= t_max) {
            t = root;
            break;
        }
    }
    if (std::isnan(t)) {
        return false;
    }

    rec.t = t;
    rec.p = r(rec.t);
    vec<Real> outward_normal = this->normal(rec.p);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr_;
    using std::abs;
    Real res = this->residual(rec.p);
    Real expected_res = this->expected_residual(rec.p);
    if (abs(res) > expected_res) {
#ifdef DEBUG
        error_count++;
        std::cerr << __FILE__ << ":" << __LINE__ << " Residual for torus intersection unexpectedly high. ";
        std::cerr << "Residual is " << res << ", but expected residual is " << expected_res << ".\n";
        std::cerr << rec << "\n";
        std::cerr << "[t_min, t_max] = [" << t_min << ", " << t_max << "]\n";
        std::cerr << "Ray: " << r << "\n";
        std::cerr << "Error count = " << error_count << "\n";
        std::cerr << "Roots are {";
        for (auto r : roots) {
            std::cerr << r << ", ";
        }
        std::cerr << "}\n";
#endif
        return false;
    }


    return true;
}

template<typename Real>
bool torus<Real>::bounding_box(aabb<Real>& output_box) const {
    vec<Real> shift(R_ + r_, R_ + r_, r_);
    output_box = aabb<Real>(center_ - shift, center_ + shift);
    return true;
}

}
#endif