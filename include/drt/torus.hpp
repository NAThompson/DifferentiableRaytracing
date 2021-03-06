#ifndef DRT_TORUS_HPP
#define DRT_TORUS_HPP
#include <utility>
#include <algorithm>
#include <cmath>
#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/roots.hpp>


namespace drt {

template<typename Real>
class torus : public hittable<Real> {
public:
    torus(vec<Real, 3> const & center, Real major_radius, Real minor_radius)
       : center_(center), R_(major_radius), r_(minor_radius)
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
        vec<Real> n(4*x*(s - R_*R_), 4*y*(s - R_*R_), 4*z*(s + R_*R_));
        normalize(n);
        return n;
    }

    Real gaussian_curvature(vec<Real, 3> const & p) const {
        using std::sqrt;
        Real zr = (p[2] - center_[2])/r_;
        Real arg = 1 - zr*zr;
        arg = std::clamp(arg, Real(0), Real(1));
        // cos(v) = cos(arcsin(z/r)) = ±√(1-z²/r²).
        Real cosv = sqrt(arg);
        // We need additional information to resolve the sign of the square root.
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];

        // x² + y² = (R+rcos(v))² = R² + r²cos(v)² + 2rRcos(v)
        if (x*x + y*y < R_*R_ + r_*r_*cosv*cosv) {
            cosv = -cosv;
        }
        Real denom = r_*(R_ + r_*cosv);
#ifdef DEBUG
        if (denom == 0) {
            std::cerr << "r(R+rcos(v)) = " << denom << "\n";
            std::cerr << "cos(v) = " << cosv << ", p = " << p << "\n";
            std::cerr << "zr = " << zr << ", 1 - (z/r)^2 = " << 1 - zr*zr << "\n";
            std::cerr << "arg = " << arg << "\n";
        }
#endif
        return cosv/denom;
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
        residual *= 20*std::numeric_limits<Real>::epsilon();
        return residual;
    }

    // x = x₀ + (R + r cos(2πv)) cos(2πu)
    // y = y₀ + (R + r cos(2πv)) sin(2πu)
    // z = z₀ + r sin(2πv)
    std::pair<Real, Real> get_uv(vec<Real> const & p) const
    {
        using std::asin;
        using std::atan2;
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        Real z = p[2] - center_[2];
        Real u = atan2(y,x)/(2*M_PI);
        if(u < 0) {
            u += 1;
        }
        // asin(y) ∈ [-π/2, π/2].
        Real v = asin(z/r_);
        // asin cannot distinguish the quadrant. Use x and y to distinguish.
        // x² + y² = (R+rcos(2πv))² = R² + r²cos(2πv)² + 2rRcos(2πv)
        if (x*x + y*y < R_*R_) {
            v = M_PI - v;
        }
        // 3rd quadrant is the only place it can be:
        if (v < 0) {
            v += 2*M_PI;
        }

        // Now v ∈ [0, 2π]. We need it in [0,1]:
        v /= 2*M_PI;
        return std::pair<Real, Real>(u,v);
    }

    // x = x₀ + (R + r cos(2πv)) cos(2πu)
    // y = y₀ + (R + r cos(2πv)) sin(2πu)
    // z = z₀ + r sin(2πv)
    vec<Real,3> operator()(Real u, Real v) const override {
        vec<Real,3> w = center_;
        w[0] += (R_ + r_*cos(2*M_PI*v))*cos(2*M_PI*u);
        w[1] += (R_ + r_*cos(2*M_PI*v))*sin(2*M_PI*u);
        w[2] += r_*sin(2*M_PI*v);
        return w;
    }
    // PBRT uses invariants; e.g. for a circle the refinement is
    // pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    // This requires a direction in which to refine, so that we can update with Newton's method:
    // f(p + δtd) ≈ f(p) + δt𝝯f·d = 0 ⇒ δt = -f(p)/𝝯f·d
    void refine_hit_point(vec<Real> & p, vec<Real> const & direction) const {
        Real x = p[0] - center_[0];
        Real y = p[1] - center_[1];
        Real z = p[2] - center_[2];
        Real length_sq = x*x + y*y + z*z;
        Real tmp = length_sq + R_*R_ - r_*r_;
        Real fp = tmp*tmp - 4*R_*R_*(x*x + y*y);
        vec<Real> nablaf(4*x*(length_sq - R_*R_ - r_*r_),
                         4*y*(length_sq - R_*R_ - r_*r_),
                         4*z*(length_sq + R_*R_ - r_*r_));

        Real nablafdd = drt::dot(nablaf, direction);
        if (nablafdd != 0)
        {
            Real dt = -fp/nablafdd;
            p[0] += dt*direction[0];
            p[1] += dt*direction[1];
            p[2] += dt*direction[2];
        }
        // Keep nans away!
        p[0] = std::clamp(p[0], center_[0] - (R_+r_), center_[0] + (R_+r_));
        p[1] = std::clamp(p[1], center_[1] - (R_+r_), center_[1] + (R_+r_));
        p[2] = std::clamp(p[2], center_[2] - r_, center_[2] + r_);
    }

    virtual ~torus() = default;

private:

    void set_fundamental_forms(hit_record<Real>& rec) const
    {
        Real tmp = R_ + r_*std::cos(2*M_PI*rec.v);
        rec.E = 4*M_PI*M_PI*tmp*tmp;
        rec.F = 0;
        rec.G = 4*M_PI*M_PI*r_*r_;

        rec.e = -4*M_PI*M_PI*tmp*std::cos(2*M_PI*rec.v);
        rec.f = 0;
        rec.g = -4*M_PI*M_PI*r_;
    }
    vec<Real, 3> center_;
    Real R_; // major radius
    Real r_; // minor radius
};

template<typename Real>
bool torus<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    using std::abs;
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

    // Descartes rule of signs. This improves performance, but not by as much as one would hope.
    // Only about 7% of queries can achieve early exit.
    bool sign_changes = false;
    for (int64_t j = 3; j >= 0; --j) {
        if (c[j] < 0) {
            sign_changes = true;
            // The break introduces another branch which makes this slower:
            //break;
        }
    }
    if (!sign_changes) {
        return false;
    }

    auto roots = quartic_roots(c[4], c[3], c[2], c[1], c[0]);
    rec.t = std::numeric_limits<Real>::quiet_NaN();
    // Roots are sorted so this gets the minimal root-or if no intersection, t stays nan.
    for (auto root : roots) {
        if (root >= t_min && root <= t_max) {
            rec.t = root;
            break;
        }
    }
    if (std::isnan(rec.t)) {
        return false;
    }

    rec.p = r(rec.t);
    // PBRT often refines hit points.
    // The conditioning of the problem "find t such that f(o+td) = 0" is different (and often worse!)
    // than "find p such that f(p) = 0".
    this->refine_hit_point(rec.p, r.direction());

    Real x = rec.p[0] - center_[0];
    Real y = rec.p[1] - center_[1];
    Real z = rec.p[2] - center_[2];
    Real s = x*x + y*y + z*z - r_*r_;
    vec<Real> n(4*x*(s - R_*R_), 4*y*(s - R_*R_), 4*z*(s + R_*R_));
    rec.gradient_magnitude = norm(n);
    vec<Real> outward_normal = n/rec.gradient_magnitude;
    rec.set_face_normal(r, outward_normal);
    std::tie(rec.u, rec.v) = get_uv(rec.p);
    set_fundamental_forms(rec);
    Real res = this->residual(rec.p);
    Real expected_res = this->expected_residual(rec.p);
    if (abs(res) > expected_res) {
#ifdef DEBUG
        std::cerr << __FILE__ << ":" << __LINE__ << " Residual for torus intersection unexpectedly high. ";
        std::cerr << "Residual is " << res << ", but expected residual is " << expected_res << ".\n";
        std::cerr << rec << "\n";
        std::cerr << "[t_min, t_max] = [" << t_min << ", " << t_max << "]\n";
        std::cerr << "Ray: " << r << "\n";
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
