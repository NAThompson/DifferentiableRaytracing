#ifndef DRT_ELLIPSOID_HPP
#define DRT_ELLIPSOID_HPP
#include <utility>
#include <algorithm>
#include <drt/hittable.hpp>
#include <drt/roots.hpp>
#include <drt/vec.hpp>

namespace drt {
template<typename Real>
class ellipsoid : public hittable<Real> {
public:
    ellipsoid(vec<Real, 3> const & center, vec<Real, 3> const & scales)
       : center_(center), scales_(scales)
    {
        if (scales_[0] <= 0) {
            std::cerr << "a <= 0 is not allowed for an ellipsoid.\n";
        }
        if (scales_[1] <= 0) {
            std::cerr << "b <= 0 is not allowed for an ellipsoid.\n";
        }
        if (scales_[2] <= 0) {
            std::cerr << "c <= 0 is not allowed for an ellipsoid.\n";
        }
    };

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    // The parametrization is:
    // σ(u,v) = (x₀ + a·cos(2πu)sin(πv), y₀ + b·sin(2πu)sin(πv), z₀ + c·cos(πv)), u,v ∈ [0,1].
    std::pair<Real, Real> get_uv(const vec<Real>& p) const {
        Real x = (p[0] - center_[0])/scales_[0];
        Real y = (p[1] - center_[1])/scales_[1];
        Real z = (p[2] - center_[2])/scales_[2];

        Real v = std::acos(z)/M_PI;
        if (v < 0) {
            v += 1;
        }
        Real u = std::atan2(y, x)/(2*M_PI);
        if (u < 0) {
            u += 1;
        }
        return std::make_pair(u, v);
    }

    // σ(u,v) = (x₀ + a·cos(2πu)sin(πv), y₀ + b·sin(2πu)sin(πv), z₀ + c·cos(πv)), u,v ∈ [0,1].
    vec<Real, 3> operator()(Real u, Real v) const override {
        vec<Real,3> w = center_;
        w[0] += scales_[0]*cos(2*M_PI*u)*sin(M_PI*v);
        w[1] += scales_[1]*sin(2*M_PI*u)*sin(M_PI*v);
        w[2] += scales_[2]*cos(M_PI*v);
        return w;
    }

    Real gaussian_curvature(const vec<Real>& p) {
        Real asq = scales_[0]*scales_[0];
        Real bsq = scales_[1]*scales_[1];
        Real csq = scales_[2]*scales_[2];
        Real ysq = (p[1] - center_[1])*(p[1] - center_[1]);
        Real zsq = (p[2] - center_[2])*(p[2] - center_[2]);
        // https://mathworld.wolfram.com/Ellipsoid.html, equation 14:
        Real numerator = asq*bsq*bsq*bsq*csq*csq*csq;
        Real sqrt_denom = csq*csq*bsq*bsq + csq*csq*(asq - bsq)*ysq + bsq*bsq*(asq-csq)*zsq;
        return numerator/(sqrt_denom*sqrt_denom);
    }

    std::pair<Real, Real> gaussian_curvature_bounds() const {
        Real asq = scales_[0]*scales_[0];
        Real bsq = scales_[1]*scales_[1];
        Real csq = scales_[2]*scales_[2];

        Real kappa_min = std::min({csq/(asq*bsq), asq/(csq*bsq), bsq/(asq*csq)});
        Real kappa_max = std::max({csq/(asq*bsq), asq/(csq*bsq), bsq/(asq*csq)});
        return std::make_pair(kappa_min, kappa_max);
    }

    virtual ~ellipsoid() = default;

private:

    // Private since rec.u, rec.v must be set before this can be called:
    void set_fundamental_forms(hit_record<Real>& rec) const {
        using std::sin;
        using std::cos;
        using std::sqrt;
        // First fundamental form:
        Real a = scales_[0];
        Real b = scales_[1];
        Real c = scales_[2];
        Real spv = sin(M_PI*rec.v);
        Real cpv = cos(M_PI*rec.v);
        Real stpu = sin(2*M_PI*rec.u);
        Real ctpu = cos(2*M_PI*rec.u);
        rec.E = 4*M_PI*M_PI*spv*spv*(a*a*stpu*stpu + b*b*ctpu*ctpu);
        rec.F = 2*M_PI*M_PI*ctpu*stpu*cpv*spv*(b*b-a*a);
        rec.G = M_PI*M_PI*(a*a*cpv*cpv*ctpu*ctpu + b*b*stpu*stpu*cpv*cpv + c*c*spv*spv);

        // Second fundamental form.
        // Computed from rescaling the information here: https://mathworld.wolfram.com/Ellipsoid.html
        Real denom = sqrt(a*a*b*b*cpv*cpv + c*c*(b*b*ctpu*ctpu + a*a*stpu*stpu)*spv*spv);
        rec.e = 4*M_PI*M_PI*a*b*c*spv*spv/denom;
        rec.f = 0;
        rec.g = M_PI*M_PI*a*b*c/denom;
    }

    vec<Real, 3> center_;
    vec<Real, 3> scales_;
};

template<typename Real>
bool ellipsoid<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> oc = r.origin() - center_;
    oc[0] /= scales_[0];
    oc[1] /= scales_[1];
    oc[2] /= scales_[2];
    vec<Real> sd = r.direction();
    sd[0] /= scales_[0];
    sd[1] /= scales_[1];
    sd[2] /= scales_[2];

    Real a = drt::squared_norm(sd);
    Real b = 2*drt::dot(oc, sd);
    Real c = drt::squared_norm(oc) - 1;
    auto opt_root = first_quadratic_root_in_range(a, b, c, t_min, t_max);
    if (!opt_root) {
        return false;
    }
    rec.t = *opt_root;
    rec.p = r(rec.t);
    vec<Real> outward_normal = (rec.p - center_);
    outward_normal[0] *= 2/(scales_[0]*scales_[0]);
    outward_normal[1] *= 2/(scales_[1]*scales_[1]);
    outward_normal[2] *= 2/(scales_[2]*scales_[2]);
    rec.gradient_magnitude = norm(outward_normal);
    outward_normal = outward_normal/rec.gradient_magnitude;
    std::tie(rec.u, rec.v) = get_uv(rec.p);
    set_fundamental_forms(rec);
    rec.set_face_normal(r, outward_normal);
    return true;
}

template<typename Real>
bool ellipsoid<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(center_ - scales_, center_ + scales_);
    return true;
}

}
#endif
