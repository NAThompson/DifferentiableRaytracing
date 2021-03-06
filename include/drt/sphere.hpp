#ifndef DRT_SPHERE_HPP
#define DRT_SPHERE_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/roots.hpp>

namespace drt {
template<typename Real>
class sphere : public hittable<Real> {
public:
    sphere() {}
    sphere(vec<Real, 3> const & center, Real radius)
       : center_(center), radius_(radius)
    {};

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~sphere() = default;

    // The parametrization is:
    // σ(u,v) = (x₀ + rcos(2πu)sin(πv), y₀ + rsin(2πu)sin(πv), z₀ + rcos(πv)), u,v ∈ [0,1].
    std::pair<Real, Real> get_uv(const vec<Real>& p) const {
        Real x = (p[0] - center_[0])/radius_;
        Real y = (p[1] - center_[1])/radius_;
        Real z = (p[2] - center_[2])/radius_;

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

    Real area() const {
        return 4*M_PI*radius_*radius_;
    }

    Real volume() const {
        return 4*M_PI*radius_*radius_*radius_/3;
    }

    // σ(u,v) = (x₀ + rcos(2πu)sin(πv), y₀ + rsin(2πu)sin(πv), z₀ + rcos(πv)), u,v ∈ [0,1].
    vec<Real,3> operator()(Real u, Real v) const override {
        vec<Real,3> w = center_;
        w[0] += radius_*cos(2*M_PI*u)*sin(M_PI*v);
        w[1] += radius_*sin(2*M_PI*u)*sin(M_PI*v);
        w[2] += radius_*cos(M_PI*v);
        return w;
    }

    vec<Real,3> operator()(Real u, Real v, matrix<Real, 3,3>& J) const override {
        vec<Real,3> w = center_;
        w[0] += radius_*cos(2*M_PI*u)*sin(M_PI*v);
        w[1] += radius_*sin(2*M_PI*u)*sin(M_PI*v);
        w[2] += radius_*cos(M_PI*v);
        // J(i,j) = ∂ⱼσᵢ.
        // J(i,0) = ∂ᵤσᵢ.
        J(0,0) = -2*M_PI*radius_*sin(2*M_PI*u)*sin(M_PI*v);
        J(1,0) = 2*M_PI*radius_*cos(2*M_PI*u)*sin(M_PI*v);
        J(2,0) = 0;

        // J(i,1) = ∂ᵥσᵢ.
        J(0,1) = M_PI*radius_*cos(2*M_PI*u)*cos(M_PI*v);
        J(1,1) = M_PI*radius_*sin(2*M_PI*u)*cos(M_PI*v);
        J(2,1) = -M_PI*radius_*sin(M_PI*v);

        J(0,2) = std::numeric_limits<Real>::quiet_NaN();
        J(1,2) = std::numeric_limits<Real>::quiet_NaN();
        J(2,2) = std::numeric_limits<Real>::quiet_NaN();

        return w;
    }

    vec<Real,3> operator()(Real u, Real v, matrix<Real, 3,3>& J, tensor<Real, 3,2,2>& H) const override {
        vec<Real,3> w = center_;
        w[0] += radius_*cos(2*M_PI*u)*sin(M_PI*v);
        w[1] += radius_*sin(2*M_PI*u)*sin(M_PI*v);
        w[2] += radius_*cos(M_PI*v);
        // J(i,j) = ∂ⱼσᵢ.
        // J(i,0) = ∂ᵤσᵢ.
        J(0,0) = -2*M_PI*radius_*sin(2*M_PI*u)*sin(M_PI*v);
        J(1,0) = 2*M_PI*radius_*cos(2*M_PI*u)*sin(M_PI*v);
        J(2,0) = 0;

        // J(i,1) = ∂ᵥσᵢ.
        J(0,1) = M_PI*radius_*cos(2*M_PI*u)*cos(M_PI*v);
        J(1,1) = M_PI*radius_*sin(2*M_PI*u)*cos(M_PI*v);
        J(2,1) = -M_PI*radius_*sin(M_PI*v);

        J(0,2) = std::numeric_limits<Real>::quiet_NaN();
        J(1,2) = std::numeric_limits<Real>::quiet_NaN();
        J(2,2) = std::numeric_limits<Real>::quiet_NaN();

        // H(i,j,k) = ∂ⱼ∂ₖσᵢ; in this case u is the 0 index, v is the 1 index.
        // H(0,0,0) = ∂ᵤ∂ᵤσ₀
        H(0,0,0) = -4*M_PI*M_PI*radius_*cos(2*M_PI*u)*sin(M_PI*v);
        H(1,0,0) = -4*M_PI*M_PI*radius_*sin(2*M_PI*u)*sin(M_PI*v);
        H(2,0,0) = 0;

        // Mixed partials:
        H(0,0,1) = -2*M_PI*M_PI*radius_*sin(2*M_PI*u)*cos(M_PI*v);
        H(0,1,0) = H(0,0,1);

        H(1,0,1) = 2*M_PI*M_PI*radius_*cos(2*M_PI*u)*cos(M_PI*v);
        H(1,1,0) = H(1,0,1);

        H(2,0,1) = 0;
        H(2,1,0) = 0;

        // ∂ᵥ∂ᵥσ₀
        H(0,1,1) = -M_PI*M_PI*radius_*cos(2*M_PI*u)*sin(M_PI*v);
        H(1,1,1) = -M_PI*M_PI*radius_*sin(2*M_PI*u)*sin(M_PI*v);
        H(2,1,1) = -M_PI*M_PI*radius_*cos(M_PI*v);
        return w;
    }

private:
    // Private since rec.u, rec.v must be set before this can be called:
    void set_fundamental_forms(hit_record<Real>& rec) const {
        // First fundamental form:
        Real dsigmadu = 2*M_PI*radius_*std::sin(M_PI*rec.v);
        rec.E = dsigmadu*dsigmadu;
        rec.F = 0;
        rec.G = M_PI*M_PI*radius_*radius_;

        // The second fundamental form of a sphere is the same as the first,
        // see Pressley, Elementary Differential Geometry, Section 8.1.
        // However, a variable radius changes this a bit, see https://mathworld.wolfram.com/Sphere.html
        rec.e = rec.E/radius_;
        rec.f = 0;
        rec.g = M_PI*M_PI*radius_;
    }

    vec<Real, 3> center_;
    Real radius_;
};

template<typename Real>
bool sphere<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> oc = r.origin() - center_;
    Real a = drt::squared_norm(r.direction());
    Real b = 2*drt::dot(oc, r.direction());
    Real c = drt::squared_norm(oc) - radius_*radius_;
    auto opt_root = first_quadratic_root_in_range(a, b, c, t_min, t_max);
    if (!opt_root) {
        return false;
    }
    rec.t = *opt_root;
    rec.p = r(rec.t);
    auto outward_normal = (rec.p - center_) / radius_;
    rec.set_face_normal(r, outward_normal);
    std::tie(rec.u, rec.v) = get_uv(rec.p);
    rec.gradient_magnitude = 2*radius_;
    set_fundamental_forms(rec);
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
