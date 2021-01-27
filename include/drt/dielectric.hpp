#ifndef DRT_DIELECTRIC_HPP
#define DRT_DIELECTRIC_HPP

#include <cmath>
#include <drt/ray.hpp>
#include <drt/material.hpp>

namespace drt {

template<typename Real>
vec<Real> refract(const vec<Real>& uv, const vec<Real>& n, Real etai_over_etat) {
    using std::sqrt;
    Real cos_theta = std::min(dot(-uv, n), Real(1));
    vec<Real> r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec<Real> r_out_parallel = -sqrt(fabs(Real(1) - squared_norm(r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}

template<typename Real>
class dielectric : public material<Real> {
public:
    dielectric(Real index_of_refraction) : ir_(index_of_refraction) {
        dis_ = std::uniform_real_distribution<Real>(0,1);
    }

    virtual bool scatter(const ray<Real>& r_in, const hit_record<Real>& rec,
                         vec<Real>& attenuation, ray<Real>& scattered) override
    {
        using std::sqrt;
        attenuation = vec<Real>(1.0, 1.0, 1.0);
        Real refraction_ratio = rec.front_face ? (1/ir_) : ir_;

        vec<Real> unit_direction = r_in.direction();
        normalize(unit_direction);

        Real cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        Real sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec<Real> direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > dis_(gen_))
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered = ray(rec.p, direction);
        return true;
    }

    virtual ~dielectric() = default;

private:

    static Real reflectance(Real cosine, Real ref_idx)
    {
        using std::pow;
        // Use Schlick's approximation for reflectance.
        auto r0 = (1-ref_idx) / (1+ref_idx);
        r0 = r0*r0;
        return r0 + (1-r0)*pow((1 - cosine), 5);
    }

    Real ir_; // Index of Refraction
    std::uniform_real_distribution<Real> dis_;
    std::mt19937 gen_;
};

}
#endif