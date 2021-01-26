#ifndef DRT_LAMBERTIAN_HPP
#define DRT_LAMBERTIAN_HPP

#include <random>
#include <drt/material.hpp>
#include <drt/ray.hpp>

namespace drt {

template<typename Real>
class lambertian : public material<Real> {
public:
    lambertian(vec<Real, 3> const & color) : albedo_(color) {
        dis_ = std::uniform_real_distribution<Real>(-1,1);
    }

    virtual bool scatter([[maybe_unused]] const ray<Real>& r_in, const hit_record<Real>& rec,
                         vec<Real, 3>& attenuation, ray<Real>& scattered) override
    {
        auto scatter_direction = rec.normal + random_in_hemisphere(rec.normal);
        if (squared_norm(scatter_direction) < sqrt(std::numeric_limits<Real>::epsilon())) {
            scatter_direction = rec.normal;
        }

        scattered = ray(rec.p, scatter_direction);
        attenuation = albedo_;
        return true;
    }

    vec<Real, 3> random_unit_vector() {
        vec<Real, 3> v;
        do {
            v[0] = dis_(gen_);
            v[1] = dis_(gen_);
            v[2] = dis_(gen_);
        } while (squared_norm(v) > 1);
        return v;
    }

    vec<Real> random_in_hemisphere(const vec<Real>& normal) {
        vec<Real> in_unit_sphere = random_unit_vector();
        if (dot(in_unit_sphere, normal) > 0) // In the same hemisphere as the normal
            return in_unit_sphere;
        else
            return -in_unit_sphere;
    }

    virtual ~lambertian() = default;

private:
    vec<Real, 3> albedo_;
    std::uniform_real_distribution<Real> dis_;
    std::mt19937 gen_;
};

}
#endif