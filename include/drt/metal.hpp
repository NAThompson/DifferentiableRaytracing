#ifndef DRT_METAL_HPP
#define DRT_METAL_HPP

#include <drt/vec.hpp>
#include <drt/material.hpp>

namespace drt {

template<typename Real>
class metal : public material<Real> {
public:
    metal(const vec<Real>& a) : albedo_(a) {}

    virtual bool scatter(const ray<Real>& r_in, const hit_record<Real>& rec,
                         vec<Real>& attenuation, ray<Real>& scattered) override
    {
        auto in = r_in.direction();
        normalize(in);
        vec<Real> reflected = reflect(in, rec.normal);
        scattered = ray(rec.p, reflected);
        attenuation = albedo_;
        return (dot(scattered.direction(), rec.normal) > 0);
    }

    virtual ~metal() = default;

private:
    vec<Real> albedo_;
};

}
#endif