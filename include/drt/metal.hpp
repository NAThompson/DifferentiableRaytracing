#ifndef DRT_METAL_HPP
#define DRT_METAL_HPP

#include <drt/vec.hpp>
#include <drt/material.hpp>

namespace drt {

template<typename Real>
class metal : public material<Real> {
public:

    metal(std::shared_ptr<texture<Real>> a) : albedo_(a) {
    }

    metal(const vec<Real, 3> & color) : metal(std::make_shared<solid_color<Real>>(color)) {}

    virtual bool scatter(const ray<Real>& r_in, const hit_record<Real>& rec,
                         vec<Real>& attenuation, ray<Real>& scattered) override
    {
        auto in = r_in.direction();
        normalize(in);
        vec<Real> reflected = reflect(in, rec.normal);
        scattered = ray(rec.p, reflected);
        attenuation = albedo_->value(rec.u, rec.v, rec.p);
        return (dot(scattered.direction(), rec.normal) > 0);
    }

    virtual ~metal() = default;

private:
    std::shared_ptr<texture<Real>> albedo_;
};

}
#endif