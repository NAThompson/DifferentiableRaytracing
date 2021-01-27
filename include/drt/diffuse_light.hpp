#ifndef DRT_DIFFUSE_LIGHT_HPP
#define DRT_DIFFUSE_LIGHT_HPP

#include <drt/vec.hpp>
#include <drt/texture.hpp>
#include <drt/material.hpp>

namespace drt {

template<typename Real>
class diffuse_light : public material<Real>  {
public:
    diffuse_light(std::shared_ptr<texture<Real>> a) : emit_(a) {}
    diffuse_light(vec<Real> c) : emit_(std::make_shared<solid_color<Real>>(c)) {}

    virtual bool scatter([[maybe_unused]] const ray<Real>& r_in, [[maybe_unused]] const hit_record<Real>& rec,
                         [[maybe_unused]] vec<Real>& attenuation, [[maybe_unused]] ray<Real>& scattered) override
    {
        return false;
    }

    virtual vec<Real> emitted(Real u, Real v, const vec<Real>& p) const override
    {
        return emit_->value(u, v, p);
    }

public:
    std::shared_ptr<texture<Real>> emit_;
};

}
#endif