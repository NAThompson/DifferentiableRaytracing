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

    virtual vec<Real> emitted(hit_record<Real> const & hr) const override
    {
        return emit_->value(hr);
    }

    virtual ~diffuse_light() = default;

public:
    std::shared_ptr<texture<Real>> emit_;
};

}
#endif