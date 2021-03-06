#ifndef DRT_ISOTROPIC_HPP
#define DRT_ISOTROPIC_HPP
#include <random>
#include <drt/material.hpp>
#include <drt/texture.hpp>
#include <drt/ray.hpp>
#include <drt/vec.hpp>

namespace drt {

template<typename Real>
class isotropic : public material<Real> {
public:
    isotropic(vec<Real> c) : isotropic(std::make_shared<solid_color<Real>>(c)) {}

    isotropic(std::shared_ptr<texture<Real>> a) : albedo_(a) {
        dis_ = std::normal_distribution<Real>(0,1);
    }

    virtual bool scatter([[maybe_unused]] const ray<Real>& r_in, const hit_record<Real>& rec,
                         vec<Real>& attenuation, ray<Real>& scattered) override
    {

        scattered = ray(rec.p, random_in_unit_sphere());
        attenuation = albedo_->value(rec);
        return true;
    }

    virtual ~isotropic() = default;

private:
    std::shared_ptr<texture<Real>> albedo_;
    std::normal_distribution<Real> dis_;
    std::mt19937_64 gen_;

    vec<Real> random_in_unit_sphere() {
        vec<Real> v;
        v[0] = dis_(gen_);
        v[1] = dis_(gen_);
        v[2] = dis_(gen_);
        Real lambda = norm(v);
        v /= lambda;
        return v;
    }
};

}

#endif
