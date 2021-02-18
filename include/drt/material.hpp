#ifndef DRT_MATERIAL_HPP
#define DRT_MATERIAL_HPP
#include <drt/ray.hpp>
#include <drt/vec.hpp>
#include <drt/hittable.hpp>


namespace drt {

template<typename Real>
class material {
public:
    virtual bool scatter(const ray<Real>& r_in, const hit_record<Real>& rec,
                         vec<Real, 3>& attenuation, ray<Real>& scattered) = 0;

    virtual vec<Real> emitted([[maybe_unused]] hit_record<Real> const &) const {
        return vec<Real>(0,0,0);
    }
};

}
#endif