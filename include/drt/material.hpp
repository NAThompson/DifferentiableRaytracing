#ifndef DRT_MATERIAL_HPP
#define DRT_MATERIAL_HPP
#include <drt/ray.hpp>

struct hit_record;


namespace drt {

template<typename Real>
class material {
public:
    virtual bool scatter(const ray<Real>& r_in, const hit_record<Real>& rec,
                         vec<Real, 3>& attenuation, ray<Real>& scattered) = 0;
};

}
#endif