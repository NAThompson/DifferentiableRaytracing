#ifndef DRT_HITTABLE_HPP
#define DRT_HITTABLE_HPP

#include <drt/ray.hpp>
#include <drt/vec.hpp>
#include <drt/aabb.hpp>


namespace drt {

template<typename Real>
class material;

template<typename Real>
struct hit_record {
    vec<Real, 3> p;
    vec<Real, 3> normal;
    Real t;
    bool front_face;
    std::shared_ptr<material<Real>> mat_ptr;

    inline void set_face_normal(const ray<Real>& r, const vec<Real, 3>& outward_normal) {
        front_face = drt::dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};

template<typename Real>
class hittable {
public:
    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const = 0;

    virtual bool bounding_box(aabb<Real>& output_box) const = 0;
};

}
#endif