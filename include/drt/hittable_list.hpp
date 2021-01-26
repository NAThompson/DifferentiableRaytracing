#ifndef DRT_HITTABLE_LIST_H
#define DRT_HITTABLE_LIST_H

#include <drt/hittable.hpp>

#include <memory>
#include <vector>

namespace drt {

template<typename Real>
class hittable_list : public hittable<Real> {
public:
    hittable_list() {}
    hittable_list(std::shared_ptr<hittable<Real>> object) { add(object); }

    void clear() { objects.clear(); }
    void add(std::shared_ptr<hittable<Real>> object) { objects.push_back(object); }

    virtual bool hit(
        const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;
public:
    std::vector<std::shared_ptr<hittable<Real>>> objects;
};

template<typename Real>
bool hittable_list<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    hit_record<Real> temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

template<typename Real>
bool hittable_list<Real>::bounding_box(aabb<Real>& output_box) const {
    if (objects.empty()) return false;

    aabb<Real> temp_box;
    bool first_box = true;
    for (const auto& object : objects) {
        if (!object->bounding_box(temp_box)) return false;
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }

    return true;
}
}
#endif