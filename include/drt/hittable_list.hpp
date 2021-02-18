#ifndef DRT_HITTABLE_LIST_H
#define DRT_HITTABLE_LIST_H

#include <drt/hittable.hpp>
#include <drt/material.hpp>
#include <memory>
#include <vector>

namespace drt {

template<typename Real>
class hittable_list {
public:
    hittable_list() {}
    hittable_list(std::shared_ptr<hittable<Real>> object, std::shared_ptr<material<Real>> mat) { add(object, mat); }

    void clear() {
        objects.clear();
        materials.clear();
    }
    void add(std::shared_ptr<hittable<Real>> object, std::shared_ptr<material<Real>> mat) {
        objects.push_back(object);
        materials.push_back(mat);
    }

    std::shared_ptr<material<Real>> hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const;

    bool bounding_box(aabb<Real>& output_box) const;

public:
    std::vector<std::shared_ptr<hittable<Real>>> objects;
    std::vector<std::shared_ptr<material<Real>>> materials;
};

template<typename Real>
std::shared_ptr<material<Real>> hittable_list<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    hit_record<Real> temp_rec;
    auto closest_so_far = t_max;
    int64_t hit_idx = -1;
    for (size_t i = 0; i < objects.size(); ++i) {
        auto object = objects[i];
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            closest_so_far = temp_rec.t;
            rec = temp_rec;
            hit_idx = i;
        }
    }
    if (hit_idx >= 0)
    {
        return materials[hit_idx];
    }
    else
    {
        return nullptr;
    }
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