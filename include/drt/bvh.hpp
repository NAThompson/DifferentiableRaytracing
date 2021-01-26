#ifndef DRT_BVH_H
#define DRT_BVH_H

#include <random>
#include <vector>
#include <memory>
#include <drt/aabb.hpp>
#include <drt/ray.hpp>
#include <drt/hittable.hpp>
#include <drt/hittable_list.hpp>

namespace drt {

template<typename Real>
class bvh_node : public hittable<Real> {
public:
    bvh_node();

    bvh_node(const hittable_list<Real>& list)
        : bvh_node(list.objects, 0, list.objects.size())
    {}

    bvh_node(const std::vector<std::shared_ptr<hittable<Real>>>& src_objects,
             size_t start, size_t end);

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~bvh_node() = default;

public:
    std::shared_ptr<hittable<Real>> left_;
    std::shared_ptr<hittable<Real>> right_;
    aabb<Real> box_;
};

template<typename Real>
bool bvh_node<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = box_;
    return true;
}

template<typename Real>
bool bvh_node<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    if (!box_.hit(r, t_min, t_max))
        return false;

    bool hit_left = left_->hit(r, t_min, t_max, rec);
    bool hit_right = right_->hit(r, t_min, hit_left ? rec.t : t_max, rec);

    return hit_left || hit_right;
}

template<typename Real>
inline bool box_compare(const std::shared_ptr<hittable<Real>> a, const std::shared_ptr<hittable<Real>> b, int axis)
{
    aabb<Real> box_a;
    aabb<Real> box_b;

    if (!a->bounding_box(box_a) || !b->bounding_box(box_b))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    return box_a.min()[axis] < box_b.min()[axis];
}

template<typename Real>
bool box_x_compare (const std::shared_ptr<hittable<Real>> a, const std::shared_ptr<hittable<Real>> b) {
    return box_compare(a, b, 0);
}

template<typename Real>
bool box_y_compare (const std::shared_ptr<hittable<Real>> a, const std::shared_ptr<hittable<Real>> b) {
    return box_compare(a, b, 1);
}

template<typename Real>
bool box_z_compare (const std::shared_ptr<hittable<Real>> a, const std::shared_ptr<hittable<Real>> b) {
    return box_compare(a, b, 2);
}

template<typename Real>
bvh_node<Real>::bvh_node(const std::vector<std::shared_ptr<hittable<Real>>>& src_objects,
                         size_t start, size_t end)
{
    auto objects = src_objects; // Create a modifiable array of the source scene objects

    std::uniform_int_distribution<int> dis(0,3);
    std::random_device rd;
    int axis = dis(rd);
    auto comparator = (axis == 0) ? box_x_compare<Real>
                    : (axis == 1) ? box_y_compare<Real>
                                  : box_z_compare<Real>;

    size_t object_span = end - start;

    if (object_span == 1)
    {
        left_ = right_ = objects[start];
    }
    else if (object_span == 2)
    {
        if (comparator(objects[start], objects[start+1]))
        {
            left_ = objects[start];
            right_ = objects[start+1];
        }
        else
        {
            left_ = objects[start+1];
            right_ = objects[start];
        }
    }
    else
    {
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span/2;
        left_ = std::make_shared<bvh_node<Real>>(objects, start, mid);
        right_ = std::make_shared<bvh_node<Real>>(objects, mid, end);
    }

    aabb<Real> box_left, box_right;

    if (  !left_->bounding_box(box_left) || !right_->bounding_box(box_right))
    {
        std::cerr << "No bounding box in bvh_node constructor.\n";
    }

    box_ = surrounding_box(box_left, box_right);
}

}
#endif