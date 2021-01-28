#ifndef DRT_BOX_HPP
#define DRT_BOX_HPP
#include <drt/vec.hpp>
#include <drt/aarect.hpp>
#include <drt/hittable_list.hpp>

namespace drt {

template<typename Real>
class box : public hittable<Real>  {
public:
    box() {}
    box(const vec<Real>& p0, const vec<Real>& p1, std::shared_ptr<material<Real>> ptr);

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override {
        output_box = aabb(box_min_, box_max_);
        return true;
    }

    virtual ~box() = default;

private:
    vec<Real> box_min_;
    vec<Real> box_max_;
    hittable_list<Real> sides_;
};

template<typename Real>
box<Real>::box(const vec<Real>& p0, const vec<Real>& p1, std::shared_ptr<material<Real>> ptr) {
    using std::make_shared;
    box_min_ = p0;
    box_max_ = p1;

    sides_.add(make_shared<xy_rect<Real>>(p0[0], p1[0], p0[1], p1[1], p1[2], ptr));
    sides_.add(make_shared<xy_rect<Real>>(p0[0], p1[0], p0[1], p1[1], p0[2], ptr));

    sides_.add(make_shared<xz_rect<Real>>(p0[0], p1[0], p0[2], p1[2], p1[1], ptr));
    sides_.add(make_shared<xz_rect<Real>>(p0[0], p1[0], p0[2], p1[2], p0[1], ptr));

    sides_.add(make_shared<yz_rect<Real>>(p0[1], p1[1], p0[2], p1[2], p1[0], ptr));
    sides_.add(make_shared<yz_rect<Real>>(p0[1], p1[1], p0[2], p1[2], p0[0], ptr));
}

template<typename Real>
bool box<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    return sides_.hit(r, t_min, t_max, rec);
}


}
#endif