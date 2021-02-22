#ifndef DRT_AABB_HPP
#define DRT_AABB_HPP

#include <algorithm>
#include <drt/vec.hpp>

namespace drt {

template<typename Real>
class aabb {
public:
    aabb() {}
    aabb(const vec<Real>& a, const vec<Real>& b) { min_ = a; max_ = b;}

    vec<Real> min() const {return min_; }
    vec<Real> max() const {return max_; }

    bool hit(const ray<Real>& r, Real t_min, Real t_max) const {
        for (int64_t i = 0; i < 3; i++) {
            auto invD = Real(1) / r.direction()[i];
            auto t0 = (min_[i] - r.origin()[i]) * invD;
            auto t1 = (max_[i] - r.origin()[i]) * invD;
            if (invD < 0)
                std::swap(t0, t1);
            t_min = t0 > t_min ? t0 : t_min;
            t_max = t1 < t_max ? t1 : t_max;
            if (t_max <= t_min)
                return false;
        }
        return true;
    }

    friend std::ostream& operator<<(std::ostream & os, const aabb<Real> & v)
    {
        os << "{" << v.min_ << ", " << v.max_ << "}";
        return os;
    }

    vec<Real> min_;
    vec<Real> max_;
};

template<typename Real>
aabb<Real> surrounding_box(const aabb<Real>&  box0, const aabb<Real>& box1)
{
    
    vec<Real> small(fmin(box0.min()[0], box1.min()[0]),
                    fmin(box0.min()[1], box1.min()[1]),
                    fmin(box0.min()[2], box1.min()[2]));

    vec<Real> big(fmax(box0.max()[0], box1.max()[0]),
                  fmax(box0.max()[1], box1.max()[1]),
                  fmax(box0.max()[2], box1.max()[2]));

    return aabb(small,big);
}
}
#endif
