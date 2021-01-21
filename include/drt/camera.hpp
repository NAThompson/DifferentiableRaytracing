#ifndef DRT_CAMERA_HPP
#define DRT_CAMERA_HPP
#include <drt/ray.hpp>
#include <drt/vec.hpp>

namespace drt {

template<typename Real>
class camera {
public:
    camera() {
        Real aspect_ratio = 16.0 / 9.0;
        Real viewport_height = 2.0;
        Real viewport_width = aspect_ratio * viewport_height;
        Real focal_length = 1.0;

        origin_ = vec<Real, 3>(0, 0, 0);
        horizontal_ = vec<Real, 3>(viewport_width, 0.0, 0.0);
        vertical_ = vec<Real, 3>(0.0, viewport_height, 0.0);
        lower_left_corner_ = origin_ - horizontal_/Real(2) - vertical_/Real(2) - vec<Real, 3>(0, 0, focal_length);
    }

    ray<Real> get_ray(Real u, Real v) const {
        return ray<Real>(origin_, lower_left_corner_ + u*horizontal_ + v*vertical_ - origin_);
    }

private:
    vec<Real, 3> origin_;
    vec<Real, 3> lower_left_corner_;
    vec<Real, 3> horizontal_;
    vec<Real, 3> vertical_;
};
}
#endif