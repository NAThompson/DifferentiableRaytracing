#ifndef DRT_CAMERA_HPP
#define DRT_CAMERA_HPP
#include <drt/ray.hpp>
#include <drt/vec.hpp>

namespace drt {

template<typename Real>
class camera {
public:
    camera(vec<Real> lookfrom,
           vec<Real> lookat,
           vec<Real> view_up,
           Real vfov, // vertical field-of-view in degrees
           Real aspect_ratio = 16.0/9.0) : origin_(lookfrom), lookat_(lookat) {
        Real theta = vfov*M_PI/180.0;
        Real h = std::tan(theta/2);
        Real viewport_height = 2.0 * h;
        Real viewport_width = aspect_ratio * viewport_height;

        auto w = lookfrom - lookat;
        normalize(w);
        auto u = cross(view_up, w);
        normalize(u);
        auto v = cross(w, u);

        horizontal_ = viewport_width*u;
        vertical_ = viewport_height*v;
        lower_left_corner_ = origin_ - horizontal_/Real(2) - vertical_/Real(2) - w;
    }

    ray<Real> get_ray(Real u, Real v) const {
        return ray<Real>(origin_, lower_left_corner_ + u*horizontal_ + v*vertical_ - origin_);
    }

    std::pair<vec<Real>,vec<Real>> backlight() const {
        vec<Real,3> o = Real(2)*origin_ - lookat_;
        vec<Real,3> d = lookat_ - origin_;
        normalize(d);
        return std::make_pair(o,d);
    }

private:
    vec<Real, 3> origin_;
    vec<Real, 3> lower_left_corner_;
    vec<Real, 3> horizontal_;
    vec<Real, 3> vertical_;
    vec<Real,3> lookat_;
};
}
#endif
