#ifndef DRT_RAY_COLOR_HPP
#define DRT_RAY_COLOR_HPP
#include <drt/ray.hpp>
#include <drt/vec.hpp>
#include <drt/hittable.hpp>

namespace drt {

template<typename Real>
drt::vec<Real, 3> ray_color(const drt::ray<Real>& r, const drt::vec<Real>& background_color, const drt::hittable<Real> & world, int depth) {
    if (depth <= 0) {
        return drt::vec<Real>(0,0,0);
    }
    drt::hit_record<Real> rec;
    if (!world.hit(r, 1000*std::numeric_limits<Real>::epsilon(), std::numeric_limits<Real>::infinity(), rec)) {
        return background_color;
    }
    drt::ray<Real> scattered;
    drt::vec<Real> attenuation;
    drt::vec<Real> emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
        return emitted;
    }
    drt::vec<Real> color = ray_color(scattered, background_color, world, depth-1);
    // Peter Shirley has an overload for componentwise vector multiplication.
    // I'm avoiding that for now.
    // It does make sense for this application: Independently attenuating each frequency.
    for (size_t i = 0; i < 3; ++i) {
        color[i] *= attenuation[i];
    }
    return emitted + color;
}
}
#endif