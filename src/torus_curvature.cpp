#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/ray.hpp>
#include <drt/sphere.hpp>
#include <drt/hittable.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>
#include <drt/lambertian.hpp>
#include <drt/dielectric.hpp>
#include <drt/metal.hpp>
#include <drt/aabb.hpp>
#include <drt/bvh.hpp>
#include <drt/texture.hpp>
#include <drt/diffuse_light.hpp>
#include <drt/aarect.hpp>
#include <drt/color_maps.hpp>
#include <drt/torus.hpp>
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>

using std::make_shared;
using namespace drt;

template<typename Real>
hittable_list<Real> torus_scene() {
    hittable_list<Real> objects;
    auto light = make_shared<diffuse_light<Real>>(vec<Real>(1, 1, 1));
    objects.add(make_shared<yz_rect<Real>>(-10, 10, -10, 10, 20, light));
    objects.add(make_shared<xy_rect<Real>>(-100, 100, -100, 100, -20, light));

    auto dummy_mat = make_shared<lambertian<Real>>(vec<Real,3>(0,0,0));
    vec<Real> center(0,0,0);
    Real R = 3.0;
    Real r = 1.0;
    torus<Real> tor(center, R, r, dummy_mat);
    std::function<vec<Real>(Real, Real, const vec<Real> &)> gaussian_curvature = [=]([[maybe_unused]] Real u, [[maybe_unused]] Real v, vec<Real> const & p) {
        auto [kappa_min, kappa_max] = tor.gaussian_curvature_bounds();
        Real kappa = tor.gaussian_curvature(p);
        if (std::isnan(kappa)) {
            return vec<Real, 3>(0,0,0);
        }
        Real scalar = (kappa - kappa_min)/(kappa_max - kappa_min);
        if (scalar < 0) {
            std::cerr << "Scalar = " << scalar << ", this bad.\n";
        }
        if (scalar > 1) {
            std::cerr << "Scalar = " << scalar << ", this bad.\n";
        }
        return viridis(scalar);
    };

    auto ptr = make_shared<decltype(gaussian_curvature)>(gaussian_curvature);
    auto ltext = drt::lambda_texture(ptr);
    auto ltext_ptr = make_shared<decltype(ltext)>(ltext);

    auto mat = make_shared<lambertian<Real>>(ltext_ptr);
    auto boundary = make_shared<torus<Real>>(center, R, r, mat);
    objects.add(boundary);
    return objects;
}


int main() {
    using Real = double;

    const Real aspect_ratio = 1.6;
    const int64_t image_width = 1200;
    const int64_t image_height = static_cast<int64_t>(image_width/aspect_ratio);
    const int64_t samples_per_pixel = 512;

    auto world = torus_scene<Real>();
    Real scale = 0.7;
    drt::vec<Real> lookfrom(scale*5, -10*scale, -15*scale);
    drt::vec<Real> lookat(0,0,0);
    drt::vec<Real> vup(0,1,0);

    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    drt::render_scene<Real>("torus_curvature.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
