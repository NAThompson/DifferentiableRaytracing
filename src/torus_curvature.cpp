#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/png.hpp>
#include <drt/ray.hpp>
#include <drt/sphere.hpp>
#include <drt/hittable.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>
#include <drt/material.hpp>
#include <drt/lambertian.hpp>
#include <drt/dielectric.hpp>
#include <drt/metal.hpp>
#include <drt/aabb.hpp>
#include <drt/bvh.hpp>
#include <drt/texture.hpp>
#include <drt/checker_texture.hpp>
#include <drt/diffuse_light.hpp>
#include <drt/aarect.hpp>
#include <drt/constant_medium.hpp>
#include <drt/box.hpp>
#include <drt/ellipsoid.hpp>
#include <drt/color_maps.hpp>
#include <drt/torus.hpp>
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>

using std::make_shared;
using namespace drt;

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 1200;
    const int64_t image_height = static_cast<int64_t>(image_width);
    const int64_t samples_per_pixel = 8;

    auto mat = make_shared<lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> center(0,0,0);
    auto world = torus<Real>(center, 3.0, 1.0, mat);

    drt::vec<Real> lookfrom(0, -10, -15);
    drt::vec<Real> lookat = center;
    drt::vec<Real> vup(0,1,0);

    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    drt::vec<Real> background(1.0, 1.0, 1.0);
    drt::render_scene<Real>("torus_curvature.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
