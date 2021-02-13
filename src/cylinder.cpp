#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/ray.hpp>
#include <drt/cylinder.hpp>
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
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>

using std::make_shared;
using namespace drt;

template<typename Real>
hittable_list<Real> cylinder_scene() {
    hittable_list<Real> objects;
    auto light = make_shared<diffuse_light<Real>>(vec<Real>(3, 3, 3));
    objects.add(make_shared<yz_rect<Real>>(-10, 10, -10, 10, 20, light));
    //objects.add(make_shared<xy_rect<Real>>(-100, 100, -100, 100, -20, light));

    Real radius = 3;
    Real z_min = -10;
    Real z_max = 10;

    auto mat = make_shared<lambertian<Real>>(vec<Real>(0,0,1));
    auto boundary = make_shared<cylinder<Real>>(radius, z_min, z_max, mat);
    objects.add(boundary);
    return objects;
}


int main() {
    using Real = double;

    const Real aspect_ratio = 1.6;
    const int64_t image_width = 1200;
    const int64_t image_height = static_cast<int64_t>(image_width/aspect_ratio);
    const int64_t samples_per_pixel = 16;

    auto world = cylinder_scene<Real>();
    Real scale = 1.0;
    drt::vec<Real> lookfrom(10, 10, 0);
    drt::vec<Real> lookat(0,0,0);
    drt::vec<Real> vup(0,1,0);

    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    drt::render_scene<Real>("cylinder.png", image_width, image_height, background, cam, world, samples_per_pixel);
}