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
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>

template<typename Real>
drt::hittable_list<Real> simple_light() {
    using std::make_shared;
    using drt::sphere;
    using drt::vec;
    using drt::lambertian;
    using drt::diffuse_light;
    using drt::constant_medium;
    drt::hittable_list<Real> objects;

    auto chtext = make_shared<drt::checker_texture<Real>>(vec<Real>(0.8, 0.0, 1.0), vec<Real>(0.0, 0.8, 1.0));
    objects.add(make_shared<sphere<Real>>(vec<Real>(0,-1000,0), Real(1000), make_shared<drt::lambertian<Real>>(chtext)));

    auto stext = make_shared<drt::solid_color<Real>>(vec<Real>(0.2, 0.4, 0.9));
    objects.add(make_shared<sphere<Real>>(vec<Real>(0,2,0), Real(2), make_shared<drt::metal<Real>>(stext)));

    auto difflight = make_shared<diffuse_light<Real>>(vec<Real>(4,4,4));
    objects.add(make_shared<drt::xy_rect<Real>>(Real(3), Real(5), Real(1), Real(3), Real(-2), difflight));

    return objects;
}

int main() {
    using Real = float;

    const Real aspect_ratio = 1.6;
    const int64_t image_width = 1800;
    const int64_t image_height = static_cast<int>(image_width / aspect_ratio);
    const int64_t samples_per_pixel = 16;

    drt::hittable_list<Real> world1 = simple_light<Real>();
    drt::bvh_node<Real> world(world1);

    drt::vec<Real> lookfrom(26,3,6);
    drt::vec<Real> lookat(0,2,0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(20), aspect_ratio);
    drt::vec<Real> background(0.01, 0.02, 0.01);
    drt::render_scene<Real>("lights.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
