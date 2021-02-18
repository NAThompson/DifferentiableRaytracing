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
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>

using std::make_shared;
using namespace drt;

template<typename Real>
hittable_list<Real> ellipsoid_scene() {
    hittable_list<Real> objects;
    auto light_mat = make_shared<diffuse_light<Real>>(vec<Real>(7, 7, 7));
    auto light_geom = make_shared<xz_rect<Real>>(123, 423, 147, 412, 700);
    objects.add(light_geom, light_mat);
    Real scale = 3.5;
    vec<Real> scales(scale*60, scale*55, scale*40);
    vec<Real> center(260, 250,45);
    auto boundary = make_shared<ellipsoid<Real>>(center, scales);
    auto e_mat = make_shared<dielectric<Real>>(1.5);
    objects.add(boundary, e_mat);
    auto txt = std::make_shared<isotropic<Real>>(vec<Real>(0.02, 0.4, 0.9));
    auto volume = make_shared<constant_medium<Real>>(boundary, 0.015);
    objects.add(volume, txt);
    return objects;
}

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 2880;
    const int64_t image_height = static_cast<int>(image_width / 1.6);
    const int64_t samples_per_pixel = 32;

    drt::hittable_list<Real> world = ellipsoid_scene<Real>();

    drt::vec<Real> lookfrom(478, 278, -600);
    drt::vec<Real> lookat(278, 278, 0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    render_scene<Real>("ellipsoid.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
