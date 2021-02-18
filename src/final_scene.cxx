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
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>

template<typename Real>
drt::hittable_list<Real> final_scene() {
    using std::make_shared;
    using drt::lambertian;
    using drt::vec;
    using drt::sphere;
    using drt::dielectric;
    using drt::constant_medium;
    using drt::metal;

    drt::hittable_list<Real> boxes1;
    auto ground = make_shared<lambertian<Real>>(vec<Real>(0.48, 0.83, 0.53));
    std::uniform_real_distribution<Real> dis(1,101);
    std::random_device rd;
    const int boxes_per_side = 20;
    for (int i = 0; i < boxes_per_side; i++) {
        for (int j = 0; j < boxes_per_side; j++) {
            Real w = 100.0;
            Real x0 = -1000.0 + i*w;
            Real z0 = -1000.0 + j*w;
            Real y0 = 0.0;
            Real x1 = x0 + w;
            Real y1 = dis(rd);
            Real z1 = z0 + w;
            auto box = make_shared<drt::box<Real>>(vec<Real>(x0,y0,z0), vec<Real>(x1,y1,z1));
            boxes1.add(box, ground);
        }
    }

    drt::hittable_list<Real> objects;
    objects.add(make_shared<drt::bvh_node<Real>>(boxes1));

    auto light = make_shared<drt::diffuse_light<Real>>(vec<Real>(7, 7, 7));
    objects.add(make_shared<drt::xz_rect<Real>>(123, 423, 147, 412, 554), light);
    auto di15 = make_shared<dielectric<Real>>(1.5);
    objects.add(make_shared<sphere<Real>>(vec<Real>(260, 150, 45), Real(50)), di15);
    auto m1 = std::make_shared<metal<Real>>(vec<Real>(0.8, 0.8, 0.9));
    auto s1 = std::make_shared<sphere<Real>>(vec<Real>(0, 150, 145), Real(50));
    objects.add(s1, m1);

    auto boundary = make_shared<sphere<Real>>(vec<Real>(360,150,145), Real(70));
    objects.add(boundary, di15);
    objects.add(make_shared<constant_medium<Real>>(boundary, 0.2, vec<Real>(0.2, 0.4, 0.9)));
    boundary = make_shared<sphere<Real>>(vec<Real>(0, 0, 0), 5000);
    auto txt = solid_color(vec<Real>(1,1,1));
    objects.add(make_shared<constant_medium<Real>>(boundary, .0001), std::make_shared(txt));

    return objects;
}

int main() {
    using Real = float;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 1800;
    const int64_t image_height = static_cast<int>(image_width / aspect_ratio);
    const int64_t samples_per_pixel = 256;

    drt::hittable_list<Real> world = final_scene<Real>();
    //drt::hittable_list<Real> world1 = final_scene<Real>();
    //drt::bvh_node<Real> world(world1);

    drt::vec<Real> lookfrom(478, 278, -600);
    drt::vec<Real> lookat(278, 278, 0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    drt::render_scene<Real>("final_scene.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
