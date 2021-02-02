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

using std::make_shared;
using drt::lambertian;
using drt::vec;
using drt::sphere;
using drt::dielectric;
using drt::constant_medium;
using drt::metal;
using drt::ellipsoid;
using drt::diffuse_light;
using drt::xz_rect;
using drt::hittable_list;

template<typename Real>
hittable_list<Real> ellipsoid_scene() {
    hittable_list<Real> objects;
    auto light = make_shared<diffuse_light<Real>>(vec<Real>(7, 7, 7));
    objects.add(make_shared<xz_rect<Real>>(123, 423, 147, 412, 700, light));
    Real scale = 3.5;
    vec<Real> scales(scale*60, scale*55, scale*40);
    vec<Real> center(260, 250,45);
    auto boundary = make_shared<ellipsoid<Real>>(center, scales, make_shared<dielectric<Real>>(1.5));
    objects.add(boundary);
    objects.add(make_shared<constant_medium<Real>>(boundary, 0.015, vec<Real>(0.02, 0.4, 0.9)));
    return objects;
}

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

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 2880;
    const int64_t image_height = static_cast<int>(image_width / 1.6);
    const int64_t samples_per_pixel = 32;

    drt::hittable_list<Real> world1 = ellipsoid_scene<Real>();
    drt::bvh_node<Real> world(world1);

    drt::vec<Real> lookfrom(478, 278, -600);
    drt::vec<Real> lookat(278, 278, 0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937_64 gen;
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    for (int64_t j = 0; j < image_height; ++j) {
        std::cerr << j << "/" << image_height << "\r";
        for (int64_t i = 0; i < image_width; ++i) {
            drt::vec<Real, 3> color(0,0,0);
            for (int64_t s = 0; s < samples_per_pixel; ++s) {
                Real u = (Real(i) + dis(gen)) / (image_width -1);
                Real v = (Real(j) + dis(gen)) / (image_height -1);
                auto r = cam.get_ray(u, v);
                color += ray_color(r, background, world, max_depth);
            }
            auto c = drt::to_8bit_rgba(color/Real(samples_per_pixel));
            int64_t idx = 4 * image_width * (image_height - 1 - j) + 4 * i;
            img[idx + 0] = c[0];
            img[idx + 1] = c[1];
            img[idx + 2] = c[2];
            img[idx + 3] = c[3];
        }
    }

    drt::write_png("ellipsoid.png", img, image_width, image_height);
}
