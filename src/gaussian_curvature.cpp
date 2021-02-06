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
#include <drt/ray_color.hpp>

using std::make_shared;
using drt::lambertian;
using drt::vec;
using drt::sphere;
using drt::dielectric;
using drt::constant_medium;
using drt::metal;
using drt::ellipsoid;
using drt::diffuse_light;
using drt::xy_rect;
using drt::viridis;
using drt::lambda_texture;

template<typename Real>
drt::hittable_list<Real> ellipsoid_scene() {

    drt::hittable_list<Real> objects;
    auto light = make_shared<diffuse_light<Real>>(vec<Real>(1, 1, 1));
    objects.add(make_shared<xy_rect<Real>>(-30, 30, -30, 30, -15, light));

    Real scale = 1.5;
    Real a = Real(scale*1.3);
    Real b = Real(scale*1);
    Real c = Real(scale*1.6);
    vec<Real> scales(a,b,c);
    vec<Real> center = vec<Real>(0, 0, 0);
    Real asq = a*a;
    Real bsq = b*b;
    Real csq = c*c;
    Real kappa_min = std::min({csq/(asq*bsq), asq/(csq*bsq), bsq/(asq*csq)});
    Real kappa_max = std::max({csq/(asq*bsq), asq/(csq*bsq), bsq/(asq*csq)});

    std::function<vec<Real>(Real, Real, const vec<Real> &)> gaussian_curvature = [=]([[maybe_unused]] Real u, [[maybe_unused]] Real v, vec<Real> const & p) {
        // https://mathworld.wolfram.com/Ellipsoid.html, equation 14:
        Real numerator = asq*bsq*bsq*bsq*csq*csq*csq;
        Real y = p[1] - center[1];
        Real z = p[2] - center[2];
        Real sqrt_denom = csq*csq*bsq*bsq + csq*csq*(asq - bsq)*y*y + bsq*bsq*(asq-csq)*z*z;
        Real kappa = numerator/(sqrt_denom*sqrt_denom);

        Real scalar = (kappa - kappa_min)/(kappa_max - kappa_min);
        vec<Real, 3> w = viridis(scalar);
        return w;
    };

    auto ptr = make_shared<decltype(gaussian_curvature)>(gaussian_curvature);
    auto ltext = drt::lambda_texture(ptr);
    auto ltext_ptr = make_shared<decltype(ltext)>(ltext);

    auto mat = make_shared<lambertian<Real>>(ltext_ptr);
    auto boundary = make_shared<ellipsoid<Real>>(center, scales, mat);
    objects.add(boundary);
    return objects;
}

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 800;
    const int64_t image_height = static_cast<int>(image_width / 1.6);
    const int64_t samples_per_pixel = 64;

    drt::hittable_list<Real> world1 = ellipsoid_scene<Real>();
    drt::bvh_node<Real> world(world1);

    drt::vec<Real> lookfrom(5, 0, -5);
    drt::vec<Real> lookat(0, 0, 0);
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

    drt::write_png("gaussian_curvature.png", img, image_width, image_height);
}
