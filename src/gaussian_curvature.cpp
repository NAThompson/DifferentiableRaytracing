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
#include <drt/render_scene.hpp>

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
    auto light_mat = make_shared<diffuse_light<Real>>(vec<Real>(1, 1, 1));
    auto light_geom = make_shared<xy_rect<Real>>(-30, 30, -30, 30, -15);
    objects.add(light_geom, light_mat);

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

    std::function<vec<Real>(const drt::hit_record<Real> &)> gaussian_curvature = [=](drt::hit_record<Real> const & hr) {
        Real kappa = hr.gaussian_curvature();

        Real scalar = (kappa - kappa_min)/(kappa_max - kappa_min);
        vec<Real, 3> w = viridis(scalar);
        return w;
    };

    auto ptr = make_shared<decltype(gaussian_curvature)>(gaussian_curvature);
    auto ltext = drt::lambda_texture(ptr);
    auto ltext_ptr = make_shared<decltype(ltext)>(ltext);

    auto mat = make_shared<lambertian<Real>>(ltext_ptr);
    auto boundary = make_shared<ellipsoid<Real>>(center, scales);
    objects.add(boundary, mat);
    return objects;
}

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 800;
    const int64_t image_height = static_cast<int>(image_width / 1.6);
    const int64_t samples_per_pixel = 64;

    drt::hittable_list<Real> world = ellipsoid_scene<Real>();

    drt::vec<Real> lookfrom(5, 0, -5);
    drt::vec<Real> lookat(0, 0, 0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    drt::render_scene<Real>("gaussian_curvature.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
