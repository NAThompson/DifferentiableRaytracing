#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/png.hpp>
#include <drt/ray.hpp>
#include <drt/aarect.hpp>
#include <drt/hittable.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>
#include <drt/material.hpp>
#include <drt/lambertian.hpp>
#include <drt/texture.hpp>
#include <drt/diffuse_light.hpp>
#include <drt/ellipsoid.hpp>
#include <drt/color_maps.hpp>
#include <drt/render_scene.hpp>

using std::make_shared;
using namespace drt;

template<typename Real>
hittable_list<Real> ellipsoid_scene() {
    hittable_list<Real> objects;
    auto light_mat = make_shared<diffuse_light<Real>>(vec<Real>(1, 1, 1));
    auto light_geom = make_shared<xy_rect<Real>>(-30, 30, -30, 30, -15);
    objects.add(light_geom, light_mat);

    Real scale = 1.5;
    Real a = Real(scale*1.3);
    Real b = Real(scale*1);
    Real c = Real(scale*1.6);
    vec<Real> scales(a,b,c);
    vec<Real> center = vec<Real>(0, 0, 0);
    auto ell = make_shared<ellipsoid<Real>>(center, scales);
    auto bounds = ell->gaussian_curvature_bounds();
    Real kappa_min = bounds.first;
    Real kappa_max = bounds.second;

    auto gaussian_curvature = [=](hit_record<Real> const & hr) {
        Real kappa = hr.gaussian_curvature();
        Real scalar = (kappa - kappa_min)/(kappa_max - kappa_min);
        return viridis(scalar);
    };

    auto texture = make_shared<lambda_texture<Real>>(gaussian_curvature);
    auto mat = make_shared<lambertian<Real>>(texture);
    objects.add(ell, mat);
    return objects;
}

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 800;
    const int64_t image_height = static_cast<int>(image_width / 1.6);
    const int64_t samples_per_pixel = 64;

    hittable_list<Real> world = ellipsoid_scene<Real>();

    vec<Real> lookfrom(5, 0, -5);
    vec<Real> lookat(0, 0, 0);
    vec<Real> vup(0,1,0);

    camera cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    vec<Real> background(0.0, 0.0, 0.0);
    render_scene<Real>("gaussian_curvature.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
