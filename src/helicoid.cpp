#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>
#include <drt/lambertian.hpp>
#include <drt/texture.hpp>
#include <drt/diffuse_light.hpp>
#include <drt/aarect.hpp>
#include <drt/color_maps.hpp>
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>
#include <drt/helicoid.hpp>
#include <drt/dielectric.hpp>
#include <drt/cylinder.hpp>
#include <drt/disk.hpp>

using std::make_shared;
using std::log;
using namespace drt;

template<typename Real>
hittable_list<Real> helicoid_scene() {
    hittable_list<Real> objects;
    Real radius = 1;
    Real speed = 3;
    auto heli = make_shared<helicoid<Real>>(radius, speed);
    auto bounds = heli->gaussian_curvature_bounds();
    Real ln_nkmin = log(-bounds.second);
    Real ln_nkmax = log(-bounds.first);
    auto gaussian_curvature = [=](hit_record<Real> const & hr) {
        Real kappa = hr.gaussian_curvature();
        Real ln_nk = log(-kappa);
        Real scalar = (ln_nk - ln_nkmin)/(ln_nkmax - ln_nkmin);
        return smooth_cool_warm(scalar);
    };

    auto texture = make_shared<lambda_texture<Real>>(gaussian_curvature);
    auto mat = make_shared<lambertian<Real>>(texture);
    objects.add(heli, mat);

    // Adding this bounding cylinder is a really nice effect with a lot of samples per pixel.
    // Not so much without.
    /*auto cyl = make_shared<cylinder<Real>>(radius, -speed/2, speed/2);
    auto emat =  make_shared<dielectric<Real>>(1.0);
    objects.add(cyl, emat);*/
    return objects;
}


int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 1200;
    const int64_t image_height = static_cast<int64_t>(image_width/aspect_ratio);
    const int64_t samples_per_pixel = 1024;

    auto world = helicoid_scene<Real>();
    drt::vec<Real> lookat(0,0,0);
    drt::vec<Real> lookfrom(0, -4, 0);
    drt::vec<Real> vup(0,1,1);
    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);

    auto [o, d] = cam.backlight();
    auto disk_ptr = make_shared<disk<Real>>(10.0, o, d);
    auto mat = make_shared<diffuse_light<Real>>(vec<Real>(2,2,2));
    world.add(disk_ptr, mat);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    drt::render_scene<Real>("helicoid.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
