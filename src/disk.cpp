#include <iostream>
#include <random>
#include <cmath>
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
#include <drt/disk.hpp>
#include <drt/dielectric.hpp>

using std::make_shared;
using std::log;
using namespace drt;

template<typename Real>
hittable_list<Real> disk_scene() {
    hittable_list<Real> objects;
    Real radius = 1;
    vec<Real,3> center(0,0,0);
    Real inv_rt2 = 1.0/sqrt(2.0);
    vec<Real,3> normal(inv_rt2,0,-inv_rt2);
    auto disk_ptr = make_shared<disk<Real>>(radius, center, normal);

    auto coord = [=](hit_record<Real> const & hr) {
        return viridis(hr.u);
    };

    auto texture = make_shared<lambda_texture<Real>>(coord);
    auto mat = make_shared<lambertian<Real>>(texture);
    objects.add(disk_ptr, mat);
    return objects;
}


int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 800;
    const int64_t image_height = static_cast<int64_t>(image_width/aspect_ratio);
    const int64_t samples_per_pixel = 256;

    auto world = disk_scene<Real>();
    drt::vec<Real> lookat(0,0,0);
    drt::vec<Real> lookfrom(0, -4, 2);
    drt::vec<Real> vup(0,0,1);
    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    auto [o, d] = cam.backlight();
    auto disk_ptr = make_shared<disk<Real>>(10.0, o, d);
    auto mat = make_shared<diffuse_light<Real>>(vec<Real>(5,5,5));
    world.add(disk_ptr, mat);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    drt::render_scene<Real>("disk.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
