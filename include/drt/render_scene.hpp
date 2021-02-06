#ifndef DRT_RENDER_SCENE_HPP
#define DRT_RENDER_SCENE_HPP
#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/camera.hpp>
#include <drt/ray_color.hpp>
#include <drt/png.hpp>

namespace drt {

void display_progress(double progress)
{
    int barWidth = 100;

    std::cout << "\033[0;32m[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "%\033[0m\r";
    std::cout.flush();
}

template<typename Real>
void render_scene(std::string filename, int64_t image_width, int64_t image_height,
                  vec<Real> const & background, camera<Real> const & cam, hittable<Real> const & world,
                  int64_t samples_per_pixel)
{
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937_64 gen(12345);
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    for (int64_t j = 0; j < image_height; ++j) {
        Real progress = Real(j)/Real(image_height);
        display_progress(progress);
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

    drt::write_png(filename, img, image_width, image_height);
}

}
#endif
