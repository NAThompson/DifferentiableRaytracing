#ifndef DRT_RENDER_SCENE_HPP
#define DRT_RENDER_SCENE_HPP
#include <iostream>
#include <chrono>
#include <drt/color_maps.hpp>
#include <drt/vec.hpp>
#include <drt/camera.hpp>
#include <drt/ray_color.hpp>
#include <drt/png.hpp>
#include <drt/pcg.hpp>
#include <drt/progress_bar.hpp>
#include <omp.h>

namespace drt {

#define DRT_COST_MAP 1


template<typename Real>
void render_scene(std::string filename, int64_t image_width, int64_t image_height,
                  vec<Real> const & background, camera<Real> const & cam, hittable_list<Real> const & world,
                  int64_t samples_per_pixel)
{
    pcg32_random_t pcg;
    pcg.state = 1235123;
    pcg.inc = 234121;
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    std::vector<double> cost(image_width*image_height, 0);
    #pragma omp parallel for
    for (int64_t j = 0; j < image_height; ++j) {
        //Real progress = Real(j)/Real(image_height);
        //display_progress(progress);
        std::cout << j << ", " << std::flush;
        for (int64_t i = 0; i < image_width; ++i) {
            #if DRT_COST_MAP
            auto start = std::chrono::high_resolution_clock::now();
            #endif
            drt::vec<Real, 3> color(0,0,0);
            for (int64_t s = 0; s < samples_per_pixel; ++s) {
                Real u = (Real(i) + pcg_real_01<Real>(&pcg)) / (image_width -1);
                Real v = (Real(j) + pcg_real_01<Real>(&pcg)) / (image_height -1);
                auto r = cam.get_ray(u, v);
                color += ray_color(r, background, world, max_depth);
            }
            auto c = drt::to_8bit_rgba(color/Real(samples_per_pixel));
            int64_t idx = 4 * image_width * (image_height - 1 - j) + 4 * i;
            img[idx + 0] = c[0];
            img[idx + 1] = c[1];
            img[idx + 2] = c[2];
            img[idx + 3] = c[3];
            #if DRT_COST_MAP
            auto end = std::chrono::high_resolution_clock::now();
            auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            cost[image_width * (image_height - 1 - j) + i] = time_span.count();
            #endif
        }
    }
    display_progress(1.0);

    drt::write_png(filename, img, image_width, image_height);
    #if DRT_COST_MAP
    double min_time = std::numeric_limits<double>::max();
    double max_time = 0;
    for (auto & t : cost) {
        if (t > max_time) {
            max_time = t;
        }
        if (t < min_time) {
            min_time = t;
        }
    }
    if (min_time < 0) {
        min_time = 0;
    }

    for (int64_t j = 0; j < image_height; ++j) {
        for (int64_t i = 0; i < image_width; ++i) {
            int64_t idx = image_width * (image_height - 1 - j) + i;
            //Real s = (log(cost[idx]) - log(min_time))/(log(max_time) - log(min_time));
            Real s = (cost[idx] - min_time)/(max_time - min_time);
            auto color = drt::inferno(s);
            auto c = drt::to_8bit_rgba(color);
            idx *= 4;
            img[idx + 0] = c[0];
            img[idx + 1] = c[1];
            img[idx + 2] = c[2];
            img[idx + 3] = c[3];
        }
    }
    std::string cost_filename = filename.substr(0, filename.size() - 4);
    cost_filename += "_cost.png";
    drt::write_png(cost_filename, img, image_width, image_height);
    #endif
}

}
#endif
