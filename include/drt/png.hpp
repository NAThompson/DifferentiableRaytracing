#ifndef DRT_PNG_HPP
#define DRT_PNG_HPP
#include <algorithm>
#include "lodepng.h"

namespace drt {

template<typename Real, size_t dimension>
std::array<uint8_t, 4> to_8bit_rgba(vec<Real, dimension> const & v)
{
    static_assert(dimension == 3 || dimension == 4, "The color must be RGB[0,1] or RGBA[0,1]");
    std::array<uint8_t, 4> pixel;
    for (size_t i = 0; i < dimension; ++i) {
        // Clamp to [0, 1] or assert when v[i] \in \mathbb{R} \setminus [0,1]? OMG questions questions questions.
        pixel[i] = 255*std::clamp(v[i], Real(0), Real(1));
    }

    if constexpr (dimension == 3) {
        pixel[3] = 255;
    }
    return pixel;
}

// In lodepng, the vector is expected to be row major, with the top row specified first.
// Note that this is a bit confusing sometimes as it's more natural to let y increase moving *up*.
unsigned write_png(std::string const & filename, std::vector<uint8_t> const & img, size_t width, size_t height)
{
    unsigned error = lodepng::encode(filename, img, width, height, LodePNGColorType::LCT_RGBA, 8);
    if(error)
    {
        std::cerr << "Error encoding png: " << lodepng_error_text(error) << "\n";
    }
    return error;
}

}
#endif