#include <drt/vec.hpp>
#include <drt/png.hpp>


int main() {
    using Real = double;

    drt::vec<Real, 3> v(1.2, 3.4, 7.5);

    std::cout << v << "\n";

    const int64_t width = 2880;
    const int64_t height = width/1.618;
    std::vector<uint8_t> img(4*width*height, 0);
    for (int64_t j = 0; j < height; ++j) {
        for (int64_t i = 0; i < width; ++i) {
            uint8_t r = 256*Real(i) / (width);
            uint8_t g = 256*Real(j) / (height);
            uint8_t b = 256/4;

            img[4 * width * j + 4 * i + 0] = r;
            img[4 * width * j + 4 * i + 1] = g;
            img[4 * width * j + 4 * i + 2] = b;
            img[4 * width * j + 4 * i + 3] = 255;
        }
    }

    drt::write_png("first.png", img, width, height);
}