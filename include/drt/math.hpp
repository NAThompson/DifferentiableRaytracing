#ifndef DRT_MATH_HPP
#define DRT_MATH_HPP
#include <cmath>

namespace drt {

// Computes ab - cd.
// See: https://pharr.org/matt/blog/2019/11/03/difference-of-floats.html
template<typename Real>
inline Real difference_of_products(Real a, Real b, Real c, Real d) {
    Real cd = c * d;
    Real err = std::fma(-c, d, cd);
    Real dop = std::fma(a, b, -cd);
    return dop + err;
}

}
#endif
