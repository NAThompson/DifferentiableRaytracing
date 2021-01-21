#ifndef COLOR_TABLES_HPP
#define COLOR_TABLES_HPP
#include <drt/vec.hpp>

namespace drt {

// Assumption: the normal is normalized!
template<typename Real>
drt::vec<Real, 3> color_by_normal(drt::vec<Real, 3> const & normal, Real attenutation = 0.5) {
    return attenuation*drt::vec<Real, 3>(normal[0]+1, normal[1]+1, normal[2]+1);
}


}
#endif