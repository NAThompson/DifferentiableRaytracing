#ifndef DRT_RAY_HPP
#define DRT_RAY_HPP
#include <drt/vec.hpp>
namespace drt {

template<typename Real>
class ray {
public:
    ray(vec<Real, 3> const & origin, vec<Real,3> const & direction) : origin_{origin}, direction_{direction}
    {}

    ray()
    {}

    vec<Real, 3> origin() const {
        return origin_;
    }

    vec<Real, 3> direction() const {
        return direction_;
    }

    vec<Real, 3> operator()(Real t) const {
        return origin_ + t*direction_;
    }

private:
    vec<Real, 3> origin_;
    vec<Real, 3> direction_;
};

}
#endif