#ifndef DRT_TEXTURE_H
#define DRT_TEXTURE_H
#include <drt/vec.hpp>

namespace drt {

template<typename Real>
class texture {
public:
    virtual vec<Real, 3> value(Real u, Real v, const vec<Real, 3>& p) const = 0;
};

template<typename Real>
class solid_color : public texture<Real> {
public:
    solid_color() {}
    solid_color(vec<Real> c) : color_value(c) {}

    solid_color(Real red, Real green, Real blue)
        : solid_color(vec<Real>(red,green,blue)) {}

    virtual vec<Real> value([[maybe_unused]] Real u, [[maybe_unused]] Real v, [[maybe_unused]] const vec<Real>& p) const override {
        return color_value;
    }

    virtual ~solid_color() = default;

private:
    vec<Real,3> color_value;
};

}
#endif
