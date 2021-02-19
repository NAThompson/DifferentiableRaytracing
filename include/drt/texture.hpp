#ifndef DRT_TEXTURE_HPP
#define DRT_TEXTURE_HPP
#include <functional>
#include <memory>
#include <drt/vec.hpp>
#include <drt/hittable.hpp>

namespace drt {

template<typename Real>
class texture {
public:
    virtual vec<Real, 3> value(hit_record<Real> const & hr) const = 0;
};

template<typename Real>
class solid_color : public texture<Real> {
public:
    solid_color() {}
    solid_color(vec<Real> c) : color_value_(c) {}

    solid_color(Real red, Real green, Real blue)
        : solid_color(vec<Real>(red,green,blue)) {}

    virtual vec<Real> value([[maybe_unused]] hit_record<Real> const & hr) const override {
        return color_value_;
    }

    virtual ~solid_color() = default;

private:
    vec<Real,3> color_value_;
};

template<typename Real>
class lambda_texture : public texture<Real> {
public:
    lambda_texture(std::function<vec<Real>(const hit_record<Real> &)> f) : f_(f) {}

    virtual vec<Real> value(const hit_record<Real> & hr) const override {
        return f_(hr);
    }

    virtual ~lambda_texture() = default;

private:
    std::function<vec<Real>(const hit_record<Real> &)> f_;
};


}
#endif
