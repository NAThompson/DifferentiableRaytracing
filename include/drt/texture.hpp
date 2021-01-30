#ifndef DRT_TEXTURE_HPP
#define DRT_TEXTURE_HPP
#include <functional>
#include <memory>
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
    solid_color(vec<Real> c) : color_value_(c) {}

    solid_color(Real red, Real green, Real blue)
        : solid_color(vec<Real>(red,green,blue)) {}

    virtual vec<Real> value([[maybe_unused]] Real u, [[maybe_unused]] Real v, [[maybe_unused]] const vec<Real>& p) const override {
        return color_value_;
    }

    virtual ~solid_color() = default;

private:
    vec<Real,3> color_value_;
};

template<typename Real>
class lambda_texture : public texture<Real> {
public:
    lambda_texture(std::shared_ptr<std::function<vec<Real>(Real, Real, const vec<Real> &)>> f_ptr) : f_ptr_(f_ptr) {}

    virtual vec<Real> value([[maybe_unused]] Real u, [[maybe_unused]] Real v, [[maybe_unused]] const vec<Real>& p) const override {
        return f_ptr_->operator()(u,v, p);
    }

    virtual ~lambda_texture() = default;

private:
    std::shared_ptr<std::function<vec<Real>(Real, Real, const vec<Real> &)>> f_ptr_;
};


}
#endif
