#ifndef CHECKER_TEXTURE_HPP
#define CHECKER_TEXTURE_HPP
#include <cmath>
#include <drt/vec.hpp>
#include <drt/texture.hpp>

namespace drt {

template<typename Real>
class checker_texture : public texture<Real> {
public:
    checker_texture() {}

    checker_texture(std::shared_ptr<texture<Real>> even, std::shared_ptr<texture<Real>> odd)
        : even_(even), odd_(odd) {}

    checker_texture(vec<Real> c1, vec<Real> c2)
        : even_(std::make_shared<solid_color<Real>>(c1)) , odd_(std::make_shared<solid_color<Real>>(c2)) {}

    virtual vec<Real> value(const hit_record<Real>& hr) const override {
        auto const & p = hr.p;
        auto sines = std::sin(10*p[0])*std::sin(10*p[1])*std::sin(10*p[2]);
        if (sines < 0)
            return odd_->value(hr);
        else
            return even_->value(hr);
    }

    virtual ~checker_texture() = default;

public:
    std::shared_ptr<texture<Real>> even_;
    std::shared_ptr<texture<Real>> odd_;
};

}
#endif