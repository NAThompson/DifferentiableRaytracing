#ifndef DRT_TORUS_HPP
#define DRT_TORUS_HPP
#include <utility>
#include <algorithm>
#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/material.hpp>

namespace drt {
template<typename Real>
class torus : public hittable<Real> {
public:
    torus(vec<Real, 3> const & center, Real R, Real r, std::shared_ptr<material<Real>> mat_ptr)
       : center_(center), R_(R), r_(r), mat_ptr_(mat_ptr)
    {
        if (R <= 0) {
            std::cerr << "Big radius must be > 0\n";
        }
        if (r <= 0) {
            std::cerr << "Little radius must be > 0\n";
        }
        if (r > R/2) {
            std::cerr << "r < R/2 is required.\n";
        }
    };

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    Real gaussian_curvature(const vec<Real>& p) {
        return std::numeric_limits<Real>::quiet_NaN();
    }

    std::pair<Real, Real> gaussian_curvature_bounds() const {
        return std::make_pair<Real,Real>(std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN());
    }

    virtual ~torus() = default;

private:
    vec<Real, 3> center_;
    Real R_; // outer radius
    Real r_; // inner radius
    std::shared_ptr<material<Real>> mat_ptr_;
};

template<typename Real>
bool torus<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    vec<Real, 3> oc = r.origin() - center_;
    vec<Real, 5> c;
    // http://blog.marcinchwedczuk.pl/ray-tracing-torus
    Real nsq = oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2];
    Real tmp = (nsq - (r_*r_ + R_*R_));
    auto d = r.direction();
    Real dsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    Real ocd = dot(oc, d);
    c[0] = tmp*tmp - 4*R_*R_*(r_*r_ - oc[1]*oc[1]);
    c[1] = 4*tmp*ocd + 8*R_*R_*oc[1]*d[1];
    c[2] = 2*dsq*tmp + 4*ocd*ocd + 4*R_*R_*d[1]*d[1];
    c[3] = 4*dsq*ocd;
    c[4] = dsq*dsq;

    return false;
    Real nan = std::numeric_limits<Real>::quiet_NaN();
    rec.t = nan;
    rec.p = r(rec.t);
    vec<Real> outward_normal(nan,nan,nan);
    normalize(outward_normal);

    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr_;
    return true;
}

template<typename Real>
bool torus<Real>::bounding_box(aabb<Real>& output_box) const {
    auto b1 = center_;
    auto b2 = center_;
    b1[0] -= (R_ + r_);
    b2[0] += (R_ + r_);
    b1[1] -= (R_ + r_);
    b2[1] += (R_ + r_);
    b1[2] -= r_;
    b2[2] += r_;
    output_box = aabb<Real>(b1, b2);
    return true;
}

}
#endif