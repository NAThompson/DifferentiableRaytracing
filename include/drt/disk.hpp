#ifndef DRT_DISK_HPP
#define DRT_DISK_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/matrix.hpp>
#include <drt/roots.hpp>

namespace drt {

template<typename Real>
class disk : public hittable<Real> {
public:

    disk(Real radius, vec<Real,3> center, vec<Real,3> normal)
       : radius_(radius), center_(center), normal_(normal)
    {
        if (radius_ <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Radius must be > 0 for a disk.\n";
        }
        Real length = squared_norm(normal);
        if (abs(length - 1) > 10*std::numeric_limits<Real>::epsilon()) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Normal must be normalized.\n";
        }

        // See: https://math.stackexchange.com/questions/61547/rotation-of-a-vector-distribution-to-align-with-a-normal-vector
        Real nz = normal[2];
        if (nz == -1) {
            for (int64_t i = 0; i < 3; ++i) {
                for (int64_t j = 0; j < 3; ++j) {
                    to_world_(i,j) = 0;
                    to_local_(i,j) = 0;
                }
            }
            to_world_(0,0) = -1;
            to_world_(1,1) = 1;
            to_world_(2,2) = -1;

            to_local_(0,0) = -1;
            to_local_(1,1) = 1;
            to_local_(2,2) = -1;
            //std::cerr << "Rotating into nz = -1!\n";
            return;
        }
        Real nx = normal[0];
        Real ny = normal[1];
        Real t = 1/(1+nz);

        to_world_(0,0) = nz + ny*ny*t;
        to_world_(0,1) = -nx*ny*t;
        to_world_(0,2) = nx;
        to_world_(1,0) = to_world_(0,1);
        to_world_(2,0) = -to_world_(0,2);
        to_world_(1,1) = nz + t*nx*nx;
        to_world_(2,2) = nz;
        to_world_(1,2) = ny;
        to_world_(2,1) = -ny;

        Real det = determinant(to_world_);
        if (abs(det - 1) > 250*std::numeric_limits<Real>::epsilon()) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Determinant of rotation matrix is "
                      << det << " = " << std::hexfloat << det << std::defaultfloat << "\n";
            std::cerr << to_world_ << "\n";
        }
        to_local_ = inverse(to_world_);
        det = determinant(to_local_);
        if (abs(det - 1) > 250*std::numeric_limits<Real>::epsilon()) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Determinant of rotation matrix is "
                      << det << " = " << std::hexfloat << det << std::defaultfloat << "\n";
        }
    };

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~disk() = default;


    // σ(u,v) = R(rvcos(2πu), rvsin(2πu), 0) + C
    vec<Real,3> operator()(Real u, Real v) const override {
        vec<Real> w;
        w[0] = radius_*v*std::cos(2*M_PI*u);
        w[1] = radius_*v*std::sin(2*M_PI*u);
        w[2] = 0;
        return to_world_*w + center_;
    }

private:
    Real radius_;
    vec<Real,3> center_;
    vec<Real,3> normal_;
    matrix<Real,3,3> to_world_;
    matrix<Real,3,3> to_local_;
};

template<typename Real>
bool disk<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    using std::sqrt;
    using std::atan2;
    vec<Real, 3> o = to_local_*(r.origin() - center_);
    auto d = to_local_*r.direction();
    if (d[2] == 0) {
        return false;
    }

    Real t = -o[2]/d[2];
    if (t < t_min || t > t_max) {
        return false;
    }

    Real x = o[0] + t*d[0];
    Real y = o[1] + t*d[1];
    Real ssq = x*x + y*y;
    if (ssq > radius_*radius_) {
        return false;
    }
    
    rec.t = t;
    rec.p = r(rec.t);
    rec.v = sqrt(ssq)/radius_;
    if (rec.v == 0) { 
        rec.u = 0;
    }
    else {
        rec.u = atan2(y,x)/(2*M_PI);
        if (rec.u < 0) {
            rec.u += 1;
        }
    }
    rec.set_face_normal(r, normal_);
    rec.E = 4*M_PI*M_PI*radius_*radius_*rec.v*rec.v;
    rec.F = 0;
    rec.G = radius_*radius_;
    rec.e = std::numeric_limits<Real>::quiet_NaN();
    rec.f = std::numeric_limits<Real>::quiet_NaN();
    rec.g = 0;

    return true;
}

template<typename Real>
bool disk<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(
        vec<Real>(-radius_, -radius_, center_[2]),
        vec<Real>(radius_, radius_, center_[2]));
    return true;
}

}
#endif
