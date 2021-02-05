#ifndef DRT_HITTABLE_HPP
#define DRT_HITTABLE_HPP
#include <limits>
#include <drt/ray.hpp>
#include <drt/vec.hpp>
#include <drt/aabb.hpp>


namespace drt {

template<typename Real>
class material;

template<typename Real>
struct hit_record {
    vec<Real, 3> p;
    vec<Real, 3> normal;
    Real t;
    // Parametric coordinates on the surface.
    // Not used by everything:
    Real u = std::numeric_limits<Real>::quiet_NaN();
    Real v = std::numeric_limits<Real>::quiet_NaN();
    bool front_face;
    std::shared_ptr<material<Real>> mat_ptr;

    inline void set_face_normal(const ray<Real>& r, const vec<Real, 3>& outward_normal) {
        front_face = drt::dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }

    friend std::ostream& operator<<(std::ostream & os, const hit_record<Real> & rec) {
        os << "Ray intersects object at " << rec.p << " and time " << rec.t << "\n";
        os << "Normal is " << rec.normal << "\n";
        os << "Parametric coordinates: (u,v) = (" << rec.u << ", " << rec.v << ")\n";
        if (rec.front_face) {
            os << "Ray intersects front face.\n";
        } else {
            os << "Ray intersects back face.\n";
        }
        return os;
    }
};

template<typename Real>
class hittable {
public:
    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const = 0;

    virtual bool bounding_box(aabb<Real>& output_box) const = 0;
};

}
#endif