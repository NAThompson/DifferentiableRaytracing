#ifndef DRT_VEC_HPP
#define DRT_VEC_HPP
#include <cmath>
#include <array>
#include <iostream>
#include <initializer_list>

namespace drt {

template<typename Real, std::size_t dimension = 3>
class vec {
public:
    vec(std::initializer_list<Real> l) : v_(l)
    {
        static_assert(dimension >= 1, "Dimension >=1 is required.");
    }

    vec(Real x, Real y, Real z)
    {
        static_assert(dimension == 3, "Must have dimension = 3 to call this constructor");
        v_[0] = x;
        v_[1] = y;
        v_[2] = z;
    }

    Real operator[](size_t i) const {
        return v_[i];
    }

    vec<Real, dimension>& operator+=(const vec<Real, dimension> &w) {
        for (size_t i = 0; i < dimension; ++i) {
            v_[i] += w.v_[i];
        }
        return *this;
    }

    vec<Real, dimension>& operator*=(const Real t) {
        for (size_t i = 0; i < dimension; ++i) {
            v_[i] *= t;
        }
        return *this;
    }

    vec<Real, dimension>& operator/=(const Real t) {
        return *this *= (1/t);
    }

    friend std::ostream& operator<<(std::ostream & os, const vec<Real, dimension> & v)
    {
        os << "[";
        for (size_t i = 0; i < dimension - 1; ++i) {
            os << v[i] << ", ";
        };
        os << v[dimension - 1] << "]";
        return os;
    }

private:
    std::array<Real, dimension> v_;
};

template<typename Real, size_t dimension>
Real squared_norm(vec<Real, dimension> const & v) {
    Real norm_sq = 0;
    for (size_t i = 0; i < dimension; ++i) {
        norm_sq += v[i]*v[i];
    }
    return norm_sq;
}

template<typename Real, size_t dimension>
Real norm(vec<Real, dimension> const & v) {
    using std::sqrt;
    return sqrt(squared_norm(v));
}

template<typename Real, size_t dimension>
Real normalize(vec<Real, dimension> & v) {
    Real t = norm(v);
    for (size_t i = 0; i < dimension; ++i) {
        v[i] /= t;
    }
}


template<typename Real, size_t dimension>
Real dot(vec<Real, dimension> const & v1, vec<Real, dimension> const & v2) {
    Real d = 0;
    for (size_t i = 0; i < dimension; ++i) {
        d += v1[i]*v2[i];
    }
    return d;
}

}
#endif