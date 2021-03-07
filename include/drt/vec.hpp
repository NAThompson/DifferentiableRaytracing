#ifndef DRT_VEC_HPP
#define DRT_VEC_HPP
#include <cmath>
#include <iostream>
#include <initializer_list>

namespace drt {

template<typename Real, int64_t dimension = 3>
class vec {
public:
    vec() {};

    vec(Real x, Real y, Real z)
    {
        static_assert(dimension == 3, "Must have dimension = 3 to call this constructor");
        v_[0] = x;
        v_[1] = y;
        v_[2] = z;
    }

    vec(Real x, Real y)
    {
        static_assert(dimension == 2, "Must have dimension = 2 to call this constructor");
        v_[0] = x;
        v_[1] = y;
    }

    Real operator[](int64_t i) const {
        return v_[i];
    }

    Real& operator[](int64_t i) {
        return v_[i];
    }

    vec<Real, dimension> operator-() const {
        vec<Real, dimension> w = *this;
        for (int64_t i = 0; i < dimension; ++i) {
            w[i] = -w[i];
        }
        return w;
    }

    vec<Real, dimension>& operator+=(const vec<Real, dimension> &w) {
        for (int64_t i = 0; i < dimension; ++i) {
            v_[i] += w.v_[i];
        }
        return *this;
    }

    vec<Real, dimension>& operator*=(const Real t) {
        for (int64_t i = 0; i < dimension; ++i) {
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
        for (int64_t i = 0; i < dimension - 1; ++i) {
            os << v[i] << ", ";
        };
        os << v[dimension - 1] << "]";
        return os;
    }

private:
    Real v_[dimension];
};

template<typename Real, int64_t dimension>
Real squared_norm(vec<Real, dimension> const & v) {
    Real norm_sq = 0;
    for (int64_t i = 0; i < dimension; ++i) {
        norm_sq += v[i]*v[i];
    }
    return norm_sq;
}

template<typename Real, int64_t dimension>
Real norm(vec<Real, dimension> const & v) {
    using std::sqrt;
    return sqrt(squared_norm(v));
}

template<typename Real, int64_t dimension>
Real max_norm(vec<Real, dimension> const & v) {
    using std::max;
    using std::abs;
    Real m = abs(v[0]);
    for (int64_t i = 1; i < dimension; ++i) {
        if (abs(v[i]) > m) {
            m = abs(v[i]);
        }
    }
    return m;
}

template<typename Real, int64_t dimension>
void normalize(vec<Real, dimension> & v) {
    Real t = norm(v);
    for (int64_t i = 0; i < dimension; ++i) {
        v[i] /= t;
    }
}

template<typename Real, int64_t dimension>
Real dot(vec<Real, dimension> const & v1, vec<Real, dimension> const & v2) {
    Real d = 0;
    for (int64_t i = 0; i < dimension; ++i) {
        d += v1[i]*v2[i];
    }
    return d;
}

template<typename Real, int64_t dimension>
inline vec<Real, dimension> operator+(const vec<Real, dimension> &u, const vec<Real, dimension> &v) {
    vec<Real, dimension> w;
    for (int64_t i = 0; i < dimension; ++i) {
        w[i] = u[i] + v[i];
    }
    return w;
}

template<typename Real, int64_t dimension>
inline vec<Real, dimension> operator-(const vec<Real, dimension> &u, const vec<Real, dimension> &v) {
    vec<Real, dimension> w;
    for (int64_t i = 0; i < dimension; ++i) {
        w[i] = u[i] - v[i];
    }
    return w;

}

template<typename Real, int64_t dimension>
inline vec<Real, dimension> operator/(vec<Real, dimension> const & v, Real t) {
    vec<Real, dimension> w;
    for (int64_t i = 0; i < dimension; ++i) {
        w[i] = v[i]/t;
    }
    return w;
}

template<typename Real, int64_t dimension>
inline vec<Real, dimension> operator*(Real t, const vec<Real, dimension> &v) {
    vec<Real, dimension> w;
    for (int64_t i = 0; i < dimension; ++i) {
        w[i] = t*v[i];
    }
    return w;
}

template<typename Real, int64_t dimension>
inline vec<Real, dimension> operator*(const vec<Real, dimension> &v, Real t) {
    return t * v;
}

template<typename Real, int64_t dimension>
vec<Real, dimension> reflect(const vec<Real, dimension>& v, const vec<Real, dimension>& n)
{
    return v - 2*dot(v,n)*n;
}

template<typename Real>
vec<Real> cross(const vec<Real>& u, const vec<Real>& v) {
    vec<Real> w;
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
    return w;
}

}
#endif
