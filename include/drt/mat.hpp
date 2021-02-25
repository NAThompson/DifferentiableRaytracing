#ifndef DRT_MAT_HPP
#define DRT_MAT_HPP
#include <cmath>
#include <array>
#include <iostream>
#include <initializer_list>

namespace drt {

template<typename Real, int64_t rows = 3, int64_t cols = 3>
class mat {
public:

    mat() {};

    Real operator()(int64_t i, int64_t j) const {
        return M_[i*cols + j];
    }

    Real& operator()(int64_t i, int64_t j) {
        return M_[i*cols + j];
    }

    // Return an x that satisfies Ax = b.
    vec<Real, cols> solve(vec<Real, rows> const & b) const;

    friend std::ostream& operator<<(std::ostream & os, const mat<Real, rows, cols> & v)
    {
        for (int64_t i = 0; i < rows; ++i) {
            os << "[";
            for (int64_t j = 0; j < cols - 1; ++j) {
                os << v.M_[i*cols + j] << ", ";
            }
            os << v.M_[(i+1)*cols - 1] << "]\n";
        }
        return os;
    }

private:
    Real M_[rows*cols];
};

template<typename Real, int64_t rows, int64_t cols>
Real determinant(mat<Real, rows, cols> const & m)
{
    static_assert(rows == cols, "Must have the same number of rows as columns to compute a determinant.");
    static_assert(rows < 4, "Determinant not implemented for n > 3.");
    if constexpr (rows == 2) {
        return m(0,0)*m(1,1) - m(0,1)*m(1,0);
    }
    if constexpr (rows == 3) {
        return m(0,0)*(m(1,1)*m(2,2) - m(1,2)*m(2,1)) - m(0,1)*(m(1,0)*m(2,2) - m(1,2)*m(2,0)) + m(0,2)*(m(1,0)*m(2,1) - m(1,1)*m(2,0));

    }
    return std::numeric_limits<Real>::quiet_NaN();
}

template<typename Real, int64_t rows, int64_t cols>
mat<Real, rows, cols> inverse(mat<Real, rows, cols> const & m)
{
    mat<Real, rows, cols> m_inv;
    if (!inverse(m, m_inv)) {
        std::cerr << "Inversion failed.\n";
    }
    return m_inv;
}

template<typename Real, int64_t rows, int64_t cols>
bool inverse(mat<Real, rows, cols> const & m, mat<Real, rows, cols>& m_inv)
{
    static_assert(rows == cols, "Must have the same number of rows as columns to compute a matrix inverse.");
    static_assert(rows < 4, "inverse not implemented for n > 3.");
    Real det = determinant(m);
    if (std::abs(det) <= std::numeric_limits<Real>::min()) {
        return false;
    }
    Real inv_det = 1/det;

    if constexpr (rows == 2) {
        m_inv(0,0) = m(1,1)*inv_det;
        m_inv(0,1) = -m(0,1)*inv_det;
        m_inv(1,0) = -m(1,0)*inv_det;
        m_inv(1,1) = m(0,0)*inv_det;
        return true;
    }
    if constexpr (rows == 3) {
        m_inv(0,0) = (m(1,1)*m(2,2) - m(1,2)*m(2,1))*inv_det;
        m_inv(0,1) = (m(0,2)*m(2,1) - m(0,1)*m(2,2))*inv_det;
        m_inv(0,2) = (m(0,1)*m(1,2) - m(0,2)*m(1,1))*inv_det;
        m_inv(1,0) = (m(1,2)*m(2,0) - m(1,0)*m(2,2))*inv_det;
        m_inv(1,1) = (m(0,0)*m(2,2) - m(0,2)*m(2,0))*inv_det;
        m_inv(1,2) = (m(0,2)*m(1,0) - m(0,0)*m(1,2))*inv_det;
        m_inv(2,0) = (m(1,0)*m(2,1) - m(1,1)*m(2,0))*inv_det;
        m_inv(2,1) = (m(0,1)*m(2,0) - m(0,0)*m(2,1))*inv_det;
        m_inv(2,2) = (m(0,0)*m(1,1) - m(0,1)*m(1,0))*inv_det;
        return true;
    }
    return false;
}

template<typename Real, int64_t rows, int64_t cols>
inline vec<Real, rows> operator*(const mat<Real, rows, cols> & M, const vec<Real, cols> &v) {
    vec<Real, rows> w;
    for (int64_t i = 0; i < rows; ++i) {
        w[i] = 0;
        for (int64_t j = 0; j < cols; ++j) {
            w[i] += M(i,j)*v[j];
        }
    }
    return w;
}

template<typename Real, int64_t rows, int64_t cols>
vec<Real, cols> mat<Real, rows, cols>::solve(vec<Real, rows> const & b) const
{
    static_assert(rows == cols, "Rows must = columns in solve.");
    static_assert(rows == 2, "Only 2x2 matrices have been implemented.");
    vec<Real, cols> v;
    if constexpr (rows == 2) {
        if (M_[2] == 0) {
            v[1] = b[1]/M_[3];
            v[0] = (b[0] - M_[1]*v[1])/M_[0];
            return v;
        }
        else {
            Real b1 = b[0] - M_[0]*b[1]/M_[2];
            v[1] = b1/(M_[1] - M_[0]*M_[3]/M_[2]);
            v[0] = b[0] - M_[1]*v[1];
            v[0] /= M_[0];
        }
    }
    else {
        std::cerr << "Not yet implemented!\n";
    }
    return v;
}

}
#endif
