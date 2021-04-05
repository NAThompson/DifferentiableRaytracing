#ifndef DRT_MATRIX_HPP
#define DRT_MATRIX_HPP
#include <cmath>
#include <array>
#include <iostream>
#include <initializer_list>
#include <drt/vec.hpp>
#include <drt/math.hpp>

namespace drt {

template<typename Real, int64_t rows = 3, int64_t cols = 3>
class matrix {
public:

    matrix() {};

    Real operator()(int64_t i, int64_t j) const {
        return M_[i*cols + j];
    }

    Real& operator()(int64_t i, int64_t j) {
        return M_[i*cols + j];
    }

    vec<Real, rows> column(int64_t j) const {
        vec<Real, rows> c;
        for (int64_t i = 0; i < rows; ++i) {
            c[i] = M_[i*cols + j];
        }
        return c;
    }

    // Return an x that satisfies Ax = b.
    vec<Real, cols> solve(vec<Real, rows> const & b) const;

    friend std::ostream& operator<<(std::ostream & os, const matrix<Real, rows, cols> & v)
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

    
    Real max_norm() const {
        using std::abs;
        Real m = abs(M_[0]);
        for (int64_t i = 1; i < rows*cols; ++i) {
            if (abs(M_[i]) > m) {
                m = abs(M_[i]);
            }
        }
        return m;
    }


private:
    Real M_[rows*cols];
};

template<typename Real, int64_t rows, int64_t cols>
inline Real determinant(matrix<Real, rows, cols> const & m)
{
    static_assert(rows == cols, "Must have the same number of rows as columns to compute a determinant.");
    static_assert(rows < 4, "Determinant not implemented for n > 3.");
    if constexpr (rows == 2) {
        return difference_of_products(m(0,0), m(1,1), m(0,1), m(1,0));
    }
    if constexpr (rows == 3) {
        return m(0,0)*(m(1,1)*m(2,2) - m(1,2)*m(2,1)) - m(0,1)*(m(1,0)*m(2,2) - m(1,2)*m(2,0)) + m(0,2)*(m(1,0)*m(2,1) - m(1,1)*m(2,0));

    }
    return std::numeric_limits<Real>::quiet_NaN();
}

template<typename Real, int64_t rows, int64_t cols>
matrix<Real, rows, cols> inverse(matrix<Real, rows, cols> const & m)
{
    matrix<Real, rows, cols> m_inv;
    if (!inverse(m, m_inv)) {
        std::cerr << "Inversion failed.\n";
    }
    return m_inv;
}

template<typename Real, int64_t rows, int64_t cols>
bool inverse(matrix<Real, rows, cols> const & m, matrix<Real, rows, cols>& m_inv)
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
inline vec<Real, rows> operator*(const matrix<Real, rows, cols> & M, const vec<Real, cols> &v) {
    vec<Real, rows> w;
    for (int64_t i = 0; i < rows; ++i) {
        w[i] = 0;
        for (int64_t j = 0; j < cols; ++j) {
            w[i] += M(i,j)*v[j];
        }
    }
    return w;
}

template<typename Real, int64_t rows, int64_t cols, int64_t cols2>
inline matrix<Real, rows, cols2> operator*(const matrix<Real, rows, cols> & A, const matrix<Real, cols, cols2> &B) {
    matrix<Real, rows, cols2> M;
    for (int64_t i = 0; i < rows; ++i) {
        for (int64_t j = 0; j < cols2; ++j)
        {
            M(i,j) = 0;
            for (int64_t k = 0; k < cols; ++k) {
                M(i,j) += A(i,k)*B(k,j);
            }
        }
    }
    return M;
}

template<typename Real, int64_t rows, int64_t cols>
vec<Real, cols> matrix<Real, rows, cols>::solve(vec<Real, rows> const & b) const
{
    static_assert(rows == cols, "Rows must = columns in solve.");
    static_assert(rows == 2 || rows == 3, "Only 2x2 or 3x3 matrices have been implemented.");
    vec<Real, cols> v;
    if constexpr (rows == 2) {
        // M(0,0) = M_[0] M(0,1) = M_[1]
        // M(1,0) = M_[2] M(1,1) = M_[3]
        // Ugh this is the real condition but . . .
        //if (M_[2]) == 0) {
        // How is this scale determined? By the failing unit test, which in turn was determined
        // by a bizarre problem found in the raytracer!!!
        if (std::abs(M_[2]) <= std::numeric_limits<Real>::epsilon()*std::numeric_limits<Real>::epsilon()) {
            v[1] = b[1]/M_[3];
            v[0] = (b[0] - M_[1]*v[1])/M_[0];
            return v;
        }
        else {
            Real rM10 = 1/M_[2];
            Real m00dm10 = M_[0]*rM10;
            Real b1 = b[0] - b[1]*m00dm10;
            v[1] = b1/(M_[1] - M_[3]*m00dm10);
            v[0] = b[1] - M_[3]*v[1];
            v[0] *= rM10;
            return v;
        }
    }
    else if constexpr (rows == 3) {
        matrix<Real, 3,3> A;
        vec<Real, 3> y;
        if (M_[0] != 0) {
            A(0,0) = 1;
            A(0,1) = M_[1]/M_[0];
            A(0,2) = M_[2]/M_[0];
            y[0] = b[0]/M_[0];
            // M(1,0) = M_[3]
            if(M_[3] != 0) {
                A(1,0) = 0;
                A(1,1) = A(0,1) - M_[4]/M_[3];
                A(1,2) = A(0,2) - M_[5]/M_[3];
                y[1] = y[0] - b[1]/M_[3];
            }
            // M(2,0) = M_[6]
            if (M_[6] != 0) {
                A(2,0) = 0;
                A(2,1) = A(0,1) - M_[7]/M_[6];
                A(2,2) = A(0,2) - M_[8]/M_[6];
                y[2] = y[0] - b[2]/M_[6];
            }

            // Now A and y are fully populated.
            if (A(1,1) != 0) {
                A(1,2) = A(1,2)/A(1,1);
                y[1] = y[1]/A(1,1);
                A(1,1) = 1;

                A(2,2) = A(1,2) - A(2,2)/A(2,1);
                y[2] = y[1] - y[2]/A(2,1);
                A(2,1) = 0;

                if (A(2,2) == 0) {
                    //std::cerr << __FILE__ << ":" << __LINE__ << " 3x3 matrix is singular.\n";
                }
                y[2] /= A(2,2);
                v[2] = y[2];
                v[1] = y[1] - A(1,2)*v[2];
                v[0] = y[0] - A(0,1)*v[1] - A(0,2)*v[2];
                return v;
            }
            else {
                if (A(1,2) == 0) {
                    //std::cerr << __FILE__ << ":" << __LINE__ << " 3x3 matrix is singular.\n";
                }
                if (A(2,1) == 0) {
                    //std::cerr << __FILE__ << ":" << __LINE__ << " 3x3 matrix is singular.\n";
                }

                v[2] = y[1]/A(1,2);
                v[1] = (y[2] - A(2,2)*v[2])/A(2,1);
                v[0] = (y[0] - A(0,1)*v[1] - A(0,2)*v[2]);
            }
        }
        else {
            // M(0,0) = 0, so make it lower triangular.
            A(0,0) = 0;
            if (M_[6] != 0) {
                A(2,1) = M_[7]/M_[6];
                A(2,2) = M_[8]/M_[6];
                y[2] = b[2]/M_[6];
                A(2,0) = 1;
                if (M_[3] != 0) {
                    A(1,0) = 0;
                    A(1,1) = A(2,1) - M_[4]/M_[3];
                    A(1,2) = A(2,2) - M_[5]/M_[3];
                    y[1] = y[2] - b[1]/M_[3];
                }
                if (A(1,1) == 0) {
                    std::cerr << __FILE__ << ":" << __LINE__ <<  ": All hell just broke loose; A(1,1) = 0.\n";
                }
                else {
                    A(1,2) = A(1,2)/A(1,1);
                    y[1] = y[1]/A(1,1);
                    A(1,1) = 1;
                }
                A(0,1) = 0;
                A(0,2) = A(1,2) - M_[2]/M_[1];
                y[0] = y[1] -b[0]/M_[1];
                if (A(0,2) == 0) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " Matrix is singular.\n";
                }
                y[0] = y[0]/A(0,2);
                A(0,2) = 1;
                v[2] = y[0];
                v[1] = y[1] - A(1,2)*v[2];
                v[0] = y[2] - A(2,1)*v[1] - A(2,2)*v[2];
                return v;
            }
            else {
                std::cerr << __FILE__ << ":" << __LINE__ <<  ": All hell just broke loose; M(2,0) = 0.\n";
            }

        }
    }
    else {
        std::cerr << "Not yet implemented!\n";
    }
    return v;
}

}
#endif
