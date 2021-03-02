#ifndef DRT_TENSOR_HPP
#define DRT_TENSOR_HPP
#include <cmath>
#include <array>
#include <iostream>
#include <initializer_list>

namespace drt {

template<typename Real, int64_t panels = 3, int64_t rows = 3, int64_t cols = 3>
class tensor {
public:

    tensor() {};

    inline Real operator()(int64_t i, int64_t j, int64_t k) const {
        return T_[i*rows*cols + j*cols + k];
    }

    inline Real& operator()(int64_t i, int64_t j, int64_t k) {
        return T_[i*rows*cols + j*cols + k];
    }

    // Tensor contraction. T is a bilinear map from ℝⁿxℝⁿ -> ℝⁿ so that T(x,y) ∈ ℝⁿ.
    vec<Real, panels> operator()(vec<Real, cols> const & x, vec<Real, rows> const & y) const {
        vec<Real, panels> z;
        for (int64_t i = 0; i < panels; ++i) {
            z[i] = 0;
            for (int64_t j = 0; j < cols; ++j) {
                for (int64_t k = 0; k < rows; ++k) {
                    z[i] += T_[i*rows*cols + j*cols + k]*x[j]*y[k];
                }
            }
        }
        return z;
    }

    friend std::ostream& operator<<(std::ostream & os, const tensor<Real, panels, rows, cols> & v)
    {
        for (int64_t i = 0; i < panels; ++i) {
            os << "{";
            for (int64_t j = 0; j < rows; ++j) {
                if (j != 0) { os << " "; }
                os << "[";
                for (int64_t k = 0; k < cols - 1; ++k) {
                    os << v.T_[i*rows*cols + j*cols + k] << ", ";
                }
                os << v.T_[i*rows*cols + j*cols + cols - 1] << "]";
                if (j != rows - 1) {
                    os << "\n";
                }
            }
            os << "}\n";
        }
        return os;
    }

private:
    Real T_[rows*cols*panels];
};


}
#endif
