#ifndef _COMPUTE_MATRIX_MAX_HPP
#define _COMPUTE_MATRIX_MAX_HPP
#include<basic_types.hpp>
#include<complex>
#include<Eigen/Dense>
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;

REAL_TYPE computeMatrixMax(ECMatrix &);
#endif
