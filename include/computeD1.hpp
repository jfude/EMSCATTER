#ifndef _COMPUTE_D1_HPP
#define _COMPUTE_D1_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<Eigen/Dense>
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
computeD1(ECMatrix &, const ECMatrix &, ECMatrix &, std::ofstream& outStream);
#endif
