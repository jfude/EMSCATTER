#ifndef _COMPUTE_D0_HPP
#define _COMPUTE_D0_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<Eigen/Dense>
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
computeD0(ECMatrix &D0, const ECMatrix &T, ECMatrix &Q12, std::ofstream& outStream);
#endif
