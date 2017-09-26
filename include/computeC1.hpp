#ifndef _COMPUTE_C1_HPP
#define _COMPUTE_C1_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                   EMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
computeC1(ECMatrix &, const ECMatrix &, const std::vector<ECMatrix> &,
	  const std::vector<AngularMultIndex> &, std::ofstream & outStream);
#endif
