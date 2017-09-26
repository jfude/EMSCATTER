#ifndef _COMPUTE_Q12_HPP
#define _COMPUTE_Q12_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                   EMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
computeQ12(ECMatrix &, ECMatrix &, std::vector<EMatrix> &, std::vector<ECMatrix> &,
	  const std::vector<AngularMultIndex> &, const REAL_TYPE, std::ofstream & outStream);
#endif
