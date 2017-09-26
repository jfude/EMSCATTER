#ifndef _COMPUTE_C0_HPP
#define _COMPUTE_C0_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<dielectric.hpp>
#include<sqrtff.hpp>
#include<boost/math/special_functions/legendre.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                   EMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
computeC0(ECMatrix &, std::vector<EMatrix> &, const std::vector<AngularMultIndex> &,
	  const std::vector<REAL_TYPE> &,const std::vector<REAL_TYPE> &,
	  const std::vector<REAL_TYPE> &,const std::vector<REAL_TYPE> &,
	  const REAL_TYPE, const REAL_TYPE, const REAL_TYPE, std::ofstream&);	  
#endif
