#ifndef _COMPUTE_T_HPP
#define _COMPUTE_T_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                   EMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
computeT(std::vector<std::vector<ECMatrix> >&, const REAL_TYPE,
	 const std::vector<REAL_TYPE> &,
	 const std::vector<std::vector<EMatrix> > &,
	 const ECMatrix &, 
	 const std::vector<AngularMultIndex> &);
#endif
