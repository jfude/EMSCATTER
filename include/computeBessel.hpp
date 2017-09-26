#ifndef _COMPUTE_BESSEL_HPP
#define _COMPUTE_BESSEL_HPP
#include<basic_types.hpp>
#include<vector>
#include<fstream>
#include<boost/math/special_functions/bessel.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic> EMatrix;
bool computeBessel(std::vector<EMatrix> &,
		   const REAL_TYPE&, const REAL_TYPE&, const size_t&, std::ofstream&);
#endif

