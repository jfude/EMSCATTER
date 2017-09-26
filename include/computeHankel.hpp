#ifndef _COMPUTE_HANKEL_HPP
#define _COMPUTE_HANKEL_HPP
#include<iostream>
#include<basic_types.hpp>
#include<vector>
#include<fstream>
#include<boost/math/special_functions/hankel.hpp>
#include<Eigen/Dense>
typedef  Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic> ECMatrix;
bool computeHankel(std::vector<ECMatrix>  &,
		   const REAL_TYPE&, const REAL_TYPE&, const size_t&, std::ofstream&);
#endif

  


