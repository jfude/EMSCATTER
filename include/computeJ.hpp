#ifndef _COMPUTEJ_HPP
#define _COMPUTEJ_HPP
#include<basic_types.hpp>
#include<vector>
#include<fstream>
#include<boost/math/special_functions/bessel.hpp>
#include<Eigen/Dense>
bool computeJ(Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic> &,
	      const REAL_TYPE&, const std::vector<REAL_TYPE>&, const size_t&, std::ofstream&);
#endif

  


