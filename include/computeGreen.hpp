#ifndef _COMPUTE_GREEN_HPP
#define _COMPUTE_GREEN_HPP
#include<basic_types.hpp>
#include<vector>
#include<fstream>
#include<complex>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic> EMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic> ECMatrix;
bool computeGreen(std::vector<ECMatrix>  &,
		  std::vector<EMatrix>   &,
		  std::vector<ECMatrix>  &,
		  const REAL_TYPE, std::ofstream&);
#endif

  


