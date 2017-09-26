#ifndef _COMPUTE_U_HPP
#define _COMPUTE_U_HPP
#include<Eigen/Dense>
#include<basic_types.hpp>
#include<vector>
#include<fstream>
#include<math.h>
#include<getGaussGrid.hpp>
#include<AngularMultIndex.hpp>
#include<dielectric.hpp>
#include<sqrtff.hpp>
#include<boost/math/special_functions/legendre.hpp>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic> EMatrix;
bool
computeU(std::vector<std::vector<std::vector<EMatrix> > > &,const REAL_TYPE &,
         const size_t&, const size_t&, const REAL_TYPE&,
	 const std::vector<AngularMultIndex>&,
         const std::vector<REAL_TYPE> &, const size_t &, std::ofstream&);
#endif

