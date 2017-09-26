#ifndef _COMPUTE_SCAT_INT_FUNCTIONS
#define _COMPUTE_SCAT_INT_FUNCTIONS
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<plegendre.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>    ECMatrix;
bool
computeScatIntFunctions(std::vector<REAL_TYPE>&, std::vector<REAL_TYPE>&, std::vector<REAL_TYPE>&,
			const std::vector<AngularMultIndex> &, const ECMatrix &, const ECMatrix &);
#endif

