#ifndef _COMPUTE_Y
#define _COMPUTE_Y
#include<basic_types.hpp>
#include<vector>
//#include<complex>
#include<AngularMultIndex.hpp>
#include<plegendre.hpp>
#include<Eigen/Dense>
//#include<boost/math/special_functions/legendre.hpp>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>   EMatrix;
bool computeY(std::vector<std::vector<std::vector<EMatrix> > > &,
	      std::vector<AngularMultIndex> &,
	      const size_t,
	      std::vector<REAL_TYPE> &,std::vector<REAL_TYPE> &,
	      std::vector<REAL_TYPE> &,std::vector<REAL_TYPE> &,
	      std::ofstream&);
#endif
