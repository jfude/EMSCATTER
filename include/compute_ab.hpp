#ifndef _COMPUTE_AB_HPP
#define _COMPUTE_AB_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<Eigen/Dense>
//typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                   EMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
compute_ab(ECMatrix &,ECMatrix &,const ECMatrix &,const std::vector<AngularMultIndex> &);
#endif
