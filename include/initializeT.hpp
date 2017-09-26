#ifndef _INITIALIZE_T_HPP
#define _INITIALIZE_T_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool initializeT(ECMatrix &); 
#endif
