#ifndef _SAVE_T_TO_FILE_HPP
#define _SAVE_T_TO_FILE_HPP
#include<basic_types.hpp>
#include<iostream>
#include<fstream>
#include<complex>
#include<Eigen/Dense>
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
saveTtoFile(REAL_TYPE, const ECMatrix &,const std::string&);
#endif
