#ifndef _RESTART_T_FROM_FILE_HPP
#define _RESTART_T_FROM_FILE_HPP
#include<basic_types.hpp>
#include<iostream>
#include<fstream>
#include<complex>
#include<Eigen/Dense>
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECMatrix;
bool 
restartTfromFile(REAL_TYPE&, ECMatrix &,const std::string &);
#endif
