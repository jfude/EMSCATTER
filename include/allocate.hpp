/*
 * allocate.hpp -- 
 *
 */
#ifndef _ALLOCATE_HPP
#define _ALLOCATE_HPP
#include<basic_types.hpp>
#include<string>
#include<vector>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                EMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic> ECMatrix;
// 1 arg Real
bool allocate(const std::string&, 
	      std::vector<EMatrix> &,
	      const size_t,const size_t,const size_t,std::ofstream&);

// 3 arg Real
bool allocate(const std::string&, 
	      std::vector<std::vector<std::vector<EMatrix> > > &,
	      const size_t,const size_t,
	      const size_t,const size_t,const size_t,std::ofstream&);


//1 arg Complex
bool allocate(const std::string&, 
	      std::vector<ECMatrix> &,
	      const size_t,const size_t,const size_t,std::ofstream&);


// raw fundamental, Real
bool allocate(const std::string&, EMatrix &,
	      const size_t,const size_t,std::ofstream&);


// raw, fundamental, Complex
bool allocate(const std::string&, ECMatrix &,
	      const size_t,const size_t,std::ofstream&);

#endif

