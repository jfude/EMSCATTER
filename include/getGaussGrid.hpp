#ifndef _GET_GAUSS_GRID
#define _GET_GAUSS_GRID
#include<basic_types.hpp>
#include<vector>
#include<fstream>
#include<getGaussWA.hpp>
bool
getGaussGrid(std::vector<REAL_TYPE>&, std::vector<REAL_TYPE>&, const size_t, const REAL_TYPE, const REAL_TYPE,
	     std::ofstream&);
#endif
