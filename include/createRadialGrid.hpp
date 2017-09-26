/**
 *
 * createRadialGrid.hpp
 *
 */
#ifndef _CREATE_RADIAL_GRID_HPP
#define _CREATE_RADIAL_GRID_HPP
#include<basic_types.hpp>
#include<iostream>
#include<vector>
#include<fstream>
bool createRadialGrid(const REAL_TYPE&,const REAL_TYPE&,const int&,
		      std::vector<REAL_TYPE>&,std::vector<REAL_TYPE>&,std::ofstream&);
#endif
