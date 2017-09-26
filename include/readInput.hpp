#ifndef _READ_INPUT_HPP
#define _READ_INPUT_HPP
/*
 *
 * readInputParams.hpp
 */

#include<stdio.h>
#include<stdlib.h>
#include<sstream>
#include<map>
#include<string>
#include<string.h>

#define NULL_CHAR '\0'

// check input represents a real number
bool checkReal(char*);
// check input represents an integer
bool checkInt(char*);
// parse line into individual words separated by spaces
int parseLine(FILE**,unsigned char*,int*,int*);
// read input parameters
bool readInputParams(REAL_TYPE&, REAL_TYPE&, 
	       coordParam&,coordParam&,coordParam&,size_t &,std::string&);


#endif
