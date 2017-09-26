#ifndef _READ_INPUT_PARAMS_HPP
#define _READ_INPUT_PARAMS_HPP
/*
 *
 * readInputParams.hpp
 */

#include<basic_types.hpp>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<string.h>
#include<sstream>
#include<fstream>
#include<map>
#include<string>
#include<coordParam.hpp>

#define NULL_CHAR '\0'

// check input represents a real number
bool checkReal(char*);
// check input represents an integer
bool checkInt(char*);
// writeTypeCheckError
void writeTypeCheckError(const std::string&,char *,std::ofstream&);
// parse line into individual words separated by spaces
int parseLine(FILE**,unsigned char*,int*,int*);
// read input parameters
bool readInputParams(std::string&, REAL_TYPE&, REAL_TYPE&, 
		     coordParam&,coordParam&,coordParam&,size_t &,
		     std::string&, const std::string&, std::ofstream&);


#endif
