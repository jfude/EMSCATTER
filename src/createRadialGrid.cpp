/**
 *
 * createRadialGrid.cpp
 */ 

#include<createRadialGrid.hpp>

bool createRadialGrid(const REAL_TYPE& r1,const REAL_TYPE& r2,const int& n,std::vector<REAL_TYPE>& r,
		      std::vector<REAL_TYPE>& w,std::ofstream& outStream) 
{
  
  // make msg 
  if(r2<r1)  return FAIL;
  if(n <= 2) return FAIL;  
  try {
    r.resize(n);
    w.resize(n);
  }
  catch(std::bad_alloc const&) {
    outStream << "createRadialGrid:"<<std::endl;
    outStream << "Allocation of memory for Radial Grid or Weights failed."<<std::endl;
    return FAIL;
  }

  
  REAL_TYPE p5  = (REAL_TYPE)(0.50);
  REAL_TYPE one = (REAL_TYPE)(1.00);
  
  REAL_TYPE dr = (r2-r1)/n;
  r[0]=r1;
  w[0]=p5;
  for(int i=1;i<n;i++) {
    r[i] = r[i-1] + dr;
    w[i] = one;
  }
  w[n-1]=p5;
  
  return SUCCESS;
}
  
