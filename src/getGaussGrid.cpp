#include<getGaussGrid.hpp>
#include<iostream>
bool
getGaussGrid(std::vector<REAL_TYPE>& v, std::vector<REAL_TYPE>& w, const size_t n,const REAL_TYPE a,const REAL_TYPE b,
	     std::ofstream& out)
{

  const REAL_TYPE TWO = (REAL_TYPE)2.00000000000000000000000L;
  // Get weights and abcissa
  if(!getGaussWA(v,w,n,out)) {return FAIL;};
  
  //for(int j=0;j<v.size();j++) {
  //  std::cout << v[j] << "    "<< w[j] << std::endl;
  //}
  
  // Adjust weights and abcissa
  REAL_TYPE diff = (b-a)/TWO;
  REAL_TYPE ave  = (b+a)/TWO;
  
  for(size_t i=0;i<n;i++) {
    v[i] = diff*v[i] + ave;
    w[i] *= diff;
  }

  return SUCCESS;
}

 
