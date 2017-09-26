#include<getTrapGrid.hpp>
#include<iostream>
bool
getTrapGrid(std::vector<REAL_TYPE>& v, std::vector<REAL_TYPE>& w, const size_t n,const REAL_TYPE a,const REAL_TYPE b)
{
  // Trapezoid firs

  
  size_t i;
  const REAL_TYPE TWO     = (REAL_TYPE)2.0000000000000000;
  REAL_TYPE h;
  h = (b-a)/(n-1);
  //std::cout << h << std::endl;
  v[0]=a;
  w[0]=h/TWO;
  for(i=1;i<(n-1);i++) {
    v[i] = v[i-1] + h;
    w[i] = h;
  }
  v[n-1]=b;
  w[n-1]=w[0];
  //std::cout << h << std::endl;

  return SUCCESS;
}
