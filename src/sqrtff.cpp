#include<sqrtff.hpp>
/**
 *  sqrtff: Computes the factor sqrt((l-m)!/(l+m)!)
 *
 **/
REAL_TYPE  sqrtff(const int l,const int m) {
    
  long double one=1.0000000000000000000L;
  long double result=one;
  int n,nstart,nend;

  if(m==0) return  (REAL_TYPE)one;
  
  nstart = l-m+1;
  nend   = l+m;
  for(n=nstart;n <= nend;n++) result *= (one/sqrt((long double)(n)));

  return (REAL_TYPE)(result);
  
}
