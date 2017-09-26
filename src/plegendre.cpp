#include<plegendre.hpp>
REAL_TYPE plegendre(const int l, const int m, const REAL_TYPE x) 
{
  
  const REAL_TYPE ZERO   = (REAL_TYPE)0.00000000000000000000000000L;
  const REAL_TYPE PI     = (REAL_TYPE)3.1415926535897932384626433832795028842L;
  const REAL_TYPE ONE    = (REAL_TYPE)1.00000000000000000000000000L;
  const REAL_TYPE TWO    = (REAL_TYPE)2.00000000000000000000000000L;
  const REAL_TYPE THREE  = (REAL_TYPE)3.00000000000000000000000000L;
  const REAL_TYPE FOUR   = (REAL_TYPE)4.00000000000000000000000000L;
  const REAL_TYPE FOURPI = FOUR*PI;

  int i,ll;
  REAL_TYPE fact,oldfact,pll,pmm,pmmp1,omx2;

  if(m > l) {return ZERO;}

  pmm=ONE;
  if(m>0) {
    omx2 = (ONE-x)*(ONE+x);
    fact=ONE;
    for(i=1;i<=m;i++) {
      pmm *= omx2*fact/(fact+ONE);
      fact += TWO;
    }
  }
  pmm = sqrt((TWO*m+ONE)*pmm/FOURPI);
  
  if(m & 1)
    pmm=-pmm;
  
  if(l==m) {
    return pmm;
  }
  else {
    pmmp1 = x*sqrt(TWO*m + THREE)*pmm;
    if(l == (m+1))
      return pmmp1;
    else {
      oldfact=sqrt(TWO*m + THREE);
      for(ll=(m+2);ll<=l;ll++) {
	fact=sqrt( (FOUR*ll*ll -ONE)/(ll*ll - m*m));
	pll=(x*pmmp1 - pmm/oldfact)*fact;
	oldfact=fact;
	pmm=pmmp1;
	pmmp1=pll;
      }
      return pll;
    }
  }

}
    
  
