#include<sqrt_factorial.hpp>
REAL_TYPE  sqrt_factorial(const int n) {
  
  // need to exc handling
  const int fac[13] = {1,1,2,6,24,120,720,5040,40320,
		       362880,3628800,39916800,479001600};
  
  return sqrt((REAL_TYPE)(fac[n]));
  
}
