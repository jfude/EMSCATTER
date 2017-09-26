
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{

  const long double ONE=1.0L;
  const long double TWO =2.0L;
  const long double wvl =1.0L;
  
  const long double PI = 3.1415926535897932384626433832795028842L;
  
  const long double TWOPI =TWO*PI;
  const long double ab_ra= TWO; 
  const long double lambda=1.0;

  long double c,b,b2,a;
  
  
  for(int i=1;i<=7;i++) {
    
    c = (long double)(i);
    
    b2 = (lambda*c)/TWOPI;
    b2 *= (b2/(ab_ra*ab_ra - ONE));

    //std::cout << b2 << std::endl;
   

    b = sqrt(b2);
    
    
    a = ab_ra * b;

    printf("%15.13Lf   %15.13Lf   %15.13Lf  %15.13Lf \n", c, lambda, a, b);
  
  }   
 
   
}
