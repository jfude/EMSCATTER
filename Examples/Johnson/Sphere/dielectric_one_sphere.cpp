/**
 *  dielectric.cpp: Dielectric function for one  uniform sphere
 *
 */
#include<dielectric.hpp>
std::complex<REAL_TYPE> 
dielectric(const REAL_TYPE& zc, const REAL_TYPE &phi,const REAL_TYPE &theta,const REAL_TYPE &r)
{

  // Constants
  const REAL_TYPE ZERO= (REAL_TYPE)0.000000000000000000000000L;
  const REAL_TYPE P5  = (REAL_TYPE)0.500000000000000000000000L;
  const REAL_TYPE ONE  = (REAL_TYPE)1.00000000000000000000000L;
  
  
  // zc = zshift
  const REAL_TYPE radius = P5; // Radius of sphere
  
  
  const std::complex<REAL_TYPE> dielectricConstant = std::complex<REAL_TYPE>(2.00L,0.000L);
  
  
  REAL_TYPE xy = r*sin(theta);
  
  REAL_TYPE z = r*cos(theta);
  REAL_TYPE dz = z-zc;
  REAL_TYPE d = sqrt(xy*xy + dz*dz);
  
  if(d < radius) return dielectricConstant;

  return std::complex<REAL_TYPE>(ONE,ZERO);

}
