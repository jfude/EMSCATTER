/**
 *  dielectric.cpp: Dielectric function for three uniform spheres
 *
 */
#include<dielectric.hpp>
std::complex<REAL_TYPE> 
dielectric(const REAL_TYPE &phi,const REAL_TYPE &theta,const REAL_TYPE &r)
{

  // Constants
  const REAL_TYPE ZERO= (REAL_TYPE)0.000000000000000000000000L;
  const REAL_TYPE ONE=  (REAL_TYPE) 1.000000000000000000000000L;
  const REAL_TYPE TWO=  (REAL_TYPE) 2.000000000000000000000000L;
  const REAL_TYPE FOUR= (REAL_TYPE) 3.000000000000000000000000L;
 
  
  const REAL_TYPE zc =  ZERO;  // Center of sphere
  const REAL_TYPE radius1 =ONE; // Radius of sphere
  const REAL_TYPE radius2 =TWO; // Radius of sphere
  const REAL_TYPE radius3 =THREE; // Radius of sphere
  
  const std::complex<REAL_TYPE> dielectricConstant1= std::complex<REAL_TYPE>(4.25L,0.000L);
  const std::complex<REAL_TYPE> dielectricConstant2= std::complex<REAL_TYPE>(3.25L,0.000L);
  const std::complex<REAL_TYPE> dielectricConstant3= std::complex<REAL_TYPE>(2.25L,0.000L);
  

  REAL_TYPE x = r*sin(theta)*cos(phi);
  REAL_TYPE y = r*sin(theta)*sin(phi);
  REAL_TYPE z = r*cos(theta);
  REAL_TYPE dz = z-zc;
  REAL_TYPE d = sqrt(x*x + y*y + dz*dz);
  
  if(d < radius1) return dielectricConstant1;
  if(d < radius2) return dielectricConstant2;
  if(d < radius3) return dielectricConstant3;
  
  return std::complex<REAL_TYPE>(ONE,ZERO);

}
