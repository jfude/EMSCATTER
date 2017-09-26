/**
 *  dielectric.cpp: Dielectric function for a single uniform cylinder
 *                  
 *
 */
#include<dielectric.hpp>
std::complex<REAL_TYPE> 
dielectric(const REAL_TYPE &phi,const REAL_TYPE &theta,const REAL_TYPE &r)
{

  // Constants
  const REAL_TYPE ZERO = (REAL_TYPE)0.000000000000000000000000L;
  const REAL_TYPE P5   = (REAL_TYPE)0.5000000000000000000000000L;
  const REAL_TYPE P25  = (REAL_TYPE)0.25000000000000000000000000L;
  const REAL_TYPE ONE=  (REAL_TYPE) 1.000000000000000000000000L;
  const REAL_TYPE ONEP5=  (REAL_TYPE) 1.500000000000000000000000L;
  const REAL_TYPE THREEP5=  (REAL_TYPE) 3.500000000000000000000000L;
  const REAL_TYPE FIVEP5=  (REAL_TYPE) 5.500000000000000000000000L;
  const REAL_TYPE TWO=  (REAL_TYPE) 2.000000000000000000000000L;
  const REAL_TYPE FOUR= (REAL_TYPE) 4.000000000000000000000000L;
 
  
  const REAL_TYPE zc =  ZERO;  // shift along z
  const REAL_TYPE radius =P5; // Radius of cylinder
  const REAL_TYPE height = FIVEP5; // half height
  const std::complex<REAL_TYPE> dielectricConstant= std::complex<REAL_TYPE>(1.7161L,0.0L);

  // REAL_TYPE x = r*sin(theta)*cos(phi);
  // REAL_TYPE y = r*sin(theta)*sin(phi);
  REAL_TYPE z = r*cos(theta);
  REAL_TYPE dz = z-zc;
  //REAL_TYPE dr = sqrt(x*x + y*y);
  REAL_TYPE dr = r*std::abs(sin(theta));
  
  if((dr < radius) && (dz> (-height)) && (dz < height) ) return dielectricConstant;
  return std::complex<REAL_TYPE>(ONE,ZERO);

}
