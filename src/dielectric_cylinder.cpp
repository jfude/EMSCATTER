/**
 *  dielectric.cpp: Dielectric function for a single uniform cylinder
 *                  Long axis (height axis) is along the z axis
 *                  
 *                  
 *
 */
#include<dielectric.hpp>
std::complex<REAL_TYPE> 
dielectric(const REAL_TYPE &zc, const REAL_TYPE &phi,const REAL_TYPE &theta,const REAL_TYPE &r)
{

  // Constants
  const REAL_TYPE ZERO = (REAL_TYPE)0.000000000000000000000000L;
  const REAL_TYPE P5   = (REAL_TYPE)0.5000000000000000000000000L;
  const REAL_TYPE ONE=  (REAL_TYPE) 1.000000000000000000000000L;
  const REAL_TYPE TWO=  (REAL_TYPE) 2.000000000000000000000000L;
  
  
  //const REAL_TYPE zc =  ZERO;  // shift along z
  const REAL_TYPE radius =3.978873577297384; // Radius of cylinder
  const REAL_TYPE height =TWO*radius; // half height
  const std::complex<REAL_TYPE> dielectricConstant= std::complex<REAL_TYPE>(2.340836L,0.02448L);

  // REAL_TYPE x = r*sin(theta)*cos(phi);
  // REAL_TYPE y = r*sin(theta)*sin(phi);
  REAL_TYPE z = r*cos(theta);
  REAL_TYPE dz = z-zc;
  //REAL_TYPE dr = sqrt(x*x + y*y);
  REAL_TYPE dr = r*abs(sin(theta));
  
  if((dr < radius) && (dz> (-height)) && (dz < height) ) return dielectricConstant;
  return std::complex<REAL_TYPE>(ONE,ZERO);

}
