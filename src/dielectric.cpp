/**
 *  dielectric.cpp: Dielectric function of the spherical Luneburg lens
 *                  eps ~ r^{2}
 */
#include<dielectric.hpp>
std::complex<REAL_TYPE>  
dielectric(const REAL_TYPE& zc,const REAL_TYPE &phi,const REAL_TYPE &theta,const REAL_TYPE &r)
{
  const REAL_TYPE ZERO= 0.000000000000000000000L;
  const REAL_TYPE ONE=  1.000000000000000000000L;
  const REAL_TYPE TWO=  2.000000000000000000000L;
  const REAL_TYPE  P5=  0.500000000000000000000L;
  const REAL_TYPE FOUR= 4.000000000000000000000L;

  
  
  REAL_TYPE x = r*sin(theta)*cos(phi);
  REAL_TYPE y = r*sin(theta)*sin(phi);
  REAL_TYPE z = r*cos(theta);

  REAL_TYPE dz = z - zc;
  REAL_TYPE d = sqrt(x*x + y*y + dz*dz);
  
  if(d < P5) {
    return std::complex<REAL_TYPE>((TWO - FOUR*d*d),ZERO);
  }
  return std::complex<REAL_TYPE>(ONE,ZERO);

}
