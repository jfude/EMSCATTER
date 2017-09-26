/**
 *
 * computeBessel.cpp -- Fills J or H matrices with spherical bessel or hankel functions
 *              -- Requires that the radial vector,r, is sorted
 */
#include<computeHankel.hpp>
bool
computeHankel(std::vector<ECMatrix> &J,const REAL_TYPE& k, const REAL_TYPE& r, 
	      const size_t& lmax, std::ofstream& outStream)
{
  size_t i,l;
  std::complex<REAL_TYPE> j_l,j_l_m1,j_l_p1;
  
  // Check J has the correct shape 
  
  try {
    for(l=0;l<lmax;l++) {
      J[l].setZero();
    }
  }
  catch(...) {
    //outStream << "computeJ: Failure in initial zeroing of J matrix "<< std::endl;
    std::cout << " Failure in initial zeroing of H matrix "<< std::endl;
    return FAIL;
  }
  
  REAL_TYPE rk = r*k;
  for(l=1;l<=lmax;l++) {
    
    REAL_TYPE l2_1c= (REAL_TYPE)(2*l +1);
    REAL_TYPE l_1c = (REAL_TYPE)(l+1)/l2_1c;
    REAL_TYPE l_c  = (REAL_TYPE)(l)/l2_1c;
    
    // Need to capture exceptions from Boost calls
    // can get all of these in the outer loop!!
    try {
      j_l    = boost::math::sph_hankel_1(l,rk);
      j_l_m1 = boost::math::sph_hankel_1(l-1,rk);
      j_l_p1 = boost::math::sph_hankel_1(l+1,rk);
    }
    catch(...) {
      //outStream << "computeJ: Failure in calculating a spherical bessel"<< std::endl;
      //outStream << "          l= "<<l<< std::endl;
      std::cout << "Failed calc of h" << std::endl;
      return FAIL;
    }
    
    //J(j  ,0)   = j_l;
    //J(j+1,1)   = l_1c*j_l_m1  -  l_c*j_l_p1;
    //J(j+2,1)   = (sqrt((REAL_TYPE)(l*(l+1)))/l2_1c) * (j_l_m1 + j_l_p1);
    J[l-1](0,0)   = j_l;
    J[l-1](1,1)   = l_1c*j_l_m1  -  l_c*j_l_p1;
    J[l-1](2,1)   = (sqrt((REAL_TYPE)(l*(l+1)))/l2_1c) * (j_l_m1 + j_l_p1);
  }
    
  return SUCCESS;
}


      
      
  
