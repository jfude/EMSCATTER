/**
 *
 * computeJ.cpp -- Fills J matrices with spherical bessel functions
 *              -- Requires that the radial vector,r, is sorted
 */
#include<computeJ.hpp>
bool
computeJ(Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>  &J,
	 const REAL_TYPE& k, const std::vector<REAL_TYPE> & rv, 
	 const size_t& lmax, std::ofstream& outStream)
{
  
  REAL_TYPE j_l,j_l_m1,j_l_p1;
  const size_t rsize=rv.size();
  
  // Check J has the correct shape 
  
  try {
    J.setZero(J.rows(),J.cols());
  }
  catch(...) {
    //outStream << "computeJ: Failure in initial zeroing of J matrix "<< std::endl;
    return FAIL;
  }
  
  
  
  for(size_t j=0,i=0;i<rsize;i++) {  
    REAL_TYPE r = rv[i]*k;
    for(size_t l=1;l<=lmax;l++, j+=3) {
      
      REAL_TYPE l2_1c= (REAL_TYPE)(2*l +1);
      REAL_TYPE l_1c = (REAL_TYPE)(l+1)/l2_1c;
      REAL_TYPE l_c  = (REAL_TYPE)(l)/l2_1c;
      
      // Need to capture exceptions from Boost calls
      try {
	j_l    = boost::math::sph_bessel(l,r);
	j_l_m1 = boost::math::sph_bessel(l-1,r);
	j_l_p1 = boost::math::sph_bessel(l+1,r);
      }
      catch(...) {
	//outStream << "computeJ: Failure in calculating a spherical bessel"<< std::endl;
	//outStream << "          l= "<<l<< std::endl;
	return FAIL;
      }
      
      J(j  ,0)   = j_l;
      J(j+1,1)   = l_1c*j_l_m1  -  l_c*j_l_p1;
      J(j+2,1)   = (sqrt((REAL_TYPE)(l*(l+1)))/l2_1c) * (j_l_m1 + j_l_p1);
    }
  }
  
  return SUCCESS;
}


      
      
  
