/*
 * computeGreen.cpp -- compose Greens function from Bessel and Hankel matrices
 *
 *
 */ 
#include<computeGreen.hpp>
bool
computeGreen(std::vector<ECMatrix> &G,
	     std::vector<EMatrix>  &J,
	     std::vector<ECMatrix> &H,
	     const REAL_TYPE k,
	     std::ofstream& outStream)
{
  size_t i,j,l; 
  const size_t lsize=J.size();
  // can check J and H are same size
  //const std::complex<REAL_TYPE> ik((REAL_TYPE)0.0,k);
  const std::complex<REAL_TYPE> ik2((REAL_TYPE)0.0,0.50*k);
  
  for(l=0;l<lsize;l++) {
    G[l] = ik2 * (H[l] * J[l].transpose()  +  J[l] * H[l].transpose());
  }
  
  // no exceptions here yet
  return SUCCESS;
  
}

