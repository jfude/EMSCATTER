#include<computeC1.hpp>
bool 
computeC1(ECMatrix &C1, const ECMatrix &C0, const std::vector<ECMatrix> &G,
	  const std::vector<AngularMultIndex> &ind,
	  std::ofstream & outStream)
{
  
  
  // A number of operations can be reduced
  // Consider tmp loads to reduce thrashing

  std::complex<REAL_TYPE> one((REAL_TYPE)1.0L,(REAL_TYPE)0.0L);
  
  ECMatrix id(3,3);
  id.setIdentity(3,3);
  ECMatrix wm(3,3);

  size_t n,np;
  size_t nsize = ind.size();
    
  for(n=0;n<nsize;n++) {
    size_t lfc = 3*n;
    for(np=0;np<nsize;np++) {
   
      size_t rfc = 3*np;
      
      C1.block(lfc,rfc,3,3).noalias() = -C0.block(lfc,rfc,3,3)*G[ind[np].l -1];
      
    }
  }
  
  size_t c1size = C1.rows();
  for(n=0;n<c1size;n++) C1(n,n) += one;

  return SUCCESS;
}
