#include<computeQ12.hpp>
bool 
computeQ12(ECMatrix &Q12, ECMatrix &Q, std::vector<EMatrix> &J,
	   std::vector<ECMatrix> &H, const std::vector<AngularMultIndex> &ind,
	   const REAL_TYPE k,std::ofstream &outStream)
{
  
  
  // A number of operations can be reduced
  // Consider tmp loads to reduce thrashing

  std::complex<REAL_TYPE> i_k((REAL_TYPE)0.0,k);
  
  //ECMatrix id(3,3);
  //id.setIdentity(3,3);

  size_t n,np;
  size_t nsize = ind.size(); 
  
  for(n=0;n<nsize;n++) {
    size_t lfc_2 = 2*n;
    size_t lfc_3 = 3*n;
    for(np=0;np<nsize;np++) {
      size_t rfc_2 = 2*np;
      size_t rfc_3 = 3*np;
      
      Q12.block(lfc_2,rfc_2,2,2) = 
	i_k * J[ind[n].l -1].transpose() * Q.block(lfc_3,rfc_3,3,3) * H[ind[np].l -1];
      
    }
  }

  return SUCCESS;
}
