#include<computeT.hpp>
bool 
computeT(std::vector<std::vector<ECMatrix> >& T, const REAL_TYPE k,
	 const std::vector<REAL_TYPE>& w,
	 const std::vector<std::vector<EMatrix> >  &J,
	 const ECMatrix &F,
	 const std::vector<AngularMultIndex> &ind)
{
    


  ECMatrix wk(2,2);
  wk.setZero();

  size_t n,np;
  size_t rsize = J.size();
  size_t nsize = ind.size();

  // Zero Tmatrix (can use for each here)
  for(n=0;n<nsize;n++) {
    for(np=0;np<nsize;np++) {
      T[n][np].setZero();
    }
  }
  
  for(size_t i_r=0;i_r<rsize;i_r++) {
    
    wk(0,0) = std::complex<REAL_TYPE>((REAL_TYPE)0.0,k*w[i_r]);
    wk(1,1)=wk(0,0);
    
    for(n=0;n<nsize;n++) {
      for(np=0;np<nsize;np++) {
	
	size_t lfc = 3*(i_r*nsize + n);
	size_t rfc = 2*np;
	

	T[n][np] +=  (wk * J[i_r][ind[n].l -1].transpose() * F.block(lfc,rfc,3,2));
	//T[n][np] +=  J[i_r][ind[n].l].transpose() * F.block(lfc,rfc,3,2);
	
      }
    }
  }
    
  return SUCCESS;
}
