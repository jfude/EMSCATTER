//#include<computeQdd.hpp>
#ifndef _COMPUTE_QDD_HPP
#define _COMPUTE_QDD_HPP
#include<basic_types.hpp>
#include<vector>
#include<complex>
#include<AngularMultIndex.hpp>
#include<Eigen/Dense>
typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                   EQMatrix;
typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic>     ECQMatrix;
template<typename JMatrix>
bool
computeQdd(ECQMatrix &Qdd, const ECQMatrix &Q, std::vector<JMatrix> &J,
	   const std::vector<AngularMultIndex> &ind,
	   const REAL_TYPE k,std::ofstream & outStream)
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
   

      Qdd.block(lfc_2,rfc_2,2,2) = i_k * J[ind[n].l -1].transpose() * Q.block(lfc_3,rfc_3,3,3) * J[ind[np].l -1];
      
    }
  }

  return SUCCESS;
}

#endif
// Explicit instantiations

/*
template<std::vector<EMatrix> > bool computeQdd(ECMatrix &, const ECMatrix &, std::vector<EMatrix> &,
				  const std::vector<AngularMultIndex> &,
				  const REAL_TYPE,std::ofstream &);

template<std::vector<ECMatrix> > bool computeQdd(ECMatrix &, const ECMatrix &, std::vector<ECMatrix> &,
				  const std::vector<AngularMultIndex> &,
				  const REAL_TYPE,std::ofstream &);

*/
