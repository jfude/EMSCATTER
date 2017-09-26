#include<computeMatrixMax.hpp>

REAL_TYPE computeMatrixMax(ECMatrix &Q) {

  REAL_TYPE qMax=(REAL_TYPE)0.00000000000000000000000L;
  REAL_TYPE tMax;
  const size_t nrows = Q.rows();
  const size_t ncols = Q.cols();

  for(size_t i=0;i<nrows;i++) {
    for(size_t j=0;j<ncols;j++) {
      tMax = std::norm(sqrt(Q(i,j)));
      if(qMax < tMax) {qMax=tMax;}
    }
  }

  return qMax;
}


    
  
