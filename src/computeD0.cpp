#include<computeD0.hpp>
bool 
computeD0(ECMatrix &D0, const ECMatrix &T, ECMatrix &Q12, std::ofstream& outStream)
{
  
  // Temporarily adjust Q12
  size_t i;
  const REAL_TYPE one = 1.0000000000000000;
  const size_t nsize2 = Q12.rows();
  for(i=0;i<nsize2;i++) {
    Q12(i,i) += one;
  }
  // can transpose then undo if necessary

  D0.noalias() = T * Q12.transpose();

  for(i=0;i<nsize2;i++) {
    Q12(i,i) -= one;
  }

  return SUCCESS;
} 
    

