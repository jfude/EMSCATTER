#include<computeD1.hpp>
bool 
computeD1(ECMatrix &D1, const ECMatrix &T, ECMatrix &Q22, std::ofstream& outStream)
{
  
  // Temporarily adjust Q12
  size_t i;
  const REAL_TYPE one = 1.0000000000000000;
  const size_t nsize2 = Q22.rows();
  
  D1.noalias() = -T * Q22.transpose();

  for(i=0;i<nsize2;i++) {
    D1(i,i) += one;
  }

  return SUCCESS;
} 
