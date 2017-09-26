#include<compute_ab.hpp>
bool 
compute_ab(ECMatrix &a, ECMatrix &b, const ECMatrix &T,
	   const std::vector<AngularMultIndex> &ind)
{
  
  const REAL_TYPE ZERO = (REAL_TYPE)(0.000000000000000);
  const REAL_TYPE ONE  = (REAL_TYPE)(1.000000000000000);
  const REAL_TYPE TWO  = (REAL_TYPE)(2.000000000000000);
  
  
  const std::complex<REAL_TYPE> i_lpl[4] = { std::complex<REAL_TYPE>(ONE,ZERO),
					     std::complex<REAL_TYPE>(ZERO,ONE),
					     std::complex<REAL_TYPE>(-ONE,ZERO),
					     std::complex<REAL_TYPE>(ZERO,-ONE)};
  
  const size_t nsize = ind.size();
  
  // try to do with const fixed size
  // Eigen::Vector<std::complex<REAL_TYPE>,Eigen::Dynamic> i_vx(2),i_vy(2);
  ECMatrix i_v(2,2),sum(2,2);
  i_v(0,0)=std::complex<REAL_TYPE>(ONE,   ZERO);
  i_v(0,1)=std::complex<REAL_TYPE>(ONE,   ZERO);
  i_v(1,0)=std::complex<REAL_TYPE>(ZERO, -ONE);
  i_v(1,1)=std::complex<REAL_TYPE>(ZERO,  ONE);
    

  //const complex<REAL_TYPE> i_v(ONE,-ONE);
  //ECMatrix i_v(2,2),sum(2,2); // make this defined at compile time <2,2>
  // i_v(0,0) = complex<REAL_TYPE>(ONE,  ZERO);
  // i_v(0,1) = complex<REAL_TYPE>(ONE.  ZERO);
  // i_v(1,0) = complex<REAL_TYPE>(ZERO,-ONE);
  // i_v(1,1) = complex<REAL_TYPE>(ZERO. ONE);


  // Initialize coefficients 
  a.setZero();b.setZero();

  int lp_l,p;
  size_t n,np;
  REAL_TYPE ratio;
  for(n=0;n<nsize;n++) {
    size_t lfc=2*n;
    sum.setZero(2,2);
    for(np=0;np<nsize;np++) {
      
      if(ind[np].m==1) {
	size_t rfc=2*np;
	lp_l = ind[np].l - ind[n].l;
	
	// std::cout << lp_l;
	if(lp_l<0) {
	  lp_l += ((-lp_l)/4  +1)*4;
	}
	
	// std::cout << "   " << lp_l%4 << std::endl;

	ratio = sqrt((TWO*ind[np].l + ONE)/(TWO*ind[n].l + ONE));
	
	p=ind[np].p;
	sum.col(p) +=  (ratio * i_lpl[lp_l%4] * T.block(lfc,rfc,2,2) * i_v.col(p));
      }
    }
    
    a(0,n)  =          sum(0,0);
    a(1,n)  =         -sum(0,1);
    b(0,n)  = i_lpl[1]*sum(1,0);
    b(1,n)  = i_lpl[3]*sum(1,1);
    
  }

  return SUCCESS;
}
    
	
    
  
  
   
  
  
