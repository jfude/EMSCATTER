#include<computeScatIntFunctions.hpp>
bool
computeScatIntFunctions(std::vector<REAL_TYPE>& I1,std::vector<REAL_TYPE>& I2,std::vector<REAL_TYPE>& thetav,
			const std::vector<AngularMultIndex> &ind,const ECMatrix &a, const ECMatrix &b)
{
  

  /* Calculate scattering intensity functions for axially symmetric problems*/  

  const REAL_TYPE P5      = 0.5000000000000000;
  const REAL_TYPE ONE     = 1.0000000000000000;
  const REAL_TYPE TWO     = 2.0000000000000000;
  const REAL_TYPE PI      = (REAL_TYPE)3.1415926535897932384626433832795028842L;
  const REAL_TYPE TWOPI   = TWO*PI;
  const REAL_TYPE FOURPI  = TWO*TWOPI;
  const REAL_TYPE sqrTWO = sqrt(TWO);
  //const size_t    nTheta  = 200;
  //const REAL_TYPE dtheta  = PI/(REAL_TYPE)(thetav.size());
  size_t ntheta;
  REAL_TYPE dtheta;

  ntheta = thetav.size();
  dtheta = PI/(REAL_TYPE)(ntheta);

  
  // Check allocation
  //thetav.resize(nTheta);
  //I1.resize(nTheta);
  //I2.resize(nTheta);
  
  //std::cout <<" dtheta = "<< dtheta << std::endl;

  int p,m,l;
  REAL_TYPE x,klm;  
  REAL_TYPE theta   = 0.000;
  REAL_TYPE radToDeg=180.0/PI;
  REAL_TYPE m_plm_sin,plm_dery;
  size_t nn;
  const size_t nsize = ind.size();

  std::complex<REAL_TYPE> S_theta,S_phi;
  
  for(nn=0,theta=0.0; nn < ntheta;++nn,theta += dtheta) {
      
   
    thetav[nn] = radToDeg*theta;
    x = cos(theta);
    S_theta=std::complex<REAL_TYPE>(0.0,0.0);
    S_phi  =std::complex<REAL_TYPE>(0.0,0.0); 
    
    //std::cout << nn << std::endl;
    
    for(size_t n=0;n<nsize;n++) {
      
      p=ind[n].p;   l=ind[n].l;   m=ind[n].m;

      REAL_TYPE ld = (REAL_TYPE)l;
      REAL_TYPE md = (REAL_TYPE)m;

      REAL_TYPE pre0         = sqrt( (TWO*ld + ONE)/(TWO*ld - ONE));
      REAL_TYPE pre_m1_l_1   = pre0;
      REAL_TYPE pre_m_1_l_1  = pre0;
      REAL_TYPE pre_m_1_l    = ONE;
      REAL_TYPE pre_m1_l     = ONE;

      
      if((l-m) >1) {
	pre_m1_l_1   *=  sqrt((ld-md)*(ld-md-ONE));
      }


      pre_m_1_l      /=  sqrt((ld+md)*(ld-md+ONE));
      if(m==0) {
	pre_m_1_l  = sqrt((ld+ONE)*ld);
      }
      
      if((l-m)>0) {
	pre_m1_l     *=  sqrt((ld+md+ONE)*(ld-md));
      }
      
      if((l+m) >1) {
	pre_m_1_l_1 /= sqrt((ld+md)*(ld+md-ONE));
      }



      if(m==1 && p==0) {
	
	m_plm_sin   =  -P5 * ( pre_m1_l_1*plegendre(l-1,m+1,x)
                               + (ld+md-ONE)*(ld+md)* pre_m_1_l_1* plegendre(l-1,m-1,x) );
	

	
	plm_dery    =  -P5 * ((ld+md)*(ld-md+ONE) * pre_m_1_l  *plegendre(l,m-1,x) 
			      - pre_m1_l*plegendre(l,m+1,x));

	
	
	
	//klm         =  (TWO*ld + ONE)/(ld*(ld+ONE));
	klm = sqrTWO/sqrt((ld*(ld+ONE)));
	
	S_theta += sqrt(TWOPI*(TWO*ld + ONE)) * klm * (a(0,n) * m_plm_sin  + b(0,n) * plm_dery);
	S_phi   -= sqrt(TWOPI*(TWO*ld + ONE)) * klm * (a(0,n) * plm_dery   + b(0,n) * m_plm_sin);
	
      }
    }

    I1[nn]  = (S_phi   * conj(S_phi)).real();
    I2[nn]  = (S_theta * conj(S_theta)).real();
    
  }

  return SUCCESS;

}





