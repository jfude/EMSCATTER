#include<computeY.hpp>
bool
computeY(std::vector<std::vector<std::vector<EMatrix> > > &Yv, 
	 std::vector<AngularMultIndex> &ind,const size_t lmax,
	 std::vector<REAL_TYPE> &theta,std::vector<REAL_TYPE> &w_theta,
	 std::vector<REAL_TYPE> &phi,std::vector<REAL_TYPE> &w_phi,
	 std::ofstream& outStream)
{
  
  const REAL_TYPE ZERO = (REAL_TYPE)0.0000000000000000000000L;
  const REAL_TYPE ONE  = (REAL_TYPE)1.0000000000000000000000L;
  const REAL_TYPE P5 = (REAL_TYPE)0.50000000000000000000000L;
  const REAL_TYPE TWO  = (REAL_TYPE)2.0000000000000000000000L;
  const REAL_TYPE PI      = (REAL_TYPE)3.1415926535897932384626433832795028842L;
  const REAL_TYPE FOURPI  = TWO*TWO*PI;  

  const size_t nsize = ind.size();
  EMatrix Y(3,3);Y.setZero(3,3);

  REAL_TYPE sgn[2];
  sgn[0]=ONE;sgn[1]=-ONE;
  //const REAL_TYPE k2 = k*k;
  REAL_TYPE csPhi[2],w,x,d,plm_dery,m_plm_sin,plm,m_plm;

  size_t nPhi   = phi.size();
  size_t nTheta = theta.size();
  size_t i_theta,i_phi;
  int m,n,l,p,np;

  //std::cout << nTheta << "   "<< nPhi << std::endl;

  //std::cout << Yv.size() <<std::endl;
  //std::cout << Yv[0].size() <<std::endl;
  //std::cout << Yv[0][0].size() <<std::endl;
  
  //exit(1);

  // check allocation

  // Should be a better way to organize this
 
  std::vector<REAL_TYPE> sll_1(lmax+1);
  for(n=1;n<=lmax;n++) sll_1[n] = sqrt((REAL_TYPE)(n*(n+1)));  
  
  for(i_phi=0;i_phi<nPhi;i_phi++) {  // Phi loop for angular integration
    
    //std::cout << "PHI" << std::endl;
    //std::cout << i_phi << "           " << phi[i_phi] << std::endl;
    
    for(i_theta=0;i_theta<nTheta;i_theta++) { // Theta loop for angular integration                              
      
      x = cos(theta[i_theta]); // argument for the Legendre polynomials        
      
      // Loop to generat Yn(theta,phi), for a single value of theta,phi                   
      for(n=0;n<nsize;n++) {
	//std::cout << n << std::endl;
	//Yv[i_phi][i_theta][n].setZero();
        p=ind[n].p;   l=ind[n].l;   m=ind[n].m;
	REAL_TYPE ld = (REAL_TYPE)l;
	REAL_TYPE md = (REAL_TYPE)m;
	
	// can reduce further by precalculating
        csPhi[0]=cos(m* phi[i_phi]);csPhi[1]=sin(m* phi[i_phi]);
        //sll_1 = sqrt((REAL_TYPE)(l*(l+1)));

	// Compute normalized legendre prefactors
	REAL_TYPE pre0         = sqrt( (TWO*ld +ONE)/(TWO*ld - ONE));
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


	// The plegendre routine does not handle the case where m<0
	// In the legendre functions below the only case where m<0 is when 
	// m=-1. Use the identity 
	// P_{l}^{-1} = (-1/(l*(l+1)) * P_{l}^{1}
	// See the line with plm_dery
	// Note we do not check the case l=0, since l=0 does not conventionally appear
	// in these calculations

	m_plm_sin = ZERO;
        d  = ONE;
	REAL_TYPE neg_prefact = -ONE/(ld*(ld+ONE));
        int m_1 = 1;
        if(m>0) {
          /*                                                                                                                                         
           * m*Plm/sin(theta) represented by recursion                                                                                               
           * \frac{m}{1-x^{2}}P_{l}^{m} = -\frac{1}{2} (P_{l-1}^{m+1} + (l+m-1)*(l+m) P_{l-1}^{m-1})                                                 
           *                                                                                                                                         
           */
	  m_1 = m-1;
          d = sqrt(TWO); // need to predefine
	  neg_prefact = ONE;
          m_plm_sin =  -P5 * ( pre_m1_l_1*plegendre(l-1,m+1,x)
                               + (ld+md-ONE)*(ld+md)* pre_m_1_l_1 *plegendre(l-1,m-1,x) );
        }

	
        // derivative fof Plm(cos(theta)) {frac{\partial Plm(cos(theta))}{\partial theta}                                                            
        // = = -sqrt(1-x^{2}) \frac{\partia Plm(x)}{\partial x}                                                          
        // = -\frac{1}{2} [ (l+m)(l-m+1) P^{m-1}_{l} - P^{m+1}_{l} ]                                                                                 
        plm_dery    =  -P5 * ((ld+md)*(ld-md+ONE) * pre_m_1_l  * neg_prefact * plegendre(l,m_1,x) - pre_m1_l*plegendre(l,m+1,x));                                                                                                                                           
        Y(1,0) =  sgn[p] * m_plm_sin * csPhi[p];
        Y(2,0) =  -plm_dery*csPhi[(p+1)%2];
        Y(1,1) =   plm_dery*csPhi[p];
        Y(2,1) =  sgn[(p+1)%2] * m_plm_sin * csPhi[(p+1)%2];
        Y(0,2) =  sll_1[l] * plegendre(l,m,x) * csPhi[p];
	
	Y     *= (d/sll_1[l]);

	
	Yv[i_phi][i_theta][n] = Y;
      
      }
      
    }
    //std::cout << " Finished --- theta"<< std::endl; 
  }

  return SUCCESS;
  
}
