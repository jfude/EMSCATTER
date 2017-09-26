/**
 *
 * computeU.cpp -- Computes U matrix
 *              -- 
 */
#include<computeU.hpp>
bool
computeU(std::vector<std::vector<std::vector<EMatrix> > > &U, const REAL_TYPE& k,
	 const size_t& nTheta, const size_t& nPhi, const REAL_TYPE& theta0,
	 const std::vector<AngularMultIndex>& ind,
	 const std::vector<REAL_TYPE> &rv, const size_t &lmax, std::ofstream& outStream)
{
  
  const REAL_TYPE ZERO    = (REAL_TYPE)0.0000000000000000;
  const REAL_TYPE ONE     = (REAL_TYPE)1.0000000000000000;
  const REAL_TYPE TWO     = (REAL_TYPE)2.0000000000000000;
  const REAL_TYPE P5      = (REAL_TYPE)0.5000000000000000;
  const REAL_TYPE TWOPI   = (REAL_TYPE)12.5663706143592/TWO;    
  const REAL_TYPE FOURPI  = (REAL_TYPE)12.5663706143592;    
  
  bool ret;
  size_t i_theta,i_phi,i_r,n,np;
  int p,l,m;
 
  const size_t rsize = rv.size();
  //Define the angular grids
  std::vector<REAL_TYPE> theta(nTheta);
  std::vector<REAL_TYPE> w_theta(nTheta);
  std::vector<REAL_TYPE> phi(nPhi);
  std::vector<REAL_TYPE> w_phi(nPhi);
  
  
  
  ret=getGaussGrid(theta,w_theta,nTheta,ZERO,theta0);
  ret=getGaussGrid(phi,  w_phi,nPhi,ZERO,TWOPI);

  
  /*
  for(i_theta=0;i_theta<nTheta;i_theta++) {
    std::cout << theta[i_theta] << "   " << w_theta[i_theta] << std::endl;
  }
  
  for(i_phi=0;i_phi<nPhi;i_phi++) {
    std::cout << phi[i_phi] << "   " << w_phi[i_phi] << std::endl;
  }
  */
  
  //exit(1);
  
  const size_t nsize = ind.size();
  

  std::cout << nsize << std::endl;
  // Define Y
  std::vector<EMatrix> Y(nsize);
  for(n=0;n<nsize;n++) {
    Y[n].resize(3,3); // check
    Y[n].setZero(3,3);
  }
  
  // Zero U
  for(i_r=0;i_r<rsize;i_r++) {
    for(n=0;n<nsize;n++) {
      for(np=0;np<nsize;np++) {
	U[i_r][n][np].setZero(3,3);
      }
    }
  }
  
  EMatrix Z(3,3);Z.setZero(3,3);
  EMatrix Yt(3,3);Yt.setZero(3,3);
  
  REAL_TYPE sgn[2];
  sgn[0]=ONE;sgn[1]=-ONE;
  const REAL_TYPE k2 = k*k;
  REAL_TYPE csPhi[2],w,r,sll_1,x,d,wr2,plm_dery,m_plm_sin,u,plm,m_plm;
  REAL_TYPE epsilon;
  
  for(i_phi=0;i_phi<nPhi;i_phi++) {  // Phi loop for angular integration
    for(i_theta=0;i_theta<nTheta;i_theta++) { // Theta loop for angular integration
      
      x = cos(theta[i_theta]); // argument for the Legendre polynomials
      // w = weight for angular integration, includes sin(theta) from angular Jacobian 
      w = w_phi[i_phi]*w_theta[i_theta]*sin(theta[i_theta]); 
      
      // Loop to generat Yn(theta,phi), for a single value of theta,phi
      for(n=0;n<nsize;n++) {    
	
	p=ind[n].p;   l=ind[n].l;   m=ind[n].m;
	csPhi[0]=cos(m* phi[i_phi]);csPhi[1]=sin(m* phi[i_phi]);
	sll_1 = sqrt((REAL_TYPE)(l*(l+1)));
	
	m_plm_sin = ZERO;
	d  = ONE;
	if(m>0) {
	  /*
	   * m*Plm/sin(theta) represented by recursion 
	   * \frac{m}{1-x^{2}}P_{l}^{m} = -\frac{1}{2} (P_{l-1}^{m+1} + (l+m-1)*(l+m) P_{l-1}^{m-1})
	   *
	   */
	  d = TWO;
	  m_plm_sin =  -P5 * ( boost::math::legendre_p(l-1,m+1,x) 
	  		       + (l+m-1)*(l+m)* boost::math::legendre_p(l-1,m-1,x) );
	}
	
	// derivative fof Plm(cos(theta)) {frac{\partial Plm(cos(theta))}{\partial theta}
	// = -sin(theta) Plm(cos(theta)) = -sqrt(1-x^{2}) \frac{\partia Plm(x)}{\partial x}
	// = -\frac{1}{2} [ (l+m)(l-m+1) P^{m-1}_{l} - P^{m+1}_{l} ]
	plm_dery    =  -P5 * ((l+m)*(l-m+1)*boost::math::legendre_p(l,m-1,x) - boost::math::legendre_p(l,m+1,x));
	
	// 
	Y[n](1,0) =  sgn[p] * m_plm_sin * csPhi[p];
	Y[n](2,0) =  -plm_dery*csPhi[(p+1)%2];
	Y[n](1,1) =   plm_dery*csPhi[p]; // negative of Y(2,0)
	Y[n](2,1) =  sgn[(p+1)%2] * m_plm_sin * csPhi[(p+1)%2];
	Y[n](0,2) =  sll_1 * boost::math::legendre_p(l,m,x) * csPhi[p];
	//std::cout << l << " " << m << "  "<< p << std::endl;
	//if(l==1 && m==0 && p==0) {
	//  std::cout << Y[n](0,2) << std::endl;
	//}
	// Multiply by Klm
	// Y[n] *= ( sqrt((d * (2*l +1))  /(FOURPI*l*(l+1)))  * (sqrt_factorial(l-m)/sqrt_factorial(l+m)));
	Y[n] *= ( sqrt( (d * (2*l +1))  /(FOURPI*l*(l+1)))  * sqrtff(l,m));
      
      } 
      
      // Loop radius
      for(i_r=0;i_r<rsize;i_r++) {
	r=rv[i_r];
	wr2=w*r*r;
	
	//Compute dielectric
	epsilon = dielectric(phi[i_phi],theta[i_theta],r);
	//std::cout << "eps = " << epsilon << std::endl;
	//epsilon = ONE;
	u       = k2*(epsilon - ONE);
	Z(0,0)  = u/epsilon;
	Z(1,1)  = u;
	Z(2,2)  = u;
	
	for(n=0;n<nsize;n++) {  
	  // Form Yn(theta,phi) loop
	  //std::cout << wr2 <<std::endl;
	  //std::cout << Z <<std::endl;
	  //std::cout << Y[n].transpose() <<std::endl;
	  
	  Yt = wr2 * Y[n].transpose() * Z;
	    //std::cout << Yt << std::endl;
	  // For each n,np,r add this (theta,phi) piece for the integration
	  for(np=0;np<nsize;np++) {
	    // if(i_r==0 && n==0 && np==0) {
	    // std::cout << Y[n](0,2)  << std::endl;
	    //  std::cout << Yt(2,0) << std::endl;
	    //  std::cout << Z       << std::endl;
	    // }
	    /*
	    if(i_r==0 && n==0 && np==0) {
	      std::cout << Yt << std::endl;
	      std::cout << Y[n] << std::endl;
	      std::cout << Z       << std::endl;
	    }
	    */
	    
	    U[i_r][n][np] += (Yt * Y[np]);
	    
	    /*
	      if(i_r==0 && n==0 && np==0) {
	    std::cout << U[0][0][0] << std::endl;
	    }
	    */
	  }
	}
	
	//std::cout << U[i_r][n][np] << std::endl;
	
      } // r loop
      
    } // theta loop
  } // phi loop
  
  //std::cout << U[0][0][0] << std::endl; 
  return SUCCESS;
}  
