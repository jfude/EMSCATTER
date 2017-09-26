#include<computeC0.hpp>
bool 
computeC0(ECMatrix& C0, std::vector<std::vector<std::vector<EMatrix> > > &Yv,
	  const std::vector<AngularMultIndex> &ind,const REAL_TYPE zc,
	  const std::vector<REAL_TYPE> &phi,  const std::vector<REAL_TYPE> &w_phi,
	  const std::vector<REAL_TYPE> &theta,const std::vector<REAL_TYPE> &w_theta,
	  const REAL_TYPE r,const REAL_TYPE w_r,const REAL_TYPE k,std::ofstream& outStream)
{
  
  //  EMatrix  tmp(3,2);
  //  ECMatrix tmpc(3,2);
  const REAL_TYPE PI     = (REAL_TYPE)3.1415926535897932384626433832795028842L;
  const REAL_TYPE ZERO   = (REAL_TYPE)0.000000000000000000L;
  const REAL_TYPE ONE    = (REAL_TYPE)1.000000000000000000L;
  const REAL_TYPE P5     = (REAL_TYPE)0.500000000000000000L;
  const REAL_TYPE TWO    = (REAL_TYPE)2.000000000000000000L;
  const REAL_TYPE FOURPI = (REAL_TYPE)TWO*TWO*PI;
  
  
  size_t nsize = ind.size();
  
  EMatrix Y(3,3);Y.setZero(3,3);
  ECMatrix Z(3,3);Z.setZero(3,3);
  ECMatrix Yt(3,3);Yt.setZero(3,3);
  
  REAL_TYPE sgn[2];
  sgn[0]=ONE;sgn[1]=-ONE;
  const REAL_TYPE k2 = k*k;
  REAL_TYPE csPhi[2],w,sll_1,x,d,plm_dery,m_plm_sin,plm,m_plm;
  std::complex<REAL_TYPE> u,epsilon;
  std::complex<REAL_TYPE> wr2;
  
  size_t nPhi   = w_phi.size();
  size_t nTheta = theta.size();
  size_t i_theta,i_phi;
  int m,n,l,p,np;


  // Precompute w_theta * sin(theta) (*use algo)
  std::vector<REAL_TYPE> wrr_sin(nTheta);
  for(i_theta=0;i_theta<nTheta;i_theta++) wrr_sin[i_theta] = w_r*r*r*w_theta[i_theta] * sin(theta[i_theta]);
    
  // Zero C0
  C0.setZero();
  
  for(i_phi=0;i_phi<nPhi;i_phi++) {  // Phi loop for angular integration
    for(i_theta=0;i_theta<nTheta;i_theta++) { // Theta loop for angular integration
      
      // Loop to generat Yn(theta,phi), for a single value of theta,phi
      //NEED TO CHANGE THIS -- SHOULD COMPUTE WITHOUT RADIAL FACTORS
      // THEN SCALE WITH RADIAL FACTORS
      wr2 = std::complex<REAL_TYPE>(w_phi[i_phi]*wrr_sin[i_theta],ZERO);
      
      //Compute dielectric
      epsilon = dielectric(zc,phi[i_phi],theta[i_theta],r);
      //std::cout << "eps = " << epsilon << std::endl;
      //epsilon = ONE;
      u       = k2*(epsilon - ONE);
      Z(0,0)  = u/epsilon;
      Z(1,1)  = u;
      Z(2,2)  = u;
      
      // multiply wr2 should be moved outside of loop
      // remove the 3* multiply in the C0 indices
      
      for(n=0;n<nsize;n++) {
	Y = Yv[i_phi][i_theta][n];
	
	Yt.noalias() = wr2 * Y.transpose() * Z;
	for(np=0;np<nsize;np++) {
	  C0.block(3*n,3*np,3,3).noalias() += (Yt * Yv[i_phi][i_theta][np]);
	}
      }
      
    }
    
  }
  
  return SUCCESS;
  
}
