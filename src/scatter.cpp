/**
 *  scatter.cpp -- EMScatter Eigen Full -- 
 *              -- Scattering code -  uses Eigen matrix package and solves the scattering problem via
 *                 interating the Tmatrix recursion relations. 
 *                 See B.R. Johnson,  "Invariant imbedding T matrix approach to
 *                 electromagnetic scattering", Applied Optics, vol. 27 (1988) p. 4861
 **/
#include<iostream>
#include<stdio.h>
#include<vector>
#include<fstream>
#include<complex>
#include<math.h>
#include<basic_types.hpp>
#include<coordParam.hpp>
#include<readInputParams.hpp>
#include<AngularMultIndex.hpp>
#include<createAngIndices.hpp>
#include<getTrapGrid.hpp>
#include<allocate.hpp>
#include<initializeT.hpp>
#include<computeBessel.hpp>
#include<computeHankel.hpp>
#include<computeGreen.hpp>
#include<computeY.hpp>
#include<computeC0.hpp>
#include<computeC1.hpp>
#include<computeD0.hpp>
#include<computeD1.hpp>
#include<computeT.hpp>
#include<computeQdd.hpp>
#include<computeQ12.hpp>
#include<computeMatrixMax.hpp>
#include<restartTfromFile.hpp>
#include<saveTtoFile.hpp>
#include<getGaussGrid.hpp>
#include<compute_ab.hpp>
#include<computeScatIntFunctions.hpp>
#include<Eigen/Dense>

int main() 
{

  // Define Numerical Constants
  const REAL_TYPE ZERO    = (REAL_TYPE)0.0000000000000000000000000000000L;
  const REAL_TYPE TWO     = (REAL_TYPE)2.0000000000000000000000000000000L;
  const REAL_TYPE PI      = (REAL_TYPE)3.1415926535897932384626433832795028842L;
  const REAL_TYPE TWOPI   = TWO*PI;
  const REAL_TYPE FOURPI  = TWO*TWOPI;
  
  // Declare the output filename and open
  std::string outFileName="emsc_output.out";
  std::string inFileName ="emsc_input.in";


  // Define most commonly used Matrix type
  typedef Eigen::Matrix<REAL_TYPE,Eigen::Dynamic,Eigen::Dynamic>                EMatrix;
  typedef Eigen::Matrix<std::complex<REAL_TYPE>,Eigen::Dynamic,Eigen::Dynamic> ECMatrix;
  
  /*~~~~~~~~~~  Declare Input Parameters ~~~~~~~~~~~~~~~~~~~~~~*/
  std::string runTitle;       // string descriptor for a specific run 
  REAL_TYPE zshift;           // shift of origin of scattering region along z-axis
  REAL_TYPE wavelength;       // wavelength in incoming plane wave
  coordParam radiusParam;     // info for integration along radial direction
  coordParam thetaParam;      // info for integration along theta direction
  coordParam phiParam;        // info for integration along phi direction
  size_t maxL;                // maximum angular momentum component (L value of spherical harmonic)
  std::string symmetry;       // symmetry of calculation: axial or none
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  
  // Define other
  bool ret;          // return value
  std::ofstream out; // output stream
  
  
  // Open output file
  out.open(outFileName.c_str());
  if(!out) {
    std::cout << "Error:: Unable to open output file... "<< outFileName << std::endl;
    return 0;
  }
  
  
  //Read input file
  ret=readInputParams(runTitle,zshift,wavelength,radiusParam,thetaParam,phiParam,maxL,symmetry, inFileName, out);

  if(!ret) {
    std::cout << "Read of input file :" << inFileName << " failed." << std::endl;
    out       << "Read of input file :" << inFileName << " failed." << std::endl;
    // close output file
    out.close(); // need check?
    return 0;
  }
    

  // Write the input to the output file
  std::cout << runTitle << std::endl;  
  out << runTitle << std::endl << std::endl << std::endl;
  out << "zshift     = " << zshift     << std::endl;
  out << "wavelength = " << wavelength << std::endl;
  out << "radiusParam.nPoints  =  "  << radiusParam.nPoints << std::endl;
  out << "radiusParam.start    =  "  << radiusParam.start   << std::endl;
  out << "radiusParam.end      =  "  << radiusParam.end     << std::endl;
  out << "thetaParam.nPoints   =  "  << thetaParam.nPoints  << std::endl;
  out << "thetaParam.start     =  "  << thetaParam.start    << std::endl;
  out << "thetaParam.end       =  "  << thetaParam.end      << std::endl;
  out << "phiParam.nPoints     =  "  << phiParam.nPoints    << std::endl;
  out << "phiParam.start       =  "  << phiParam.start      << std::endl;
  out << "phiParam.end         =  "  << phiParam.end        << std::endl;
  out << "maxL                 =  "  << maxL                << std::endl;
  out << "symmetry             =  "  << symmetry            << std::endl;
  out << std::endl << std::endl;
  //


  //
  // const bool saveT = false;
  // const std::string saveTFile    = "restart.dat";
  // const bool restartT = false;
  // const std::string restartTFile = "restart.dat"; 
  //
  
  
  /* Define the ang momentum index */
  std::vector<AngularMultIndex> angInd;
  if(!createAngIndices(maxL,angInd,symmetry,out)) return 0;
  
  const size_t nsize  = angInd.size();
  const REAL_TYPE k = TWOPI/wavelength;
  

  /* ~~~~~~ Define and allocate matrix container objects upfront ~~~~~~~~~~~~~~~~~~~~ */

  // Bessel Matrix, diagonal in (l,l') indices for a given r. 
  // Fundamental matrix is 3x2. Index layout: J[l](,)
  std::vector<EMatrix>   J;
  
  // Hankel Matrix, diagonal in (l,l') indices for a given r. 
  // Fundamental matrix is 3x2. Index layout: H[l](,)
  std::vector<ECMatrix>  H;

  // Radial greens function matrix composed from products of
  // J and H. G is diagonal in (l,l'), this is diagonal component of g(r,r')
  // for a given r
  // Fundamental matrix is 3x3, Index layout: G[l](,) 
  std::vector<ECMatrix>  G;
  
  // U matrix: Fundamental matrix is 3x3, Index layout 
  // U[n][n'](,)
  //std::vector<std::vector<EMatrix> > U;
  
  // C0 = w*U at a given r, raw matrix (3*nsize,3*nsize)
  //ECMatrix C0;
  // can maybe get away with EMatrix
  ECMatrix C0;
  
  // C1 = 1 - C0*g, raw matrix (3*nsize,3*nsize)
  ECMatrix C1;

  // Q matrix, raw matrix, (2*nsize,2*nsize) solution of C1*Q=C0;
  ECMatrix Q;
  
  // T matrix, raw matrix (2*nsize,2*nsize)
  ECMatrix T;

  // Q11 = Jt * Q * J  matrix nsize2 x nsize2 
  // Q12 = (Jt * Q * H) + I  matrix nsize2 x nsize2 
  // Q22 = Ht * Q * H  matrix nsize2 x nsize2 
  ECMatrix Q11,Q12,Q22; 

  // D0 = T(r_{i-1}) *[I + Q12.transpose())], nsize2 x nsize2
  ECMatrix D0;

  // D1 = I - T(r_{i-1})*Q22
  ECMatrix D1;

  // q matrix as in D1*q=D0;, (2*nsize, 2*nsize)
  ECMatrix q;

  // check matrices
  ECMatrix QCcheck,QDcheck;
  
  // Y matrices
  std::vector<std::vector<std::vector<EMatrix> > > Y;

  // a and b coeffs
  ECMatrix a_coeff,b_coeff;
  
  size_t i,j,pp;
  
  
  
  /* Compute radial angular grids */
  //should wrap allocations in try
  
  std::vector<REAL_TYPE>   r(radiusParam.nPoints);
  std::vector<REAL_TYPE> w_r(radiusParam.nPoints);
  std::vector<REAL_TYPE>   theta(thetaParam.nPoints);
  std::vector<REAL_TYPE> w_theta(thetaParam.nPoints);
  std::vector<REAL_TYPE>   phi(phiParam.nPoints);
  std::vector<REAL_TYPE> w_phi(phiParam.nPoints);


  if(!getTrapGrid( r    ,w_r    ,radiusParam.nPoints,radiusParam.start  ,radiusParam.end))           return 0;
  if(!getGaussGrid(theta,w_theta,thetaParam.nPoints ,thetaParam.start,thetaParam.end,out))          return 0;
  if(!getGaussGrid(phi  ,w_phi  ,phiParam.nPoints   ,phiParam.start  ,phiParam.end,out))            return 0;
    
  
  // Write grids to output file
  /*
  std::cout << "rsize = " << r.size() << std::endl;
  for(i=0;i<nTheta;i++) {
    std::cout << theta[i] << "   "<< w_theta[i] << std::endl;
  }
  for(i=0;i<nPhi;i++) {
    std::cout << phi[i] << "   "<< w_phi[i] << std::endl;
  }
  */
  
  
  
  // Allocate other matrices used in the T cycle loop
  if(!allocate("J Bessel Matrix",         J,maxL,3,2,out)) return 0;
  if(!allocate("H Hankel Matrix",         H,maxL,3,2,out)) return 0;
  if(!allocate("G Greens Function Matrix",G,maxL,3,3,out)) return 0;
  if(!allocate("C0  Matrix",      C0,3*nsize,3*nsize,out)) return 0;
  if(!allocate("C1  Matrix",C1,3*nsize,3*nsize,out))   return 0;  
  if(!allocate("Q  Matrix  ",Q,  2*nsize,2*nsize,out)) return 0;
  if(!allocate("Q11  Matrix",Q11,2*nsize,2*nsize,out)) return 0;
  if(!allocate("Q12  Matrix",Q12,2*nsize,2*nsize,out)) return 0;
  if(!allocate("Q22  Matrix",Q22,2*nsize,2*nsize,out)) return 0;
  if(!allocate("D0   Matrix",D0 ,2*nsize,2*nsize,out)) return 0;
  if(!allocate("D1   Matrix",D1 ,2*nsize,2*nsize,out)) return 0;  
  if(!allocate("q    Matrix",q,2*nsize,2*nsize,out))   return 0;
  if(!allocate("T    Matrix",T,2*nsize,2*nsize,out))   return 0;
  if(!allocate("a_coeff Matrix",a_coeff,2,nsize,out))  return 0;
  if(!allocate("b_coeff Matrix",b_coeff,2,nsize,out))  return 0;
  // Qcheck matrices
  if(!allocate("QC check Matrix",QCcheck,3*nsize,3*nsize,out))   return 0;
  if(!allocate("QD check Matrix",QDcheck,2*nsize,2*nsize,out))   return 0;
  
  if(!allocate("Y matrix",Y,phiParam.nPoints,thetaParam.nPoints,nsize,3,3,out))  return 0;
 

  /* ~~~~~~~~~~~~~~~~~~~  Loop over radial increments  ~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  // Might be able to combine J and H compute, like compute<sph_bessel>                  

  // Iniitialize T matrix here (for restarting)
  /******************************************************************************
  if(restartT) {
    std::cout << "Restarting T matrix from File"<<std::endl;
    if(!restartTfromFile(r_begin,T,restartTFile)) return 0;
    r_begin += r_begin_shift;
  }
  else {
    // initialize at zero
    std::cout << "Initializing Tmatrix at zero.."<< std::endl;
    
    // reset rbegin here
  }
  *******************************************************************************/
  
  //Initialize T
  if(!initializeT(T)) {
    // check print failure
    out.close();
    return 0;
  }
  
  // Compute Y matrices
  if(!computeY(Y,angInd,maxL,theta,w_theta,phi,w_phi,out)) {
    // check print failure
    out.close();
    return 0;
  }

 
  bool computeError=false;
  REAL_TYPE qMax;
  int jj;
  
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~ T  C y c l e   L o o p ~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for(i=0;(!computeError) && i<radiusParam.nPoints;i++) {
    
    // COMPUTE BESSEL MATRIX
    if(!computeBessel(J,k,r[i],maxL,out))    {
      computeError=true;continue;
    }
    
    // COMPUTE HANKEL MATRIX
    if(!computeHankel(H,k,r[i],maxL,out))    {
      computeError=true;continue;
    }
    
    // COMPUTE GREENS FUNCTION MATRIX
    if(!computeGreen(G,J,H,k,out)) {
      computeError=true;continue;
    }

    // Compute C0 Matrix , C0 = w_{n}*U(r_{n}), C0 is my notation 
    // w*U is from Eq. (84) in Johnson
    if(!computeC0(C0,Y,angInd, zshift, phi,w_phi,theta,w_theta,r[i],w_r[i],k,out)) {
      computeError=true;continue;
    }    
    
    
    // Compute C1 matrix C1 = 1 - C0*g = 1 - w*U*g from Eq. (84)
    if(!computeC1(C1,C0,G,angInd,out)) {
      computeError=true;continue;
    }
    
    
    // Solve Eq. (84)
    // need to catch specific exceptions
    try {
      Q  = C1.partialPivLu().solve(C0);
    }
    catch(...) {
      std::cout << "Problem solving for Q" << std::endl;
      out << "Problem solving for Q" << std::endl;
      computeError=true;continue;
    }
    
    
    // Compute solution difference
    //QCcheck.noalias() = C1*Q - C0;
    //qMax = computeMatrixMax(QCcheck);
    //std::cout << "max_norm(C1*Q - C0)  = " << qMax << std::endl;
    
    
    // Compute Q11, Eq. 85
    if(!computeQdd<EMatrix> (Q11,Q,J,  angInd,k,out)) {
      computeError=true;continue;
    }
    // Compute Q12, Eq. 86
    if(!computeQdd<ECMatrix>(Q22,Q,H,  angInd,k,out)) {
      computeError=true;continue;
    }
    // Compute Q22, Eq. 87
    if(!computeQ12(Q12,Q,J,H,angInd,k,out)) {
      computeError=true;continue;
    }
   
    
    // Compute D0 = T * (1 + Q12)
    if(!computeD0(D0,T,Q12,out)) {computeError=true;continue;} 
    // Compute D1 = 1 - T*Q22
    if(!computeD1(D1,T,Q22,out)) {computeError=true;continue;}
    
    
    // By combining Eqs. 94 & 95, one gets D1*q = D0
    // 
    // Here we solve this equation for q
    // need to catch all the exceptions
    try {
      q  = D1.partialPivLu().solve(D0);
    }
    catch(...) {
	std::cout << "Problem solving for Q..." << std::endl;
	return 0;
    }   
  
    
    //std::cout << "Computed q  Matrix"<< std::endl;
    //QDcheck.noalias() = D1*q - D0;
    //qMax = computeMatrixMax(QDcheck);
    //std::cout << "max_norm(D1*Q - D0)  = " << qMax << std::endl;

    
    
    // Compute next value of T, this is Eq. 97
    T.noalias() = Q11 + q + Q12*q;
    
    std::cout << "T cycle ri = " << i << std::endl;
    
  }

  // need to take down most memory here

  
  if(computeError) {
    out.close();
    return 0;
  }

  // Save T matrix if requested
  // if(saveT) {
  //  if(!saveTtoFile(r[nRadPts-1],T,saveTFile)) return 0;
  //}
  
  
  // Print T matrices to output file

  out << "~~~~~~~~~~~~~~~~~~  T  M a t r i c e s  ~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;   
  for(i=0;i<nsize;i++) {
    size_t lfc = 2*i;
    for(j=0;j<nsize;j++) {
      size_t rfc = 2*j;
      
      out << "n  =   " << angInd[i].l << "   " << angInd[i].m << "   " << angInd[i].p << std::endl;
      out << "np =   " << angInd[j].l << "   " << angInd[j].m << "   " << angInd[j].p << std::endl;
      out << T.block(lfc,rfc,2,2) << std::endl;
      out << std::endl << std::endl;
    }
  }
  out << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;   
  out << std::endl << std::endl;  
  
  // need to flush the output

  
  // Compute a & b coefficients (see Johnson Eqs. 108 - 112)
  if(!compute_ab(a_coeff,b_coeff,T,angInd)) {
    out.close();
    return 0;
  }
    
  
  // Print a & b coefficients
  out << "a & b coefficients --------------------------------" << std::endl;
  out << "n           ax             ay" << std::endl;
  for(i=0;i<nsize;i++) {
    out << angInd[i].l << angInd[i].m << angInd[i].p << "   " << a_coeff(0,i)  << "    " << a_coeff(1,i) << std::endl;
  }
  
  out << "n           bx             by" << std::endl; 
  for(i=0;i<nsize;i++) {
    out << angInd[i].l << angInd[i].m << angInd[i].p << "    "<<  b_coeff(0,i)  << "    " << b_coeff(1,i) << std::endl;
  }
  out << "---------------------------------------------------" << std::endl;
  out << std::endl << std::endl;  
  


  // For axially symmetric problems  
  std::vector<REAL_TYPE> Itheta,I1,I2;
  const size_t nItheta = 1000;
  
  REAL_TYPE p11;
  if(symmetry=="axial") {
    Itheta.resize(nItheta);
    I1.resize(nItheta);I2.resize(nItheta); // should be allocated before T cycle loop, checked
    if(!computeScatIntFunctions(I1,I2,Itheta,angInd,a_coeff,b_coeff)) {
      out.close();
      return 0;
    }
    
    out << std::endl << "Scattering Intensity Functions -----------------------" << std::endl;
    out << "Theta            I1             I2              P11  " << std::endl;
    
    
    for(i=0;i<Itheta.size();i++) {   
      out.width(15);
      out.precision(8);
      p11 = (I1[i] + I2[i])/TWO;
      out << Itheta[i] << "         " << I1[i] << "         " << I2[i] << "        " << p11 << std::endl;
    }
    
  }

  
  out       << "Completed Run:Normal Exit..." << std::endl;
  std::cout << "Completed Run:Normal Exit..." << std::endl;
  
  out.close();
  return 0;
  
}
