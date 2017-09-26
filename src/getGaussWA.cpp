#include<getGaussWA.hpp>
bool
getGaussWA(std::vector<REAL_TYPE>& v, std::vector<REAL_TYPE>& w, const size_t level, std::ofstream& out)
{
  
  const long double Zero = 0.000000000000000000000000000L;
  const long double One  = 1.000000000000000000000000000L;
  const long double P5   = 0.500000000000000000000000000L;
  const long double P25  = 0.250000000000000000000000000L;
  const long double Pi   = 3.14159265358979323846264338327950288L;
  const long double Two  = 2.000000000000000000000000000L;
  const long double Eps1  = 6e-14L;

  // submitted vectors for abscissas and weights
  // check allocation
  
  // long double precision
  std::vector<long double> vl;
  std::vector<long double> wl;
  

  
  try {
    v.resize(level);
    w.resize(level);
    vl.resize(level+1);
    wl.resize(level+1);
  }
  catch(std::bad_alloc const&) {
    out << "Error: Allocation of vectors in gaussian grid generation routine  getGausWA  failed."<<std::endl;
    out << "Allocating two ld vectors of size " << level << " and two of  "<< (level+1) << std::endl; 
    std::cout << "Error: Allocation of vectors in gaussian grid generation routine  getGausWA  failed."<<std::endl;
    std::cout << "Allocating two ld vectors of size " << level << " and two of  "<< (level+1) << std::endl; 
    return FAIL;
  }

  
  int i,j;
  size_t n=level;
  const size_t m = (n+1)/2;
  long double z,z1,p1,p2,p3,pp;
  
  
  for(i=1;i<=m;i++) {
    
    z = cos(Pi * ((long double)(i) - P25)/((long double)(n) + P5));
    /*printf("%25.21Lf\n",z);*/
    do {
      p1=One;
      p2=Zero;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((Two * (long double)(j) - One)*z*p2 - (j - One)*p3)/(long double)(j);
      }
      pp= (long double)(n) * (z*p1-p2)/(z*z-One);
      z1=z;
      z=z1- p1/pp;
      
    } while (std::abs(z-z1) > Eps1);
    //printf("%25.21Lf  %25.21Lf  %25.21Lf\n",z,z1,std::abs(z-z1));
    

    vl[i]    = P5*(One - z);
    vl[n+1-i]= P5*(One + z);
    wl[i]=One/((One -z*z)*pp*pp);
    wl[n+1-i]=wl[i];
    
  }
  
  for (i=1;i<=n;i++) {
    vl[i] -= P5;
    vl[i] *= Two;
    wl[i] *= Two;
    
    v[i-1]=(REAL_TYPE)vl[i];
    w[i-1]=(REAL_TYPE)wl[i];
  }
  
  return SUCCESS;

}

