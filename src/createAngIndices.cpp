/*
 *
 * createAngIndices.cpp ----
 *
 */
#include<basic_types.hpp>
#include<createAngIndices.hpp>
#include<iostream>
bool
createAngIndices(const size_t& maxL,std::vector<AngularMultIndex>& ind, const std::string& symmetry,
		 std::ofstream& out)
{
    
  size_t lmax;
  if(symmetry=="axial") {
    lmax=maxL;
  }
  else if(symmetry=="none") {
    lmax=maxL*(maxL+2);
  }
  else {
    out << "Error: Value for keyword symmetry " << symmetry << "unrecognized." << std::endl;
    std::cout << "Error: Value for keyword symmetry " << symmetry << "unrecognized." << std::endl;
    return FAIL;
  }
  
  
  try{
    ind.resize(lmax);
  }
  catch(std::bad_alloc const&) {
    out       << "Error: createAngIndices:"<<std::endl;
    out << "Allocation of memory for Angular Multi Index failed."<<std::endl;
    return FAIL;
  }
  

  size_t n,m,p,l;

  if(symmetry=="axial") {
    for(n=0,l=1;l<=maxL;++l,++n) {
      ind[n].l=l;ind[n].m=1;ind[n].p=0;
    }
  }
  

  if(symmetry=="none") {
    for(n=0,l=1;l<=maxL;++l) {
      for(m=0;m<=l;++m) {
	for(p=0;p<2;++p) {
	  if(m==0 && p==1) {continue;}
	  ind[n].l=l;ind[n].m=m;ind[n].p=p;
	  n++;  
	}
      }
    }
  }
  
  return SUCCESS;  
}


