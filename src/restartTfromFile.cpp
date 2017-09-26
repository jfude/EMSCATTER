#include<restartTfromFile.hpp>

bool
restartTfromFile(REAL_TYPE& r, ECMatrix &T,const std::string& inFileName) {
  
  
  std::ifstream in;
  size_t rank=0;
  
  std::string inName = inFileName;
  in.open(inName.c_str(),std::ios_base::binary);
  if(!in) {
    std::cout << "Stream "<<inName<< " failed to open..."<<std::endl;
    return FAIL;
  }

  // Do sanity check on matrix size
  
  in>>r;
  std::cout << "r = " << r << std::endl; 
  in>>rank;
  if(rank != T.rows()) {
    std::cout << "restartTFromFile:" << std::endl;
    std::cout << "rankFromFile = "<<rank << " does not match rank of allocated Tmatrix"
	      <<std::endl;
    std::cout << "rank = "<< T.rows() << std::endl;
    return FAIL;
  }
  
  // slow read, maybe not necessary
  for(size_t i=0;i<rank;i++) {
    for(size_t j=0;j<rank;j++) {in>>T(i,j);}
  }
  
  in.close(); // any checking?

  return SUCCESS;
}
  
  
  
  
