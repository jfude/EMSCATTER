#include<saveTtoFile.hpp>

bool
saveTtoFile(REAL_TYPE r, const ECMatrix &T,const std::string& outFileName) {
  
  
  std::ofstream out;
  
  std::string outName = outFileName;
  out.open(outName.c_str(), std::ios_base::binary);
  if(!out) {
    std::cout << "Stream "<<outFileName<< " failed to open..."<<std::endl;
    return FAIL;
  }

  size_t rank = T.rows();

  out<<r;
  out<<rank;
  
  // slow write, maybe not necessary
  for(size_t i=0;i<rank;i++) {
    for(size_t j=0;j<rank;j++) {out<<T(i,j);}
  }
  
  out.close(); // any checking?

  return SUCCESS;

}
  
  
  
  
