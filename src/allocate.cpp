/**
 *  allocate.cpp -- allocate space for Eigen matrix
 *               -- need to template some of these of course
 *
 **/
#include<allocate.hpp>
#include<iostream>

// REAL, 1 index
bool
allocate(const std::string& desc,
	 std::vector<EMatrix> &M,
	 const size_t n1,
	 const size_t nrows,const size_t ncols,std::ofstream& out)
{
  try {
    M.resize(n1);
    for(size_t i=0;i<n1;i++) {
      M[i].resize(nrows,ncols);
    }
  }
  
  catch(...) { // Get correct specification
    // out << "ERROR: Failure allocating  "<<desc<<std::endl; 
    // out << "nrows = "<<nrows <<" ;   ncols = "<<ncols << std::endl;
    std::cout << "Allocation failed..."<<std::endl;
    return FAIL;
  }
  return SUCCESS;
}

// REAL, 3 index
bool
allocate(const std::string& desc,
	 std::vector<std::vector<std::vector<EMatrix> > > &M,
	 const size_t n1,
	 const size_t n2,
	 const size_t n3,
	 const size_t nrows,const size_t ncols,std::ofstream& out)
{
  try {
    M.resize(n1);
    for(size_t i=0;i<n1;i++) {
      M[i].resize(n2);
      for(size_t j=0;j<n2;j++) {
	M[i][j].resize(n3);
	for(size_t k=0;k<n3;k++) {
	  M[i][j][k].resize(nrows,ncols);
	}
      }
    }
  }
  catch(...) { // Get correct specification
    // out << "ERROR: Failure allocating  "<<desc<<std::endl; 
    // out << "nrows = "<<nrows <<" ;   ncols = "<<ncols << std::endl;
    std::cout << "Allocation failed..."<<std::endl;
    return FAIL;
  }
  return SUCCESS;
}




// COMPLEX, 1 index
bool
allocate(const std::string& desc,
	 std::vector<ECMatrix>  &M,
	 const size_t n1,
	 const size_t nrows,const size_t ncols,std::ofstream& out)
{ 
  try {
    M.resize(n1);
    for(size_t i=0;i<n1;i++) {
      M[i].resize(nrows,ncols);
    }
  }
  catch(...) { // Get correct specification
    // out << "ERROR: Failure allocating  "<<desc<<std::endl; 
    // out << "nrows = "<<nrows <<" ;   ncols = "<<ncols << std::endl;
    std::cout << "Allocation failed..."<<std::endl;
    return FAIL;
  }
  return SUCCESS;
}



// RAW, Fundamental, real_type
bool
allocate(const std::string& desc,
	 EMatrix &M,
	 const size_t nrows,const size_t ncols,std::ofstream& out)
{ 
  try {
    M.resize(nrows,ncols);
  }
  catch(...) { // Get correct specification
    // out << "ERROR: Failure allocating  "<<desc<<std::endl; 
    // out << "nrows = "<<nrows <<" ;   ncols = "<<ncols << std::endl;
    std::cout << "Allocation failed..."<<std::endl;
    return FAIL;
  }
  return SUCCESS;
}





// RAW, Fundamental, complex
bool
allocate(const std::string& desc,
	 ECMatrix &M,
	 const size_t nrows,const size_t ncols,std::ofstream& out)
{ 
  try {
    M.resize(nrows,ncols);
  }
  catch(...) { // Get correct specification
    // out << "ERROR: Failure allocating  "<<desc<<std::endl; 
    // out << "nrows = "<<nrows <<" ;   ncols = "<<ncols << std::endl;
    std::cout << "Allocation failed..."<<std::endl;
    return FAIL;
  }
  return SUCCESS;
}
