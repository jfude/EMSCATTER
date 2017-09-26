/*
 * readInputParams.cpp : routine for reading simple keyworded input
 *                       using C style approach which is very very old-fasioned,
 *                       could replace with Boost Tokenizer. In fact, I was using 
 *                       the tokenizer, but since this  is the only place where it was used
 *                       I dropped it so that the program would have one less dependency.
 */

#include<readInputParams.hpp>

bool checkReal(char *word) {

  char *p;
  p=word;
  
  int deccnt=0;
  while(p[0] != NULL_CHAR) {
    if(!isdigit(p[0])) {
      if(p[0] != '.') {
        return 0;
      }
      ++deccnt;
    }
    ++p;
  }

  if(deccnt<2) {return true;}
  return false;  
}


bool checkInt(char *word) {
  
  char *p;
  p=word;
  
  while(p[0] != NULL_CHAR) {
    if(!isdigit(p[0])) {
      return false;
    }
    ++p;
  }
  return true;  
}


void writeTypeCheckError(const std::string& key,char *word,std::ofstream& out)
{
 
  std::ostringstream tmps;
  tmps.str("");
  tmps << word;
  
  out << "ERROR: Input for the keyword " << key << "   " << tmps.str() << std::endl
      << " does not appear to be of the correct type. " << std::endl;
  
  std::cout << "ERROR: Input for the keyword " << key << "   " << tmps.str() << std::endl
	    << " does not appear to be of the correct type. " << std::endl;
  
}


// parse line into individual words, which are separated by spaces
// Return values:
//                0      --  comment line or empty line
//                -1     --  found end of file
//                n (>0) --  number of words found on line
int parseLine(FILE **fp, unsigned char *line,int *lineStart,int *lineLength) {

  int i=0;
  int n=0;

  while( ((line[i++]=fgetc(*fp)) != '\n') && feof(*fp)==0);
  if(feof(*fp)) return -1;
  line[--i]=NULL_CHAR;

  i=0;
  while(1) {
    while(isspace(line[i])) {i++;}
    if(n==0 && line[i]=='#') {
      return 0; //comment
    }
    
    if(line[i]==NULL_CHAR)   return n;
    
    lineStart[n]=i;
    while( (line[++i] != NULL_CHAR) && !isspace(line[i]) );
    lineLength[n]=i - lineStart[n];
    //printf("wordLength = %d\n",lineLength[n]);
    ++n;
  }
}


//read input parameters
bool readInputParams(std::string& runTitle,REAL_TYPE& zshift, REAL_TYPE& wavelength, 
		     coordParam& radius, coordParam& theta, coordParam& phi,
		     size_t & maxL, std::string& symmetry, const std::string& inputFileName, 
		     std::ofstream& out)
{
  
  

  FILE *fp;
  std::ostringstream tmps; 

  
  fp=fopen(inputFileName.c_str(),"r");
  if(!fp) {
    std::cout << "Failed to open the input file " << inputFileName << std::endl;
    out       << "Failed to open the input file " << inputFileName << std::endl;
    return false;
  }

  // initialize isSet map
  std::map<std::string,bool>  isSet;
  isSet["RunTitle:"]           =false;
  isSet["Zshift:"]             =false;
  isSet["Wavelength:"]         =false;
  isSet["Radius_NumPoints:"]   =false;
  isSet["Radius_Start:"]       =false;
  isSet["Radius_End:"]         =false;
  isSet["Theta_NumPoints:"]    =false;
  isSet["Theta_Start:"]        =false;
  isSet["Theta_End:"]          =false;
  isSet["Phi_NumPoints:"]      =false;
  isSet["Phi_Start:"]          =false;
  isSet["Phi_End:"]            =false;
  isSet["Max_L:"]              =false;
  isSet["Symmetry:"]           =false;  
  
    
  

  // start read, match to keywords
  int r;
  int wordStart[6],wordLength[6];
  char word1[80],word2[80];
  unsigned char line[100];

  bool readError=false;bool endOfFile=false;
  while(!readError && !endOfFile) {
    
    // Read lines until we find something other than empty or comment line
    // Line should only have 2 elements, key and value
    while((r=parseLine(&fp,&line[0],&wordStart[0],&wordLength[0]))==0);
    
    // found keyword and value
    if(r==2) {
      strncpy(word1,(const char *)&line[wordStart[0]],wordLength[0]);
      word1[wordLength[0]]=NULL_CHAR;
      strncpy(word2,(const char *)&line[wordStart[1]],wordLength[1]);
      word2[wordLength[1]]=NULL_CHAR;

      
      /**********  Real Types **********************************/
      if(!strcmp(word1,"Zshift:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Zshift:",word2,out);
	  readError=true;
	  continue;
	}
	
   
	zshift = (REAL_TYPE)strtold(word2,NULL);
	isSet["Zshift:"]=true;
	continue;
      }
      
      
      // wavelength
      if(!strcmp(word1,"Wavelength:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Wavelength:",word2,out);
	  readError=true;
	  continue;
	}

	wavelength = (REAL_TYPE)strtold(word2,NULL);
	isSet["Wavelength:"]=true;
	continue;
      }
      

      // radius integration start
      if(!strcmp(word1,"Radius_Start:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Radius_Start:",word2,out);
	  readError=true;
	  continue;
	}	

	radius.start = (REAL_TYPE)strtold(word2,NULL);
	isSet["Radius_Start:"]=true;
	continue;
      }

      // radius integration end
      if(!strcmp(word1,"Radius_End:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Radius_End:",word2,out);
	  readError=true;
	  continue;
	}
	
	radius.end = (REAL_TYPE)strtold(word2,NULL);
	isSet["Radius_End:"]=true;
	continue;
      }

      
      // theta start
      if(!strcmp(word1,"Theta_Start:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Theta_Start:",word2,out);
	  readError=true;
	  continue;
	}
	
	theta.start = strtod(word2,NULL);
	isSet["Theta_Start:"]=true;
	continue;
      }


      // theta end
      if(!strcmp(word1,"Theta_End:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Theta_End:",word2,out);
	  readError=true;
	  continue;
	}
	
	theta.end = strtod(word2,NULL);
	isSet["Theta_End:"]=true;
	continue;
      }


      // phi start
      if(!strcmp(word1,"Phi_Start:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Phi_Start:",word2,out);
	  readError=true;
	  continue;
	}
	
	phi.start = strtod(word2,NULL);
	isSet["Phi_Start:"]=true;
	continue;
      }


      // Phi End
      if(!strcmp(word1,"Phi_End:")) {
	if(!checkReal(word2)) {
	  writeTypeCheckError("Phi_End:",word2,out);
	  readError=true;
	  continue;
	}
	
	phi.end = strtod(word2,NULL);
	isSet["Phi_End:"]=true;
	continue;
      }
      
      /***************************************************************************/

      /************** Int Types **************************************************/

      // radius integration num points
      if(!strcmp(word1,"Radius_NumPoints:")) {
	if(!checkInt(word2)) {
	  writeTypeCheckError("Radius_NumPoints:",word2,out);
	  readError=true;
	  continue;
	}
	
	radius.nPoints = (size_t)strtoul(word2,NULL,0);
	isSet["Radius_NumPoints:"]=true;
	continue;
      }

      // theta integration num points
      if(!strcmp(word1,"Theta_NumPoints:")) {
	if(!checkInt(word2)) {
	  writeTypeCheckError("Theta_NumPoints:",word2,out);
	  readError=true;
	  continue;
	}
	
	theta.nPoints = (size_t)strtoul(word2,NULL,0);
	isSet["Theta_NumPoints:"]=true;
	continue;
      }

      // phi integration num points
      if(!strcmp(word1,"Phi_NumPoints:")) {
	if(!checkInt(word2)) {
	  writeTypeCheckError("Phi_NumPoints:",word2,out);
	  readError=true;
	  continue;
	}
	
	phi.nPoints = (size_t)strtoul(word2,NULL,0);
	isSet["Phi_NumPoints:"]=true;
	continue;
      }
      
      // number of ang. moment components
      if(!strcmp(word1,"Max_L:")) {
	if(!checkInt(word2)) {
	  writeTypeCheckError("Max_L:",word2,out);
	  readError=true;
	  continue;
	}
	
	maxL = (size_t)strtoul(word2,NULL,0);
	isSet["Max_L:"]=true;
	continue;
      }

      /********************************************************************/

      /************************ String Types ******************************/
      
      // Run Title
      if(!strcmp(word1,"RunTitle:")) {
	
	tmps.str(""); //clear stringstream
	tmps << word2; // word2 is null-terminated
	runTitle = tmps.str();
	isSet["RunTitle:"]=true;
	continue;
      }
      



      // Symmetry
      if(!strcmp(word1,"Symmetry:")) {

	tmps.str("");
	tmps << word2; 
	symmetry = tmps.str();
	
	isSet["Symmetry:"]=true;
	continue;
      }
      
      /************************************************************************/

    }
    
    
    // Found end of file
    else if(r==-1) {
      endOfFile=true;
    }
    
    else {  
      std::cout << "Number of words found on line was not 2" << std::endl;
      readError=true;
    }

  }
  
  // clear keyTypes, keyValues
  // close file
  fclose(fp); 
  if(readError) return false;
  
  // check if all parameters are set
 
  for(std::map<std::string,bool>::const_iterator s= isSet.begin(); s != isSet.end(); s++) {
    
    if(s->second==false) {
      std::cout << " ERROR: Keyword " << s->first << "  was not set..." << std::endl;
      return false;
    }
      
  }
    
  return true;

}
