PROG  = emscatter

OBJS  = scatter.o \
	createAngIndices.o \
	createRadialGrid.o \
	allocate.o \
	computeBessel.o \
	computeHankel.o \
	computeGreen.o  \
	getTrapGrid.o   \
	getGaussWA.o \
	getGaussGrid.o  \
	plegendre.o \
	dielectric.o \
	initializeT.o \
	computeY.o   \
	computeC0.o \
	computeC1.o \
	computeQ12.o \
	computeT.o \
	computeD0.o \
	computeD1.o \
	compute_ab.o \
	computeScatIntFunctions.o \
	saveTtoFile.o \
	restartTfromFile.o \
	computeMatrixMax.o \
	readInputParams.o \
 



# include directories
BOOST_INCDIR_        =   -I../boost/
EIGEN_INCDIR_        =   -I../Eigen/
EMSC_INCDIR_         =   -I../include/
ALLINCDIR_      =    $(EMSC_INCDIR_)  $(BOOST_INCDIR_) $(EIGEN_INCDIR_)  

# lib directories and libs
ALLLIB_  =  


# can add -O optimization here
#CPPLNKFLAGS    =  -LANG:std -64 -mips4 -OPT:Olimit=0
#CPPCMPFLAGS    =  -64 -mips4

CPP            = /usr/bin/g++
CPPLNKFLAGS    =  
#CPPCMPFLAGS    = -c -O3 -DBESSEL -DHANKEL -DGREEN -DC0_MATRIX 
CPPCMPFLAGS    = -c -O
all: $(PROG)
$(PROG): $(OBJS)
	$(CPP) $(CPPLNKFLAGS) -o $@ $(OBJS) $(ALLLIB_)  
	mv $(PROG) ../bin/

.SUFFIXES: 

.SUFFIXES: .o .cpp 


.cpp.o:  
	$(CPP)  $(CPPCMPFLAGS) $(ALLINCDIR_) $< 

clean:
	-\rm -rf $(PROG) $(OBJS) 
cleanall:
	make clean 
