README.TXT

E M S C A T T E R  v0.2 -- Initial Release

This is the initial release of this code, "emscatter" 

This program is provided in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program 
makes use of unmodified GPL and (more weakly) publicly licensed software. 


0) PREAMBLE     
-----------

The EMSCATTER code calculates the T-matrices and scattering profiles of a dielectric scatterer 
encountering a monochromatic plane wave. It is a direct implementation of the method detailed in 
the paper 

B.R. Johnson,  "Invariant imbedding T matrix approach to electromagnetic scattering", 
Applied Optics, vol. 27 (1988) p. 4861


The code is written in C++. It has two external dependencies, the Eigen Matrix/ Linear Algebra
package and the popular Boost routines.

http://eigen.tuxfamily.org
http://www.boost.org


The matrix classes/representations and linear solve routines provided by Eigen are used throughout 
emscatter. One nice feature of both the Boost and Eigen  packages (at least, the parts we are using)  
is that they are provided entirely in header (*.h,*.hpp) files. They are  therefore not compiled into 
separate libraries and then linked.  

  
The EMSCATTER tarball directory structure is as follows.
    	

	0.3/
	    README.TXT -- this document
	    .../bin/
	    .../src/
	    .../include/
            .../Eigen/
            .../boost/
	    
	    

1) INSTALLATION
---------------

There is very little configuration for this initial release. The only adjustments that needed are in the
 simple makefile,  

 	0.3/src/Makefile


basically setting the compiler and related flags

	a) C++ compiler :              CPP = /usr/bin/g++ 
	b) C++ compile flags :         CPPCMPFLAGS = -c -O3




Run "make" at the prompt

     $> make   


The code should compile without issue and copy the executable "emscatter" to the directory 
0.3/bin/. 

The default dielectric structure compiled into the code is a sphere (centered at the origin)
of radius 0.5 and constant dielectric 2.0 (see below).

The code uses "long double" (64 bit) for the real floating point type (as default). This is defined in the "basic_types.hpp" header. 



2) RUNNING THE CODE
-------------------

Before the running the code, one needs to define the dielectric structure and to setup an input file.

       a) Dielectric structure: For this release, the dielectric structure needs to be defined at compile time
          and compiled into the code. The routine defining the dielectric structure is "dielectric.cpp". The prototype is
	  
	  std::complex<REAL_TYPE> dielectric(const REAL_TYPE &zc, const REAL_TYPE &phi,const REAL_TYPE &theta,const REAL_TYPE &r)
	 
	  in plain speak
    
          epsilon = dielectric(zshift,phi,theta,r)

          All arguments to the dielectric function are inputs. zshift allows one to shift the origin of the scatterer along the
          z-axis. For this release, please only use 0.0. The shift has been tested, but special instructions are needed for 
          the input file, and those will be provided in the next release. phi, theta, and r are just the spherical coordinates 
	  of the point in the dielectric to be evaluated. The function returns a complex dielectric value at this point. 
	  The function should also return ONE (vacuum dielectric) outside the dielectric structure.

          A number of examples for dielectric functions are provided in the ./Examples/ directory.

	
       b) Input file: The name of the input file for the code is "emsc_input.in". Lines beginning with a '#' (comment line)
          or blank lines are ignored. Each line with valid input should form a (key,value) pair. Eack keyword should be given 
	  a single value. At the moment, the RunTitle  cannot contain spaces. 
	  
	  
	  i) Valid keywords with example value. (Keywords are case sensitive, the ':' is considered part of the keyword).
	  

	 	# INPUT FILE : emsc_input.in
		# The RunTitle is just a string that describes the run.
		RunTitle: Prolate_spheroid__c=3__b=0.2756644477109__a=2*b

                # The shift of the scatterer away from the origin along the z-axis 
		Zshift:   0.0           
		
		# Wavelength of the incoming plane wave 
		Wavelength:   1.0

                # Definition of integration region
		# Radial integration
		# Number of points for integration along radial direction
                Radius_NumPoints:      160
                # Integrate from Radius Start to Radius_End
		# Radius_Start should not be zero!!! (see below)
		# Radius_End should be at or just beyond edge or greatest extent of scatterer
 		Radius_Start:          0.0004
                Radius_End:            0.5513289
		#
		#		
		# Theta  integration, for zshift=0 this would normally be 0 to PI
		Theta_NumPoints:       240
		Theta_Start:           0.000000000000000
		Theta_End:             3.141592653589793
		#
		#
		# Phi  integration,   for zshift=0 this would normally be 0 to 2*PI
		Phi_NumPoints:         38
		Phi_Start:             0.000000000000000
		Phi_End:               6.283185307179586
		#
		# Max_L, the number of spherical harmonics used for the calculation
		Max_L:                 14
		#
		# Symmetry, one of two valid values: 'axial' or 'none' (see below)
		Symmetry:              axial


	  ii) The Radius_Start value cannot be set exactly to zero, since the Hankel functions used are
	      singular at r=0. We start with a value just off of zero, and this works surprisingly well. It can
              be shown that the T matrices --> 0 as r->0 for well behaved dielectrics. It is surprising however that 
              asymptotic expressions do not seem to be needed.

	  iii) The angular momentum components included in a given calculation are given by
	  
	  l: 1<= l <= Max_L
	  m: 0<= m <= l
	  p: 0<= p <= 1, for m => 1 ; p =0 for m=0
	  
	  If Symmetry: is set to 'none', all these components are included, there are 
	  Max_L *(Max_L +2) components. 
	      
	  If Symmetry: is set to 'axial' , then for a plane wave scatterer, the only angular momentum components
 	  that need to be included in the calculation are those for m=1. There are therefore only Max_L components 
	  in this case.    
	    
        


To run the code, copy the input file to the ./bin directory where the executable is located and type

   $> ./emscatter

 

3) OUTPUT
---------

For this release (0.3), when the code is run it prints the RunTitle to stdout, and then a statement like

T cycle ri= n 

where n is the radial integration step 0 <= n < Radial_NumPoints. Finally, it prints a statement at the end 
of the run to let you know it has completed.

All other output is written to the output file "emsc_output.out". The format of the file is as follows:



    emsc_output.out
    =======================================================================================================
    
    Input Parameters -- prints all of the input parameters in the "emsc_input.in" file
    
    T-matrices    ------------------- T_{n,np} , Tmatrices are 2x2 matrices subscripted by 2 
                                                 angular multi-indices
    n     =  l   m   p   
    np    =  lp  mp  pp
    
   ( T11     T12 )
   ( T21     T22 )_{n,np}    
 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
   a & b coefficients --------------------------------
   l,m,p           ax             ay
   
   l,m,p           bx             by

   ---------------------------------------------------  a & coefficient definitions are given in
   							Johnson, Eqs. 108-112

   Scattering Intensity Functions -----------------------
   Theta            I1             I2

   ----------------------------------------------------- These intensity functions are only well defined for and only computed
                                                         for axially symmetric problems.
   
   Completed Run:Normal Exit...
   =======================================================================================================


4) EXAMPLES
-----------

	Currently, only axially symmetric examples are included, although non-symmetric runs have been tested successfully
        against another code that implements the calculation as one large matrix solve. In addition, it was verified 
	that the T matrices possess the proper symmetry for dielectrics with symmetry and  when choosing Symmetry: 'none'.
	Along with sample output, plots of the scattering intensity functions are included and may be compared to those of 
	Johnson. (The agreement with his results is excellent.)	Non-axially symmetric problems take much longer to run, obviously.

	
	./Examples/Johnson/Sphere    -- Sphere of constant dielectric
	./Examples/Johnson/Luneburg  -- Sphere with quadratic luneburg dielectric, eps = 2 - (r/r_{0})^{2}
	./Examples/Johnson/ProlateSpheroid  -- Results for 7 different prolate spheroids 
	./Examples/Cylinder/         -- Results for 2 different cylinders, radius=0.5, semi-height = 0.5,1.5 


	
	


 

 

  