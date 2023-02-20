#include <stdlib.h>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include <string>
#include <iostream>
#include <fstream>

#include "model_param.h"

#include "pseudo_pdf_lattice.h"

using namespace std;



int main(int argc, char** argv){

  double a_kzpos, aP_kzpos, norm_kzpos, a_kzall, aP_kzall, norm_kzall ;
  
  double kz, kzmin, dkz, kzmax, Pz, Pzmax, p ;

  string vertex (argv[1]) ;      // the vertex can be dipole or point 
  
  string category (argv[2]) ;    // the category is passed as input from the command line: no_regge, regge, reggekz
  
  string pole (argv[3]) ;        // the category is passed as input from the command line: no_regge, regge, reggekz
  
  string filnm (argv[4]) ;       // the datafile

  
  int Npts = 500 ; //no. of points for plotting, no. of kz for each Pz

  
  FILE* kzfile = fopen(filnm.c_str(), "w") ;
  
  // #1 Pz  #2 kz  #3 kz/Pz  #4 fu norm 1 kzpos  #5 fu norm Pz kzpos  #6 fu norm 1 kzall  #7 fu norm Pz kzall
  
  Pzmax = 65. ;

  Pz = 0. ;
  
  while( Pz < Pzmax)
    {

      // calculating the normalization factor

      
      norm_kzpos = Norm_u(vertex, category, pole, Pz, 0.) ;                          /* <---------------------- */ // kz +ve only
      
      norm_kzall = Norm_u(vertex, category, pole, Pz, kzmin_c) ;                   /* <---------------------- */ // kz all
      
   
      kzmax = 3. * Pz ;
      kzmin = -5. * Pz ;
      
      dkz = (kzmax - kzmin)/(static_cast<double>(Npts) - 1.) ;       
   
      kz = kzmin ;

      
      // calculating the kz dependence

      
      for(int i = 0; i < Npts; i++){

	if ( kz < 0. ) 
	  a_kzpos = aP_kzpos = 0. ;
	else{
	  a_kzpos = H_u_kzPz(category, pole, fabs(kz), Pz, norm_kzpos) ;       /* <---------------------- */ // kz +ve only, norm 1
	  aP_kzpos = Pz*a_kzpos ;                                                                                   // norm Pz
	}

	if (category == "noRegge" ){
	  a_kzall = H_u_kzPz(category, pole, kz, Pz, norm_kzall) ;             /* <---------------------- */ // kz all, norm 1
	  aP_kzall = Pz*a_kzall ;                                                                                   // norm Pz
	}
	else
	  a_kzall = aP_kzall = 0. ;
	
	if(Pz == 0.)
	  p = 10000. ;
	else 
	  p = kz/Pz ;

	//**** Printing ***//
	
	fprintf(kzfile, "%e,%e,%e,%e,%e,%e,%e\n", Pz, kz, p, a_kzpos, aP_kzpos, a_kzall, aP_kzall) ;  
	
	/******************/
	
	kz += dkz ;
      }

      
      //**** Printing ***//
       
      fprintf(kzfile, "\n\n\n") ;

      /******************/

      //printf("Pz %f \t j %f\n", Pz, jd) ;

      if(Pz < 2.) // scans Pz in steps of 0.25 till Pz = 8, then Pz = 16, 32, 64 ... till nearest power of 2 smaller than Pzmax 
	Pz += 0.25 ;
      else
	Pz *=  2. ;

      
  };

  
  fclose(kzfile) ;

  
  return 0 ;
  
}


