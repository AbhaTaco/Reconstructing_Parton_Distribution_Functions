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

  double Pz, nu, Pzmax ;

  double cosr_kzpos, sinr_kzpos, cosr_kzall, sinr_kzall ;

  
  int j = 0 ;

  clock_t t1 = clock() ; //for measuring how long it takes to run the whole program
  

  string vertex (argv[1]) ;      // the vertex can be dipole or point 
  
  string category (argv[2]) ;    // the category is passed as input from the command line: no_regge, regge, reggekz
  
  string pole (argv[3]) ;        // the pole can be diquark, quark or both
  
  string filnm (argv[4]) ;       // the datafile


  
  FILE* nufile = fopen(filnm.c_str(), "w") ;

  
  double* cos0_kzpos = (double*) calloc(Nfft, sizeof(double)) ;  // cos, Pz = 0, kz +ve only

  double* cos0_kzall = (double*) calloc(Nfft, sizeof(double)) ;  // cos, Pz = 0, kz all

  double* cos_kzpos = (double*) calloc(Nfft, sizeof(double)) ;   // cos, kz +ve only

  double* cos_kzall = (double*) calloc(Nfft, sizeof(double)) ;   // cos, kz all
 
  double* sin_kzpos = (double*) calloc(Nfft, sizeof(double)) ;   // sin, kz +ve only

  double* sin_kzall = (double*) calloc(Nfft, sizeof(double)) ;   // sin, kz all
 
 
  //  #1 nu = P.z  #2 Pz (GeV)  #3 z  #4 cos (GeV) kzpos  #5 sin (GeV) kzpos  #6 cos/cos0 kzpos  #7 sin/cos0 kzpos  #8 cos (GeV) kzall  #9 sin (GeV) kzall  #10 cos/cos0 kzall  #11 sin/cos0 kzall

  //*  M(z,0) calculated by taking Pz=0  *//
  
  pseudo_u(vertex, category, pole, cos0_kzpos, sin_kzpos, 0., 0.) ;           /* <------------------------- */ // kz +ve only

  if (category == "noRegge")
    pseudo_u(vertex, category, pole, cos0_kzall, sin_kzall, 0., kzmin_c) ;      /* <------------------------- */ // kz all


  
  //**********//

  
  Pzmax = 65. ;

  Pz = 0. ;

  
  //*  calculating and printing M(z,Pz)   *//
  
  while( Pz < Pzmax)
    {

      pseudo_u(vertex, category, pole, cos_kzpos, sin_kzpos, Pz, 0.) ;           /* <------------------------- */ // kz +ve only

      if (category == "noRegge")
	pseudo_u(vertex, category, pole, cos_kzall, sin_kzall, Pz, kzmin_c) ;      /* <------------------------- */ // kz all


      for(int i = 0; i < Nfft; i++)
	{
	  nu = 2.*Pi*((double)(i - Nfft/2))/sqrt(Nfft) ;	  
	  if(fabs(nu)<15.)
	    {
	      
	      cosr_kzpos = cos_kzpos[i]/cos0_kzpos[i] ;   sinr_kzpos = -sin_kzpos[i]/cos0_kzpos[i] ;
	      
	      if (category == "noRegge")
		cosr_kzall = cos_kzall[i]/cos0_kzall[i] ;   sinr_kzall = -sin_kzall[i]/cos0_kzall[i] ;
	      else
		cosr_kzall = sinr_kzall = 0. ;
 	      

	      /** Printing **/
	      
	      fprintf(nufile, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", Pz*nu, Pz, nu, cos_kzpos[i], -sin_kzpos[i],  cos_kzall[i], -sin_kzall[i], cosr_kzpos, sinr_kzpos,  cosr_kzall, sinr_kzall) ;

	      
	      /**************/
	      
	    }
	}

      /** Printing **/
      
      fprintf(nufile, "\n\n\n") ;

      
      /**************/

      /** Print to Terminal **/
      
      if(j%5 == 0)
	printf("pseudo routine Pz %f \t j %d\n", Pz, j) ; 
      

      if(Pz < 2.) // scans Pz in steps of 0.25 till Pz = 8, then Pz = 16, 32, 64 ... till nearest power of 2 smaller than Pzmax 
	Pz += 0.25 ;
      else
	Pz *=  2. ;
      
      
      j++ ;
      
    } ;
 

  t1 = clock() - t1 ;  //for measuring how long it takes to run the program

  
  /** Print to Terminal **/
  
  printf("\n Time for whole program %f minutes \n\n", ((double)t1)/(60.*CLOCKS_PER_SEC)) ;


  free((void*) cos0_kzpos) ;
  free((void*) cos0_kzall) ;
  
  free((void*) cos_kzpos) ;
  free((void*) sin_kzpos) ;
  free((void*) cos_kzall) ;
  free((void*) sin_kzall) ;
 
  fclose(nufile) ;
  
  return 0 ;
  
}

