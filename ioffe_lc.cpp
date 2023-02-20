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

  double nu ;

  string s1 = argv[1] ; // the output filename is passed as input from the command line

  double* zcosar = (double*) calloc(Nfft, sizeof(double)) ; 
  double* zsinar = (double*) calloc(Nfft, sizeof(double)) ;

  FILE* zfile = fopen(s1.c_str(), "w") ;
  
  ///*//*/fprintf(zfile, "# Fourier transform of Hu(x) calculated in the model at initial scale \n\n") ;
  ///*//*/fprintf(zfile, "# P.z (nu = ioffe time) \t cos (GeV) \t sin (GeV) \n") ;

  ioffe_lcu(zcosar, zsinar) ;   /* <------------------------- */

  for(int i = 0; i < Nfft; i++){
    nu = 2.*Pi*((double)(i - Nfft/2))/sqrt(Nfft) ;
    if(fabs(nu)<30.)
      /*//*/fprintf(zfile, "%f,%f,%f \n", nu, zcosar[i], zsinar[i]) ;
  }
  
  free((void*)zcosar) ;
  free((void*)zsinar) ;

  fclose(zfile) ;

  return 0 ;
  
}

