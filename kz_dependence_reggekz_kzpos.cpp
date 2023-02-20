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

  float a, a2, b, kz, kp, kzmin, dkz, kzmax, norm, Pz ;

  string s1 = argv[1] ; // the output filename is passed as input from the command line

  a = a2 = b = kz = kp = dkz = kzmin = 0. ;

  int NPz = 160 ;
  int Npts = 5000 ; 
 
  FILE* kzfile = fopen(s1.c_str(), "w") ; 
  kzmin = 0. ;
  
  // fprintf(kzfile, "# Pz GeV \t kz GeV \t kz/Pz \t fu(kz, Pz) \t fu'(kz, Pz) fu and fu' normalized to 1 and Pz resp when integrated over all kz,  \n\n") ;

  for(int j = 0; j < NPz; j++){

    Pz = 0.25*(float)j ;
    
    // fprintf(kzfile, "#index %d \t Pz %f GeV  \n", j, Pz) ;

    kzmax = 150. ;

    
    norm = Norm_regkz_u(Pz, 0.) ;                      /* <---------------------- */
    
    dkz = (kzmax - kzmin)/((float)Npts - 1.) ;   
   
    kz = kzmin ;
    
    for(int i = 0; i < Npts; i++){    
      a = H_regkz_u_kzPz(kz, Pz, norm) ;            /* <---------------------- */
      a2 = Pz*a ;
      
      if(j == 0)
	b = 10000. ;
      else 
	b = kz/Pz ;
      
      ///*//*/ fprintf(kzfile, "%f \t %f \t %f \t %f \t %f \n", Pz, kz, b, a, a2) ;
      fprintf(kzfile, "%f,%f,%f,%f,%f\n", Pz, kz, b, a, a2) ;
      kz += dkz ;
    }
    //fprintf(kzfile, "\n\n\n") ;
  }

  fclose(kzfile) ;
    
  return 0 ;
  
}
