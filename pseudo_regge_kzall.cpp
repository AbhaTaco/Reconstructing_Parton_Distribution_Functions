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

  float Pz, z ;

  string s1 = argv[1] ; // the output filename is passed as input from the command line

  clock_t t1 = clock() ; //for measuring how long it takes to run the whole program
  int NPz = 130 ;

  float* z0cosar = (float*) calloc(Nfft, sizeof(float)) ; 
  float* zcosar = (float*) calloc(Nfft, sizeof(float)) ; 
  float* zsinar = (float*) calloc(Nfft, sizeof(float)) ;

  FILE* zfile = fopen(s1.c_str(), "w") ;

  ///*//*/fprintf(zfile, "#1 Pz.z \t 2 Pz (GeV) \t 3 z (GeV^-1) \t 4 cos (GeV) \t 5 sin (GeV) \t 6 cos/cos0 \t 7 sin/cos0 \t 8 Pz*cos (GeV^2) \t 9 Pz*sin (GeV^2) \t 10 Pz*cos/cos0 (GeV) \t 11 Pz*sin/cos0 (GeV)\n") ;

  pseudo_reg_u(zcosar, zsinar, 0., kzmin_c) ;          /* <------------------------- */
  
  for(int i = 0; i < Nfft; i++)
    z0cosar[i] = zcosar[i] ;

  for(int j = 1; j < NPz; j++){

    Pz = 0.25*(float)j ;
    
    ///*//*/fprintf(zfile, "#index %d \t Pz = %f GeV  \n", j, Pz) ;
  }

  ///*//*/fprintf(zfile, "\n\n") ;
  
  for(int j = 0; j < NPz; j++){

    Pz = 0.25*(float)j ;
    
    ///*//*/fprintf(zfile, "#index %d \t Pz = %f GeV  \n", j, Pz) ;

    pseudo_reg_u(zcosar, zsinar, Pz, kzmin_c) ;             /* <------------------------- */

    for(int i = 0; i < Nfft; i++){
      z = 2.*Pi*((float)(i - Nfft/2))/sqrt(Nfft) ;
      if(fabsf(z)<10.)
	/*//*/fprintf(zfile, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", Pz*z, Pz, z, zcosar[i], zsinar[i], zcosar[i]/z0cosar[i], zsinar[i]/z0cosar[i], Pz*zcosar[i], Pz*zsinar[i], Pz*zcosar[i]/z0cosar[i], Pz*zsinar[i]/z0cosar[i]) ;
	//	/*//*/fprintf(zfile, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", Pz*z, Pz, z, zcosar[i]*2.*Pi, zsinar[i]*2.*Pi, zcosar[i]/z0cosar[i], zsinar[i]/z0cosar[i], Pz*zcosar[i]*2.*Pi, Pz*zsinar[i]*2.*Pi, Pz*zcosar[i]/z0cosar[i], Pz*zsinar[i]/z0cosar[i]) ;
    }
    
    ///*//*/fprintf(zfile, "\n\n\n") ;

    if(j%30 == 0)
      /****/ printf("pseudo routine Pz %f \t index %d\n", Pz, j) ;
  }

  t1 = clock() - t1 ;  //for measuring how long it takes to run the program
  
  /***/ printf("\n Time for whole program %f minutes \n\n", ((float)t1)/(60.*CLOCKS_PER_SEC)) ;


  free((void*)z0cosar) ;
  free((void*)zcosar) ;
  free((void*)zsinar) ;

  fclose(zfile) ;
  
  return 0 ;
  
}

