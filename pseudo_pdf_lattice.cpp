#include <stdlib.h>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "model_param.h"
#include "pseudo_pdf_lattice.h"




/********************************************************************


string drop_one(string in)
{
  int l = 0 ;
  int i = 0 ;
    
  istringstream buf1(in) ;
  istringstream buf2(in) ;
  
  for(string word; buf1 >> word; )
    l++ ;

  string v[l];

  
  for(string word; buf2 >> word; )
    {
      v[i] = word ;
      i++ ;
    }
  
  string s (v[0]) ;
        
  
  for(i = 0 ; i < l - 1 ; i++)     
     
    if(i == l - 2)
      s += v[i] ;
    else
      s += v[i] + " " ;

  return s ;
    
}


********************************************************************


string drop_two(string in)
{
  int l = 0 ;
  int i = 0 ;
    
  istringstream buf1(in) ;
  istringstream buf2(in) ;
  
  for(string word; buf1 >> word; )
    l++ ;

  string v[l];

  
  for(string word; buf2 >> word; )
    {
      v[i] = word ;
      i++ ;
    }
  
  string s (v[0]) ;
        
  
  for(i = 0 ; i < l - 2 ; i++)     
     
    if(i == l - 3)
      s += v[i] ;
    else
      s += v[i] + " " ;

  return s ;
    
}


********************************************************************/


void gauss(int npts, int job, double a, double b, double* x, double* w)
{    
/*     npts     number of points                                       */
/*     job = 0  rescaling uniformly between (a,b)                      */
/*           1  for integral (0,b) with 50% points inside (0, ab/(a+b))*/
/*           2  for integral (a,inf) with 50% inside (a,b+2a)          */
/*     x, w     output grid points and weights.                        */ 

      int     m, i, j; 
      double  t, t1, pp, p1, p2, p3;
      double  pi = 3.1415926535897932385E0;
      double  eps = 3.e-10;			/* limit for accuracy */
      
      m = (npts+1)/2;
      for(i=1; i<=m; i++)
      {  
         t  = cos(pi*(i-0.25)/(npts+0.5));
         t1 = 1;
         while((fabs(t-t1))>=eps)
         { 
            p1 = 1.0;
            p2 = 0.0;
            for(j=1; j<=npts; j++)
            {
               p3 = p2;
               p2 = p1;
               p1 = ((2*j-1)*t*p2-(j-1)*p3)/j;
	    }
            pp = npts*(t*p1-p2)/(t*t-1);
            t1 = t;
            t  = t1 - p1/pp;
         }   
         x[i-1] = -t;
         x[npts-i] = t;
         w[i-1]    = 2.0/((1-t*t)*pp*pp);
         w[npts-i] = w[i-1];
      } 

      if(job==0) 
      {
         for(i=0; i<npts ; i++)
         {
            x[i] = x[i]*(b-a)/2.0+(b+a)/2.0;
            w[i] = w[i]*(b-a)/2.0;
         }
      }
      if(job==1) 
      {
         for(i=0; i<npts; i++)
         {
            x[i] = a*b*(1+x[i]) / (b+a-(b-a)*x[i]);
            w[i] = w[i]*2*a*b*b /((b+a-(b-a)*x[i])*(b+a-(b-a)*x[i]));
         }
      }
      if(job==2) 
      {
         for(i=0; i<npts; i++)
         {
            x[i] = (b*x[i]+b+a+a) / (1-x[i]);
            w[i] =  w[i]*2*(a+b)  /((1-x[i])*(1-x[i]));
         }
      }
}


//********************************************************************//


double k0_u(string pole, double kz, double Pz, double kp)
{

  double a, P0 ;
  
  a = P0 = 0. ;

  P0 = sqrt(M*M + Pz*Pz) ;

  if (pole ==  "diquark")
    
    a = P0 - sqrt((Pz-kz)*(Pz-kz) + kp*kp + Mx*Mx) ;
  
  else if (pole ==  "quark")
    
    a = sqrt(kz*kz + kp*kp + Lambda*Lambda) ;
      
  return a ;
  
}


//********************************************************************//


double kdotP_u(string pole, double kz, double Pz, double kp)
{

  double a, P0 ;
  
  a = P0 = 0. ;

  P0 = sqrt(M*M + Pz*Pz) ;

  a = k0_u(pole, kz, Pz, kp)*P0 - kz*Pz ;

  return a ;

}


//********************************************************************//


double ksq_u(string pole, double kz, double Pz, double kp)
{

  double a, k0 ;
  
  a = k0 = 0. ;

  k0 = k0_u(pole, kz, Pz, kp) ;

  a = k0*k0 - kz*kz - kp*kp ;

  return a ;
  
}


//********************************************************************//


double H_u_kpkzPz(string vertex, string category, string pole, double kz, double Pz, double kp)
{

  double a, D, ksq, kdotP, k0, P0 ;
  
  a = D = ksq = kdotP = k0 = P0 = 0. ; 
  
  P0 = sqrt(M*M + Pz*Pz) ;

  ksq = ksq_u(pole, kz, Pz, kp) ;
  kdotP = kdotP_u(pole, kz, Pz, kp) ;
  k0 = k0_u(pole, kz, Pz, kp) ;


  if( vertex == "dipole")
    {
      if( pole == "diquark")
	{      
	  D = (ksq -Lambda*Lambda)*(ksq -Lamifbda*Lambda) ;
      	  
	  a = 4.*2.*Pi*(2.*(kdotP + m*M)*k0 - (ksq - m*m)*P0)/(2.*fabs(k0 -P0)*D*D) ;
	
	  if (category == "Regge")
	    { 
	      double reg = (double)pow(fabs((k0 + kz)/(P0 + Pz)), -alpha) ;
	      a *= reg ;
	    }
	}
      else if (pole == "quark")
	{
	  double N0, N1, N2, N3, p1, p2, p3, p4 ;
	  
	  D = (M*M + ksq - 2.*kdotP - Mx*Mx) ;  
	  
	  N0 = 2.*(kdotP + m*M)*k0 - (ksq - m*m)*P0 ; //the numerator
	  N1 = 2.*(kdotP + m*M) ;                     //first derivative of the numerator
	  N2 = 2.*P0 ;
	  N3 = 0. ;
        
	  p1 = N3 - 2.*N2*(k0 - P0) ;
	  p2 = -(4.*N2*(k0 - P0) + 6.*N1) ;
	  p3 = 16.*(N1*(k0 - P0) + N0)*(k0 - P0) ;
	  p4 = -48.*N0*(k0 - P0)*(k0 - P0)*(k0 - P0) ;
	  
	  a = 4.*2.*Pi*(p1/D + p2/(D*D) + p3/(D*D*D) +p4/(D*D*D*D))/(2.*fabs(k0)) ;
	}
    }

  return a ;

}


//********************************************************************//



double H_u_kzPz(string vertex, string category, string pole, double kz, double Pz, double norm){

  double a, p, a1, p1, a2, p2 ;
  
  int Ng = 51 ;
  
  a = p = a1 = p1 = a2 = p2 = 0.00000 ;
  
  double* y = (double*) calloc(Ng, sizeof(double)) ; 
  double* w = (double*) calloc(Ng, sizeof(double)) ;

  gauss(Ng, 0, 0., kTmax, y, w) ;

  if(pole == "both")
    {
      for(int i = 0 ; i < Ng ; i++)
	{
	  p1 = H_u_kpkzPz(vertex, category, "quark", kz, Pz, y[i]) ;
	  p2 = H_u_kpkzPz(vertex, category, "diquark", kz, Pz, y[i]) ;
	
	  if(category == "Reggekz")
	    {   
	      double P0 = sqrt(M*M + Pz*Pz) ;
	      double reg = pow(sqrt(2.)*fabs(kz/(P0 +Pz)), -alpha) ;
	    
	      p1 *= reg ;
	      p2 *= reg ;
	    }
	
	  a1 += w[i]*y[i]*p1 ;
	  a2 += w[i]*y[i]*p2 ;
	}
 
      a = a1+a2 ;
    }
  
  else
    
      for(int i = 0 ; i < Ng ; i++)
	{
	  p = H_u_kpkzPz(vertex, category, pole, kz, Pz, y[i]) ;
	  
	  if(category == "Reggekz")
	    {   
	      double P0 = sqrt(M*M + Pz*Pz) ;
	      double reg = pow(sqrt(2.)*fabs(kz/(P0 +Pz)), -alpha) ;
	      
	      p *= reg ;	      
	    }
	
	  a += w[i]*y[i]*p ;	  
	}
 
  
  a *= norm ;
  
  free((void*)y) ;
  free((void*)w) ;
  
  return a ;

}


//********************************************************************//



double Norm_u(string vertex, string category, string pole, double Pz, double kzmin){

  double a, p, b, kzmax ;
  
  int Ng = 101 ;
  
  b = 0.00000 ;
  
  double* y = (double*) calloc(Ng, sizeof(double)) ; 
  double* w = (double*) calloc(Ng, sizeof(double)) ;

  if(Pz < 2.)
    kzmax = 5000. ;
  else
    kzmax = 1000. ;
  
  gauss(Ng, 0, kzmin, kzmax, y, w) ;
  
  for(int i=0; i<Ng; i++){
    p = H_u_kzPz(vertex, category, pole,  y[i], Pz, 1.) ;
    b += w[i]*p ;
  }
  
  free((void*)y) ;
  free((void*)w) ;

  a = 2./b ;
  
  return a ;
  
}


//********************************************************************//


void pseudo_u(string vertex, string category, string pole, double* zcosar, double* zsinar, double Pz, double kzmin){

  double norm, kz, normu ;
  fftw_complex *in_fft, *out_fft ;
  fftw_plan p ;

  normu = Norm_u(vertex, category, pole, Pz, kzmin) ;
  norm = sqrt((double)Nfft) ;
  
  in_fft = (fftw_complex*) fftw_malloc(Nfft * sizeof(fftw_complex)) ;
  out_fft = (fftw_complex*) fftw_malloc(Nfft * sizeof(fftw_complex)) ;
  
  p = fftw_plan_dft_1d(Nfft, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE) ;

  for(int i = 0 ; i < Nfft ; i++){
    kz = (double)(i - Nfft/2)/sqrt((double)Nfft) ;

    if(kz < kzmin)
      in_fft[i][0] = 0. ;     
    else
      in_fft[i][0] = H_u_kzPz(vertex, category, pole, kz, Pz, normu) ;//exp(-(Pz+0.5)*kz*kz) ; 
    
    in_fft[i][1] = 0. ;      
  }
  
  fftw_execute(p) ;

  for(int i = 0 ; i < Nfft ; i++){
    out_fft[i][0] = pow(-1,i)*out_fft[i][0]/norm ;
    out_fft[i][1] = pow(-1,i)*out_fft[i][1]/norm ;
  }

  for(int i = 0 ; i < Nfft ; i++){
    zcosar[i] = out_fft[(i+Nfft/2)%Nfft][0] ;
    zsinar[i] = out_fft[(i+Nfft/2)%Nfft][1] ;
  }

  fftw_destroy_plan(p) ;
 
  fftw_free(in_fft) ;
  fftw_free(out_fft) ;

}


//********************************************************************//

//********************************************************************//


/***** On the lightcone *****/


double H_lcu_kp(double x, double kp){

  double a, p, D, D2, ak, bk, ksq ;
  double Lambda, Mx, alpha, m, Norm, alphap, M_ ;

  a = p = D = D2 = ak = bk = ksq = 0. ;

  Lambda = 1.018 ;
  Mx     = 0.604 ;
  alpha  = 0.210 ;
  m      = 0.420 ;
  Norm   = 2.043 ;

  alphap = 2.448 ;
  M_     = 0.938 ;
  

  ksq = -(kp*kp +x*Mx*Mx)/(1. -x) +x*M_*M_ ; 
  D   = (ksq -Lambda*Lambda)*(ksq -Lambda*Lambda) ;  
   
  a   = 2.*3.14*((m + x*M_)*(m + x*M_) + kp*kp)/(D*D) ;
  a  *= pow(x,-alpha)/(1. - x) ;

  a  *= Norm ;


  return a ;

}


double H_lcu(double x){

  double a, p, eps ;
  
  int Ng = 51 ;
  
  a = 0.00000 ;
  eps = 0.0001 ;
    
  double* y = (double*) calloc(Ng, sizeof(double)) ; 
  double* w = (double*) calloc(Ng, sizeof(double)) ;

  gauss(Ng, 0, 0., kTmax, y, w) ;

  if(x < eps)
    x = eps ;
  
  for(int i=0; i<Ng; i++){
    
    p = H_lcu_kp(x, (double)y[i]) ;
    a += (double)w[i]*y[i]*p ;
  }
 
  free((void*)y) ;
  free((void*)w) ;
  
  return a ;

}


void ioffe_lcu(double* zcosar, double* zsinar){

  double norm, x ;
  fftw_complex *in_fft, *out_fft ;
  fftw_plan p ;
  
  norm = sqrt((double)Nfft) ;
  
  in_fft = (fftw_complex*) fftw_malloc(Nfft * sizeof(fftw_complex)) ;
  out_fft = (fftw_complex*) fftw_malloc(Nfft * sizeof(fftw_complex)) ;
  
  p = fftw_plan_dft_1d(Nfft, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE) ;
    
  for(int i = 0 ; i < Nfft ; i++){
    x = (double)(i - Nfft/2)/sqrt((double)Nfft) ;
    
    if(x > 1.|| x < 0.)
      in_fft[i][0] = 0. ;
    else
      in_fft[i][0] = H_lcu(x) ;
    
    in_fft[i][1] = 0. ;      
  }

  fftw_execute(p) ;

  for(int i = 0 ; i < Nfft ; i++){
    out_fft[i][0] = pow(-1,i)*out_fft[i][0]/norm ;
    out_fft[i][1] = pow(-1,i)*out_fft[i][1]/norm ;
  }

  for(int i = 0 ; i < Nfft ; i++){
    zcosar[i] = (out_fft[(i+Nfft/2)%Nfft][0]) ;
    zsinar[i] = (out_fft[(i+Nfft/2)%Nfft][1]) ;
  }

  fftw_destroy_plan(p) ;
 
  fftw_free(in_fft) ;
  fftw_free(out_fft) ;

}



