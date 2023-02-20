#ifndef PSEUDO_PDF_LATTICE_H
#define PSEUDO_PDF_LATTICE_H

#include <iostream>


const double kTmax = 10. ;

const double kzmin_c = -2000. ;


//**********//


void   gauss (int npts, int job, double a, double b, double* x, double* w) ;


double k0_u (std::string pole, double kz, double Pz, double kp) ;

double kdotP_u (std::string pole,double kz, double Pz, double kp) ;

double ksq_u (std::string pole, double kz, double Pz, double kp) ;


double H_u_kpkzPz (std::string vertex, std::string category, std::string pole, double kz, double Pz, double kp) ;

double H_u_kzPz (std::string vertex, std::string category, std::string pole, double kz, double Pz, double norm) ;

double Norm_u (std::string vertex, std::string category, std::string pole, double Pz, double kzmin) ;


void   pseudo_u (std::string vertex, std::string category, std::string pole, double* zcosar, double* zsinar, double Pz, double kzmin) ;


double H_lcu_kp (double x, double kp) ;

double H_lcu (double x) ;

void   ioffe_lcu (double* zcosar, double* zsinar) ;



#endif
