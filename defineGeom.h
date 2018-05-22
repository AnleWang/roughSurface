#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <cstring>
#include <vector>

using namespace std;

typedef ptrdiff_t Lint;
typedef std::complex<double> Complex;

#ifndef _PI_
  #define _PI_
  const double PI = 4.*atan(1.), TWOPI = 2*PI;
#endif

//global variables
double lengthX, lengthY;
Lint nx, ny;
double dqx, dqy, dqx2, dqy2;
int randSeed;

//variable of rough surface
double lambdaS, lambdaR, Hurst;
double rmsValueSet; //in micros

//variable of step profile
double stepHeight;

//variable of dump out
bool dumpStep, dumpSpecSmooth, dumpSpecStep, dumpACF;

Lint nxH, nyHP1, sizeCompl, sizeReal, nxny;
double dx, dy;
double nqxMax, nqyMax;
double *equilPos, *equilStep; // smooth and step self-affine surface

//inline functions
Lint iqx, iqy, irx, iry;
inline Lint getLinRIndex(Lint ix, Lint iy){return(ix*ny+iy);}
inline void get_irxiry(Lint k) {iry = k/ny; iry = k%ny;}
inline Lint getLinCIndex(Lint ix, Lint iy){return(ix*nyHP1+iy);}
inline void get_iqxiqy(Lint k){iqx=k/nyHP1; iqy= k%nyHP1;}
double getQ(Lint);
double rrand(){return(rand()*1./RAND_MAX);}

//subroutines
void initParams();
void addFractal();
double sqrtSpec(double, double);
void continuumAnalysis(double*, double*,double*);
void dumpFractal();
void dumpSpectra();
void dumpHHACF();
