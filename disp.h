#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "math2.h"

double uAtt(double r, double sig, double T, double rc);

void c1DispSpherical(double *c1DispVec, double *rVec, double *denVec, double R,  double T, double rc); 
void c1DispSphericalSolute(double *c1DispVec, double *rVec, double *denVec, double R,double RSolute,  double T, double rc); 
void c1ChainDispFlat(double **c1DispMat,double *rVec,double **denMat,double R,double M, double T,double rc);
void FDispSpherical(double *FDispVec, double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc); 
void FDispSphericalSolute(double *FDispVec, double *c1DispVec, double *rVec, double *denVec, double R,double RSolute, double dr, double sig, double T, double rc); 
