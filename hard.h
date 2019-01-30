#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "math2.h"


double fphi2(double n2,double n3,double nv2);
double fphi3(double n2,double n3,double nv2);
double fphiv2(double n2,double n3,double nv2);

void c1HSFlat(double *c1HSVec, double *rVec, double *denVec, double R, double rc);
void c1HSSpherical(double *c1HSVec, double *rVec, double *denVec, double R, double rc);
void c1HSSphericalSolute(double *c1HSVec, double *rVec, double *denVec, double R, double RSolute, double d, double rc,int FHSFlag);

double FHSFlat(double *rVec, double *denVec, double R, double d, double rc, int FHSFlag); 
void FHSSpherical(double *FHSVec, double *rVec, double *denVec, double R, double d, double rc, int FHSFlag);
void FHSSphericalSolute(double *FHSVec, double *rVec, double *denVec, double R, double RSolute, double d, double rc, int FHSFlag); 
