#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "math2.h"


double fphi2(double n2,double n3,double nv2,double d, int FHSFlag);
double fphi3(double n2,double n3,double nv2,double d, int FHSFlag);
double fphiv2(double n2,double n3,double nv2,double d, int FHSFlag);

void c1HSFlat(double *c1HSVec, double *rVec, double *denVec, double R, double d, double dr, int FHSFlag);
void c1HSSpherical(double *c1HSVec, double *rVec, double *denVec, double R, double d, double dr, int FHSFlag);
void c1HSSphericalSolute(double *c1HSVec, double *rVec, double *denVec, double R, double RSolute, double d, double dr, int FHSFlag);

double FHSSpherical(double *rVec, double *denVec, double R, double d, double dr, int FHSFlag);
