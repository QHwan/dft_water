#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "math2.h"

double uAtt(double r, double sig, double T, double rc);

void c1DispFlat(double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc);
void c1DispSpherical(double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc);

double FDispSpherical(double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc);
