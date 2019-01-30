#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "math2.h"



double phip_chain_n2(double n2, double n3, double v2, double M, double d, int FHSFlag);
double phip_chain_n3(double n2, double n3, double v2, double M, double d, int FHSFlag);
double phip_chain_v2(double n2, double n3, double v2, double M, double d, int FHSFlag);



void c1ChainFlat(double *c1ChainVec, double *rVec, double *denVec, double R, double M, double d, double dr, int FHSFlag);
double FChainFlat(double *rVec, double *denVec, double R, double M, double d, double dr, double rc, int FHSFlag); 
