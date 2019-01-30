#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "math2.h"

double get_D(double n2,double n3,double v2,double K,double T_assoc,double d);
double get_Dp2(double n2,double n3,double v2,double K,double T_assoc,double d);
double get_Dp3(double n2,double n3,double v2,double K,double T_assoc, double d);
double get_Dpv2(double n2,double n3,double v2,double K,double T_assoc,double d);

double get_X(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d);
double get_Xp2(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d);
double get_Xp3(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d);
double get_Xpv2(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d);


double phip_assoc_n2(double n2, double n3, double v2, double Ma, double K, double T_assoc, double d, int FHSFlag);
double phip_assoc_n3(double n2, double n3, double v2, double Ma, double K, double T_assoc, double d, int FHSFlag);
double phip_assoc_v2(double n2, double n3, double v2, double Ma, double K, double T_assoc, double d, int FHSFlag);


void c1AssocSpherical(double *c1AssocVec, double *rVec, double *denVec, double R, double Ma,double K,double T_assoc,double d,double dr, double rc, int FHSFlag, int wallFlag); 
void c1AssocSphericalSolute(double *c1AssocVec, double *rVec, double *denVec, double R, double RSolute,double Ma,double K,double T_assoc,double d,double dr, double rc, int FHSFlag, int wallFlag); 
void FAssocSpherical(double *FAssocVec, double *rVec, double *denVec, double R, double Ma, double K, double T_assoc, double d, double dr, double rc, int FHSFlag);

void FAssocSphericalSolute(double *FAssocVec, double *rVec, double *denVec, double R,double RSolute, double Ma, double K, double T_assoc, double d, double dr, double rc, int FHSFlag); 
