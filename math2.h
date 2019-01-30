#include <stdio.h>
#include <stdlib.h>
double get_ghs(double n2, double n3, double v2, double d);
double get_ghsp2(double n2,double n3,double v2,double d);
double get_ghsp3(double n2,double n3, double v2,double d);
double get_ghspv2(double n2,double n3,double v2,double d);

double trapz(double *func, double rMin, double rMax, int N);
void derivative(double *func, double *deriFunc, double rMin, double rMax, double dr);

