#include <stdlib.h>
#include <math.h>
#include "math2.h"
#include "constant.h"

double FIdealSpherical(double *denVec, double R, double d, double dr, double rhob) {
	int i;
	int nGrids = round(R/dr)+1;
	double r, den;
	double F = 0;
	double *FFunc;
	FFunc = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		FFunc[i] = 0;
	}

	for (i=0; i<nGrids; i++) {
		r = i*dr;
		den = denVec[i];
		if (den < 1E-8) {
			FFunc [i] += 0;
		} else {
			FFunc[i] += 4*PI*r*r*den*(log(den/rhob)-1);
		}
		//printf("%f\t%f\t%f\t%f\n",r,den,rhob, FFunc[i]);
	}

	F = trapz(FFunc, 0, R, nGrids);

	free(FFunc);

	return F;
}
