#include <stdlib.h>
#include <math.h>
#include "math2.h"
#include "constant.h"
#include "memory.h"

double FIdealFlat(double *denVec, double R, double d, double dr) {
	int i;
	int nGrids = round(R/dr)+1;
	double den;
	double F = 0;
	double *FFunc;

	FFunc = dVector(nGrids);

	for (i=0; i<nGrids; i++) {
		den = denVec[i];
		if (den==0) {
			FFunc[i] = 0;
		} else {
			FFunc[i] += den*(log(den)-1);
		}
	}

	F = trapz(FFunc, 0, R, nGrids);

	free_dVector(FFunc);

	return F;
}


void FIdealSpherical(double *FIdealVec, double *denVec, double R, double d, double dr) {
	int i;
	int nGrids = round(R/dr)+1;
	double r, den;
	double F = 0;

	for (i=0; i<nGrids; i++) {
		r = i*dr;
		den = denVec[i];
		FIdealVec[i] = den*(log(den)-1);
	}


}

void FIdealSphericalSolute(double *FIdealVec, double *denVec, double R, double RSolute,double d, double dr) {
	int i;
	int nGrids = (int)((R-RSolute)/dr)+1;
	double r, den;
	double F = 0;

	for (i=0; i<nGrids; i++) {
		r = RSolute+i*dr;
		den = denVec[i];
		if (den==0) {
			FIdealVec[i] = 0;
		} else {
			FIdealVec[i] = den*(log(den)-1);
		}
	}

}

