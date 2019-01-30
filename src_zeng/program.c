#include <stdio.h>
#include <math.h>
#include "math2.h"

int checkTol(double *tolVec, double *denVec, double *denOutVec, int nGrids, double tol) {
	int i, tolFlag;
	double r;
	double max;
	for (i=0; i<nGrids; i++) { tolVec[i] = 0; }
	for (i=0; i<nGrids; i++) {
		tolVec[i] = fabs(denOutVec[i]-denVec[i]);
	}

	double *func;
	func = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		r = i*0.01;
		func[i] = 4*3.1415*r*r*tolVec[i];
	}
	free(func);
	max = trapz(func,0,0.01*(nGrids-1),nGrids);

	if (max>tol) {
		tolFlag = 0;
	} else {
		tolFlag = 1;
	}

	printf("Max is %f.\n", max);

	return tolFlag;

}

void updateDen(double *denVec, double *denInVec, double *denOutVec, int nGrids, double q) {
	int i;
	for (i=0; i<nGrids; i++) {
//		denInVec[i] = denVec[i]*(1.0-q) + q*denOutVec[i];
		denInVec[i] = denOutVec[i];
	}
}
