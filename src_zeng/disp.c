#include <stdlib.h>
#include <math.h>

#include "constant.h"
#include "math2.h"

double uAtt(double r, double sig, double T, double rc) {
	double rpMin;
	rpMin = pow(2,(1./6.))*sig;
	double u = 0.;
	if (r<rpMin) {
		u = -1.*(1/T);
	}
	else if (r>=rpMin && r<rc) {
		u = 4*(1/T)*(pow(sig/r,12)-pow(sig/r,6));
	}
	else {
		u = 0.;
	}
	return u;
}

void c1DispFlat(double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc) {
	int i, j, k;
	int rpMaxIndex, rpMinIndex, rppMaxIndex, rppMinIndex, Np, Npp, nGrids;
	double r, rp, rpMax, rpMin, rpp, rppMax, rppMin;
	double cAtt, dist;
	double den, denp, u, xi, c;
	double *cAttVec;

	nGrids = round(R/dr)+1;
	cAttVec = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		c1DispVec[i] = 0.;
	}

	
	for (i=round(sig*0.5/dr); i<nGrids-round(sig*0.5/dr); i++) {
		r = rVec[i];

		for (j=0; j<nGrids; j++) {
			cAttVec[j] = 0;
		}
		for (j=round(sig*0.5/dr); j<nGrids-round(sig*0.5/dr); j++) {
			rp = rVec[j]-rVec[nGrids/2];
			denp = denVec[j];
			dist = abs(j-i)*dr;

			if (dist > rc) {
				cAttVec[j] = 0;
			} else if (dist > sig) {
				if (abs(i+j)*dr > rc) {
					cAttVec[j] = pow(sig,6)/pow(rc,4.) - 0.4*pow(sig,12)/pow(rc,10.) - pow(sig,6)/pow(dist,4.) + 0.4*pow(sig,12)/pow(dist,10.);
				} else {
					cAttVec[j] = pow(sig,6)/pow((i+j)*dr,4) - 0.4*pow(sig,12)/pow((i+j)*dr,10.) - pow(sig,6)/pow(dist,4.) + 0.4*pow(sig,12)/pow(dist,10.);
				}
			} else {
				if ((i+j)*dr > rc) {
					cAttVec[j] = pow(sig,6)/pow(rc,4.) - 0.4*pow(sig,12)/pow(rc,10.) - pow(sig,6)/pow(sig,4.) + 0.4*pow(sig,12)/pow(sig,10.);
				} else if ((i+j)*dr > sig) {
					cAttVec[j] = pow(sig,6)/pow((i+j)*dr,4.) - 0.4*pow(sig,12)/pow((i+j)*dr,10.) - pow(sig,6)/pow(sig,4.) + 0.4*pow(sig,12)/pow(sig,10.);
				} else {
					cAttVec[j] = 0;
				}
			}

			cAttVec[j] *= denp;
		}
		c1DispVec[i] = trapz(cAttVec,rVec[0],rVec[nGrids-1],nGrids);
		c1DispVec[i] *= 2*PI*(1/T);
	}
	for (i=0; i<nGrids; i++) {
		printf("%f\n", c1DispVec[i]);
	}
	exit(1);
	free(cAttVec);
}

void c1DispSpherical(double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc) {
	int i, j, k;
	int rpMaxIndex, rpMinIndex, rppMaxIndex, rppMinIndex, Np, Npp, nGrids;
	double r, rp, rpMax, rpp, rppMax, rppMin;
	double cAtt, dist;
	double den, denp, u, xi, c;
	double *cAttVec;
	double rpMin = pow(2,(1./6.))*sig;

	double d = 1;

	nGrids = (int)(R/dr)+1;


	// Why error?
	cAttVec = (double*)malloc(nGrids*sizeof(double));

	for (i=0; i<nGrids; i++) {
		cAttVec[i] = 0.;
	}

	
	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		for (j=0; j<nGrids; j++) {
			cAttVec[j] = 0;
		}
	//	for (j=1; j<nGrids-round(sig*0.5/dr); j++) {
		for (j=1; j<nGrids; j++) {
			rp = rVec[j];
			denp = denVec[j];
			dist = abs(j-i)*dr;

			if (dist > rc) {
				cAttVec[j] = 0;
			} else if (dist > rpMin) {
				if ((i+j)*dr > rc) {
					cAttVec[j] = pow(sig,6)/pow(rc,4.) - 0.4*pow(sig,12)/pow(rc,10.) - pow(sig,6)/pow(dist,4.) + 0.4*pow(sig,12)/pow(dist,10.);
				} else {
					cAttVec[j] = pow(sig,6)/pow((i+j)*dr,4) - 0.4*pow(sig,12)/pow((i+j)*dr,10.) - pow(sig,6)/pow(dist,4.) + 0.4*pow(sig,12)/pow(dist,10.);
				}
			} else {
				if ((i+j)*dr > rc) {
					cAttVec[j] = pow(sig,6)/pow(rc,4.) - 0.4*pow(sig,12)/pow(rc,10.) - pow(sig,6)/pow(rpMin,4.) + 0.4*pow(sig,12)/pow(rpMin,10.) -0.5*(rpMin*rpMin)+0.5*(dist*dist);
				} else if ((i+j)*dr > rpMin) {
					cAttVec[j] = pow(sig,6)/pow((i+j)*dr,4.) - 0.4*pow(sig,12)/pow((i+j)*dr,10.) - pow(sig,6)/pow(rpMin,4.) + 0.4*pow(sig,12)/pow(rpMin,10.) -0.5*(rpMin*rpMin)+0.5*(dist*dist);
				} else {
					cAttVec[j] = -0.5*(i+j)*dr*(i+j)*dr + 0.5*dist*dist;
				}
			}

			cAttVec[j] *= rp*denp;
		}
		c1DispVec[i] = trapz(cAttVec,rVec[0],rVec[nGrids-1],nGrids);
		c1DispVec[i] *= (2*PI/r)*(1/T);
	}
	c1DispVec[0] = c1DispVec[1];
	/*
	for (i=0; i<nGrids; i++) {
		printf("%f\n", c1DispVec[i]);
	}
	exit(1);
	*/
	free(cAttVec);
}


/**************************************************************/

double FDispSpherical(double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc) {
	int i;
	int nGrids = round(R/dr)+1;
	double r, den;
	double F = 0;
	double *FFunc;
	
	FFunc = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		FFunc[i] = 0;
	}

	nGrids = round(R/dr)+1;
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		den = denVec[i];
		FFunc[i] = 2*PI*r*r*den*c1DispVec[i];
	}

	F = trapz(FFunc,0,R,nGrids);

	free(FFunc);

	return F;
}

