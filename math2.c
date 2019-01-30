#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constant.h"

double get_ghs(double n2, double n3, double v2, double d){
	double zeta;
	zeta = 1. - (v2*v2/n2/n2);
	return 1./(1.-n3) + (1./4.)*(zeta*n2)/((1.-n3)*(1.-n3)) + (1./72.)*((zeta*n2*n2)/((1.-n3)*(1.-n3)*(1.-n3)));
}

double get_ghsp2(double n2,double n3,double v2,double d){
	double zeta;
	zeta = 1. - (v2*v2/n2/n2);
	return (1./4.)*1./((1-n3)*(1-n3))*(1.+v2*v2/n2/n2) + (1./36.)*(n2/((1.-n3)*(1.-n3)*(1.-n3)));
}

double get_ghsp3(double n2,double n3, double v2,double d){
	double zeta;
	zeta = 1. - (v2*v2/n2/n2);
	return 1./((1.-n3)*(1.-n3)) + (1./2.)*(zeta*n2)/((1-n3)*(1-n3)*(1-n3)) + (1./24.)*((zeta*n2*n2)/((1-n3)*(1-n3)*(1-n3)*(1-n3)));
}

double get_ghspv2(double n2,double n3,double v2,double d){
	double zeta;
	zeta = 1 - (v2*v2/n2/n2);
	return -(1./2.)*(v2)/(n2*(1-n3)*(1-n3)) - (1./36.)*((v2)/((1-n3)*(1-n3)*(1-n3)));
}

double trapz(double *func, double rMin, double rMax, int N) {
	unsigned int i;
	double o = 0.;
	double h;
	h = (rMax-rMin)/(N-1);

	if (N==1) {
		o = 0.;
	} else {
		for (i=0; i<N; i++) {
			o += func[i]*h;
		}
		o -= 7./12.*h*func[0];
		o += 1./12.*h*func[1]; 
		o -= 7./12.*h*func[N-1];
		o += 1./12.*h*func[N-2];
	}

	return o;
}

void derivative(double *func, double *deriFunc, double rMin, double rMax, double dr) {
	unsigned int i;
	int rMinIndex, rMaxIndex;
	int N;

	rMinIndex = round(rMin/dr);
	rMaxIndex = round(rMax/dr);
	N = (rMaxIndex-rMinIndex)+1;

	if (N==1) {
		printf("We can't define any derivatves in this function\n");
		exit(1);
	} else {
		for (i=0; i<N; i++) {
			if (i==0) {
				deriFunc[rMinIndex+i] = (func[rMinIndex+1]-func[rMinIndex+0])/dr;
			} else if (i==N-1) {
				deriFunc[rMinIndex+i] = (func[rMinIndex+i]-func[rMinIndex+i-1])/dr;
			} else {
				deriFunc[rMinIndex+i] = (func[rMinIndex+i+1]-func[rMinIndex+i-1])/dr;
			}
		}
	}
}

