#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double trapz(double *func, double rMin, double rMax, int N) {
	unsigned int i;
	double o = 0.;

	if (N==1) {
		o = 0.;
	} else {
		for (i=0; i<N; i++) {
			o += 2*func[i];
		}
		o -= func[0];
		o -= func[N-1];
		o *= (rMax-rMin)/(2*(N-1));
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

