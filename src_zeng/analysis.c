#include <stdlib.h>
#include "constant.h"
#include "program.h"
#include "hard.h"
#include "math2.h"


void pressureSpherical(double *pNVec, double *denOutVec, double R, double d, double dr, double rhob) {
	double *denVec;
	double *deriDenVec;
	double den, deriDen;
	double gContact;
	int nGrids = round(R/dr) + 1;
	int i;

	denVec = (double*)malloc(nGrids*sizeof(double));
	deriDenVec = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		denVec[i] = denOutVec[i];
		deriDenVec[i] = 0;
	}
	
	derivative(denVec,deriDenVec,0,(nGrids-round(d*0.5/dr)-1)*dr,dr);
	/*
	for (i=0; i<nGrids; i++) {
		printf("%f\n", deriDenVec[i]);
	}exit(1);
	*/

	// Calculate g(r,rhob)
	double RSolute = 0.5;
	double RSystem = RSolute+20.; 
	double *rVec;
	double *vec;
	double *iVec;
	double *oVec;
	double *c1HSVec;
	double *tolVec;
	int n = round(RSystem/dr)+1;
	double eta = (PI*rhob)/6;
	double mu = eta*(8-9*eta+3*eta*eta)/(pow(1-eta,3.));
	double q = 0.01;
	double tol = 1E-3;
	int tolFlag;

	rVec = (double*)malloc(n*sizeof(double));
	vec = (double*)malloc(n*sizeof(double));
	iVec = (double*)malloc(n*sizeof(double));
	oVec = (double*)malloc(n*sizeof(double));
	c1HSVec = (double*)malloc(n*sizeof(double));
	tolVec = (double*)malloc(n*sizeof(double));
	for (i=0; i<n; i++) {
		rVec[i] = i*dr;
		vec[i] = 0;
		iVec[i] = 0;
		oVec[i] = 0;
		c1HSVec[i] = 0;
	}

	for (i=round(RSolute/dr)+round(d*0.5/dr); i<n-round(d*0.5/dr); i++) {
			iVec[i] += rhob;
	}

	int again = 0;
	while (again < 1) {
		for (i=0; i<n; i++) {
			vec[i] = iVec[i];
			iVec[i] = 0;
		}

		c1HSSphericalSolute(c1HSVec, rVec, vec, RSystem, RSolute, d, dr, 1);
		for (i=round(RSolute/dr)+round(d/2/dr); i<n-round(d/2/dr); i++) {
					oVec[i] = rhob*exp(mu-c1HSVec[i]);
		}
	/*	
		for (i=0; i<n; i++) {
			printf("%f\n", oVec[i]);
		}exit(1);
		*/
		

		tolFlag = checkTol(tolVec, vec, oVec, n, tol);
		if (tolFlag == 1){
			again = 1;
		}
		else {
			updateDen(vec, iVec, oVec, n, q);
		}

	}
	for (i=0; i<n; i++) {
		oVec[i] /= rhob;
	}
	free(rVec);
	free(vec);
	free(iVec);
	free(oVec);
	free(c1HSVec);
	free(tolVec);

	// We obtain g(r) (oVec), so can obtain contact value
	i = round(RSolute/dr)+round(d/2/dr);
	gContact = oVec[i];
	/*
	for (i=0; i<n; i++) {
		printf("%f\n",oVec[i]);
	}
	printf("%f\n",gContact);
	*/

	// OK, let's calculate pressure
	printf("gContact = %f\n", gContact);
	for (i=0; i<nGrids-round(d*0.5/dr); i++) {
		den = denVec[i];
		deriDen = deriDenVec[i];
//		pNVec[i] = denVec[i] + (2*PI/3*den*den*gContact*pow(d,3)) - (PI/10*deriDen*deriDen*gContact*pow(d,5)) - (PI/12*den*deriDen*gContact*pow(d,4));
		pNVec[i] = denVec[i] + (2*PI/3*den*den*gContact*pow(d,3)) - (PI/10*deriDen*deriDen*gContact*pow(d,5));
	}

	free(denVec);
	free(deriDenVec);
}
