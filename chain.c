#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "constant.h"
#include "math2.h"


double phip_chain_n2(double n2, double n3, double v2, double M, double d, int FHSFlag) {
	double fphi2 = 0.0;
	double n0, zeta, ghs, ghsp2;
	ghs = get_ghs(n2, n3,v2,d);
	ghsp2 = get_ghsp2(n2,n3,v2,d);
	n0 = n2/PI/d/d;
	zeta = 1-(v2*v2/n2/n2);

	fphi2 = (1/PI)*((1-M)/M)*(log(ghs)*(1+v2*v2/n2/n2) + (n2-v2*v2/n2)*(ghsp2/ghs));
	return fphi2;
}

double phip_chain_n3(double n2, double n3, double v2, double M, double d, int FHSFlag) {
	double fphi3 = 0.0;
	double n0, zeta, ghs, ghsp3;
	ghs = get_ghs(n2, n3,v2,d);
	ghsp3 = get_ghsp3(n2,n3,v2,d);
	n0 = n2/PI/d/d;
	zeta = 1-(v2*v2/n2/n2);

	fphi3 = (1/PI)*((1-M)/M)*(n2-v2*v2/n2)*(ghsp3/ghs);
	return fphi3;
}

double phip_chain_v2(double n2, double n3, double v2, double M, double d, int FHSFlag) {
	double fphiv2 = 0.0;
	double n0, zeta, ghs, ghspv2;
	ghs = get_ghs(n2, n3,v2,d);
	ghspv2 = get_ghspv2(n2,n3,v2,d);
	n0 = n2/PI/d/d;
	zeta = 1-(v2*v2/n2/n2);

	fphiv2 = (1/PI)*((1-M)/M)*(log(ghs)*(-2*v2/n2) + (n2-v2*v2/n2)*(ghspv2/ghs));

	return fphiv2;
}

/*********************************************/

void c1ChainFlat(double *c1ChainVec, double *rVec, double *denVec, double R, double M, double d, double dr, int FHSFlag) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *phi2Vec, *phi3Vec, *phiv2Vec;
	double *n2Func, *n3Func, *v2Func;

	nGrids = (int)(R/dr)+1;
	n2Vec = (double*)malloc(nGrids*sizeof(double));
	n3Vec = (double*)malloc(nGrids*sizeof(double));
	v2Vec = (double*)malloc(nGrids*sizeof(double));
	phi2Vec = (double*)malloc(nGrids*sizeof(double));
	phi3Vec = (double*)malloc(nGrids*sizeof(double));
	phiv2Vec = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		n2Vec[i] = 0;
		n3Vec[i] = 0;
		v2Vec[i] = 0;
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}

	//for (i=(int)(d*0.5/dr)+1; i<nGrids-(int)(d*0.5/dr)-1; i++) 
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		rpMin = r-d*0.5; 
		rpMax = r+d*0.5;

		if (rpMin < 0) {
			rpMin = 0;
		}
		if (rpMax > R) {
			rpMax = R; 
		}

		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n2Func[j] = den;
			n3Func[j] = den*((0.25*(d*d))-((rp-r)*(rp-r)));
			v2Func[j] = den*(rp-r);
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] += n2*(PI*d);
		n3Vec[i] += n3*(PI);
		v2Vec[i] += v2*(2*PI);
	}

	// phi2Vec time
	for (i=0; i<nGrids; i++) {
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}

	//for (i=(int)(0.5/dr); i<nGrids-(int)(0.5/dr)-1; i++) {
	for (i=0; i<nGrids; i++) {
		r = rVec[i];

		rpMin = r-d*0.5; 
		rpMax = r+d*0.5;

		if (rpMin < 0) {
			rpMin = 0;
		}
		if (rpMax > R) {
			rpMax = R; 
		}


		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));

		for (j=0; j<N; j++) {
			rp = rVec[rpMinIndex+j];
			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];
			n2Func[j] = phip_chain_n2(n2,n3,v2,M,d,FHSFlag);
			n3Func[j] = ((d*d*0.25)-((rp-r)*(rp-r)))*phip_chain_n3(n2,n3,v2,M,d,FHSFlag);
			v2Func[j] = (rp-r)*phip_chain_v2(n2,n3,v2,M,d,FHSFlag);
		}
		phi2 = trapz(n2Func,rpMin,rpMax,N);
		phi3 = trapz(n3Func,rpMin,rpMax,N);
		phiv2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);
			
		phi2Vec[i] += phi2*(PI*d);
		phi3Vec[i] += phi3*(PI);
		phiv2Vec[i] += phiv2*(2*PI);
	}

	//for (i=(int)(0.5/dr); i<nGrids-(int)(0.5/dr)-1; i++) {
	for (i=0; i<nGrids; i++) {
		c1ChainVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(phi2Vec);
	free(phi3Vec);
	free(phiv2Vec);
}

/***************************************/

double FChainFlat(double *rVec, double *denVec, double R, double M, double d, double dr, double rc, int FHSFlag) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *n2Func, *n3Func, *v2Func;
	double *FFunc;
	double F = 0;
	double n0, zeta, X, D, ghs;

	nGrids = (int)(R/dr)+1;
	n2Vec = (double*)malloc(nGrids*sizeof(double));
	n3Vec = (double*)malloc(nGrids*sizeof(double));
	v2Vec = (double*)malloc(nGrids*sizeof(double));
	FFunc = (double*)malloc(nGrids*sizeof(double));

	for (i=0; i<nGrids; i++) {
		n2Vec[i] = 0;
		n3Vec[i] = 0;
		v2Vec[i] = 0;
		FFunc[i] = 0;
	}


	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		rpMin = r-d*0.5; 
		rpMax = r+d*0.5;
		if (rpMin < 0) {
			rpMin = 0;
				}
		if (rpMax > R) {
			rpMax = R; }

		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n2Func[j] = den;
			n3Func[j] = den*((0.25*(d*d))-((rp-r)*(rp-r)));
			v2Func[j] = den*(rp-r);
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] += n2*(PI*d);
		n3Vec[i] += n3*(PI);
		v2Vec[i] += v2*(2*PI);

	}

	for (i=0;i<(int)(rc*200);i++) {
		n2Vec[i] = n2Vec[(int)(rc*200)];
		n3Vec[i] = n3Vec[(int)(rc*200)];
		v2Vec[i] = v2Vec[(int)(rc*200)];
	}
	for (i=nGrids-(int)(rc*200);i<nGrids;i++) {
		n2Vec[i] = n2Vec[nGrids-(int)(rc*200)-1];
		n3Vec[i] = n3Vec[nGrids-(int)(rc*200)-1];
		v2Vec[i] = v2Vec[nGrids-(int)(rc*200)-1];
	}



	// FFunc time
	for (i=0; i<nGrids; i++) {
		n2 = n2Vec[i];
		n3 = n3Vec[i];
		v2 = v2Vec[i];

		n0 = n2/PI;
		zeta = 1 - (v2*v2/n2/n2);

		ghs = get_ghs(n2, n3, v2, d);

		FFunc[i] = ((1-M)/M)*n0*zeta*(log(ghs));

	}


	for (i=0;i<(int)(rc*200);i++) {
		FFunc[i] = FFunc[(int)(rc*200)];
	}
	for (i=nGrids-(int)(rc*200);i<nGrids;i++) {
		FFunc[i] = FFunc[nGrids-(int)(rc*200)-1];
	}


	F = trapz(FFunc, 0, R, nGrids);

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(FFunc);

	return F;
}
/*******************************************************/


