#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "constant.h"
#include "math2.h"
#include "memory.h"

//#define CUBE(x) ((x)*(x)*(x))

double fphi2(double n2,double n3,double nv2) {
	double fphi2 = 0.0;

	if (n3==0) {
		fphi2 = 0.;
	}
	else {
//		fphi2 = (-log(1-n3)/(PI)) + n2/(PI*(1-n3)) + ((n2*n2-nv2*nv2)/(12*PI*(n3*n3*n3)))*((n3*log(1-n3))+(n3*n3/((1-n3)*(1-n3))));
		fphi2 = -(1/PI)*log(1-n3) + (1/PI)*n2/(1-n3) + (n2*n2-nv2*nv2)/(12*PI*n3*n3*n3)*(n3*log(1-n3) + n3*n3/(1-n3)/(1-n3));
	}
	return fphi2;
}

double fphi3(double n2,double n3,double nv2) {
	double fphi3 = 0.0;
	
	if (n3==0) {
		fphi3 = 0.; 
	} else {
//		fphi3 = (n2)/((1-n3)*(PI)) + (n2*n2-nv2*nv2)*0.5/((PI)*((1-n3)*(1-n3))) + ((-2*log(1-n3)/(n3*n3*n3))+((-2+5*n3-n3*n3)/(n3*n3*(((1-n3)*(1-n3)*(1-n3))))))*((n2*n2*n2-3*n2*nv2*nv2)/(36*PI));
		fphi3 = (1./PI)*n2/(1.-n3) + (1/2./PI)*(n2*n2-nv2*nv2)/(1-n3)/(1-n3) + (n2*n2*n2-3*n2*nv2*nv2)/36./PI*(-2.*log(1.-n3)/n3/n3/n3 + (-n3*n3+5.*n3-2.)/n3/n3/(1-n3)/(1-n3)/(1-n3));

	} 
	return fphi3;
}

double fphiv2(double n2,double n3,double nv2) {
	double fphiv2 = 0.0;
	if (n3==0) {
		fphiv2 = 0.;
	}
	else {
		//fphiv2 = (-nv2)/((PI)*(1-n3)) - (n2*nv2)/(6*PI*(n3*n3*n3))*(n3*log(1-n3)+((n3*n3)/((1-n3)*(1-n3))));
		fphiv2 = -(1./PI)*nv2/(1-n3) - n2*nv2/6./PI/n3/n3/n3*(n3*log(1.-n3) + n3*n3/(1-n3)/(1-n3));
	}
	return fphiv2;
}

/*********************************************/

void c1HSFlat(double *c1HSVec, double *rVec, double *denVec, double R, double rc) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *phi2Vec, *phi3Vec, *phiv2Vec;
	double *n2Func, *n3Func, *v2Func;
	double dr;

	dr = rVec[1]-rVec[0];

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

	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		rpMin = r-0.5; 
		rpMax = r+0.5;

		if (rpMin < 0) {
			rpMin = 0;
		}
		if (rpMax > R) {
			rpMax = R; 
		}

		rpMinIndex = (round)(rpMin/dr);
		rpMaxIndex = (round)(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n2Func[j] = den;
			n3Func[j] = den*((0.25)-((rp-r)*(rp-r)));
			v2Func[j] = den*(rp-r);
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] = n2*(PI);
		n3Vec[i] = n3*(PI);
		v2Vec[i] = v2*(2*PI);
	}

	// phi2Vec time
	for (i=0; i<nGrids; i++) {
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}

	for (i=(int)(0.5/dr); i<nGrids-(int)(0.5/dr)-1; i++) {
	//for (i=0; i<nGrids; i++) {
		r = rVec[i];

		rpMin = r-0.5; 
		rpMax = r+0.5;

		if (rpMin < 0) {
			rpMin = 0;
		}
		if (rpMax > R) {
			rpMax = R; 
		}


		rpMinIndex = (round)(rpMin/dr);
		rpMaxIndex = (round)(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));

		for (j=0; j<N; j++) {
			rp = rVec[rpMinIndex+j];
			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];
			n2Func[j] = fphi2(n2,n3,v2);
			n3Func[j] = ((0.25)-((rp-r)*(rp-r)))*fphi3(n2,n3,v2);
			v2Func[j] = (rp-r)*fphiv2(n2,n3,v2);
		}
		phi2 = trapz(n2Func,rpMin,rpMax,N);
		phi3 = trapz(n3Func,rpMin,rpMax,N);
		phiv2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);
			
		phi2Vec[i] += phi2*(PI);
		phi3Vec[i] += phi3*(PI);
		phiv2Vec[i] += phiv2*(2*PI);
	}

	for (i=((int)round(0.5/dr)); i<nGrids-(int)round(0.5/dr)-1; i++) {
	//for (i=0; i<nGrids; i++) {
		c1HSVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(phi2Vec);
	free(phi3Vec);
	free(phiv2Vec);
}

/***************************************/

void c1HSSpherical(double *c1HSVec, double *rVec, double *denVec, double R, double rc) {
	int i, j, range;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *phi2Vec, *phi3Vec, *phiv2Vec;
	double *n2Func, *n3Func, *v2Func;
	double dr;


	dr = rVec[1]-rVec[0];

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


	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		if (r >= 0.5) {
			rpMin = r-0.5; 
		} else {
			rpMin = 0.5-r; 
		}
		rpMax = r+0.5;
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

			n2Func[j] = den*rp;
			n3Func[j] = den*rp*((0.25)-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*((0.25)-((r-rp)*(r-rp)));
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] += n2*(PI/r);
		n3Vec[i] += n3*(PI/r);
		v2Vec[i] += v2*(PI/(r*r));

	}


	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		if (r >= 0.5) {
			rpMin = r-0.5; 
		} else {
			rpMin = 0.5-r; 
		}
		rpMax = r+0.5;
		if (rpMax > R) {
			rpMax = R; }

		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		v2Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			v2Func[j] = den*rp*(r-rp);
		}

		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(v2Func);

		v2Vec[i] += v2*(2*PI/(r));

	}

	// n3Vec calculate alone
	for (i=0; i<round(0.5/dr); i++) {
		r = rVec[i];
		rpMin = 0.;
		rpMax = 0.5 - r;
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n3Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n3Func[j] = 4*PI*(rp*rp)*den;
		}
		n3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		n3Vec[i] += n3;
	}

	n2Vec[0] = n2Vec[1];
	n3Vec[0] = n3Vec[1];
	v2Vec[0] = v2Vec[1];



	// phi2Vec time
	for (i=0; i<nGrids; i++) {
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}

	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		rpMin = r-0.5;
		if (rpMin < 0) {
			rpMin = 0.5-r; }
		rpMax = r+0.5;
		if (rpMax > R) {
			rpMax = R; }
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));

		for (j=0; j<N; j++) {

			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];

			n2Func[j] = rp*fphi2(n2,n3,v2);
			n3Func[j] = rp*((0.25)-((r-rp)*(r-rp)))*fphi3(n2,n3,v2);
			v2Func[j] = rp*(0.25-(r-rp)*(r-rp))*fphiv2(n2,n3,v2);
		}
		phi2 = trapz(n2Func,rpMin,rpMax,N);
		phi3 = trapz(n3Func,rpMin,rpMax,N);
		phiv2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);
			
		phi2Vec[i] += phi2*(PI/r);
		phi3Vec[i] += phi3*(PI/r);
		phiv2Vec[i] += phiv2*(PI/(r*r));
	}
	

	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		rpMin = r-0.5;
		if (rpMin < 0) {
			rpMin = 0.5-r; }
		rpMax = r+0.5;
		if (rpMax > R) {
			rpMax = R; }
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		v2Func = (double*)malloc(N*sizeof(double));

		for (j=0; j<N; j++) {

			v2 = v2Vec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];

			v2Func[j] = rp*(r-rp)*fphiv2(n2,n3,v2);
		}
		phiv2 = trapz(v2Func,rpMin,rpMax,N);

		free(v2Func);
			
		phiv2Vec[i] += phiv2*(2*PI/(r));
	}

	// phi3Vec calculate alone
	for (i=0; i<round(0.5/dr); i++) {
		r = rVec[i];
		
		rpMin = 0;
		rpMax = 0.5 - r;
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n3Func = (double*)malloc(N*sizeof(double));

		for (j=0; j<N; j++) {
			rp = rVec[rpMinIndex+j];
			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];

			n3Func[j] = 4.*PI*(rp*rp)*fphi3(n2,n3,v2);
		}
		phi3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		phi3Vec[i] += phi3;
	}
	
	phi2Vec[0] = phi2Vec[1];
	phi3Vec[0] = phi3Vec[1];
	phiv2Vec[0] = phiv2Vec[1];


	for (i=0; i<nGrids; i++) {
		c1HSVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}


	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(phi2Vec);
	free(phi3Vec);
	free(phiv2Vec);

}

/***************************************/

void c1HSSphericalSolute(double *c1HSVec, double *rVec, double *denVec, double R, double RSolute, double d,double rc, int FHSFlag) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *phi2Vec, *phi3Vec, *phiv2Vec;
	double *n2Func, *n3Func, *v2Func;
	double dr;

	dr = rVec[1]-rVec[0];

	nGrids = (int)((R-RSolute)/dr)+1;
	n2Vec = dVector(nGrids);
	n3Vec = dVector(nGrids);
	v2Vec = dVector(nGrids);
	phi2Vec = dVector(nGrids);
	phi3Vec = dVector(nGrids);
	phiv2Vec = dVector(nGrids);
	

	for (i=0; i<nGrids; i++) {
	//for (i=1; i<nGrids; i++) {
		r = rVec[i];
		rpMin = r - 0.5; 
		rpMax = r + 0.5;
		if (rpMin < RSolute) {
			rpMin = RSolute;
		}
		if (rpMax > R) {
			rpMax = R; 
		}

		rpMinIndex = i-(int)(0.5/dr);
		rpMaxIndex = i+(int)(0.5/dr);
		if (rpMinIndex < 0) {
			rpMinIndex = 0;
		}
		if (rpMaxIndex > nGrids) {
			rpMaxIndex = nGrids;
		}
		/*
		rpMinIndex = (int)((rpMin-RSolute)/dr);
		rpMaxIndex = (int)((rpMax-RSolute)/dr);
		*/
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = dVector(N);
		n3Func = dVector(N);
		v2Func = dVector(N);
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n2Func[j] = den*rp;
			n3Func[j] = den*rp*(0.25 - ((r-rp)*(r-rp)));
			v2Func[j] = den*rp*(r*r - rp*rp + 0.25);
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free_dVector(n2Func);
		free_dVector(n3Func);
		free_dVector(v2Func);

		n2Vec[i] += n2*(PI/r);
		n3Vec[i] += n3*(PI/r);
		v2Vec[i] += v2*(PI/(r*r));
	}



	// phi2Vec time
	//
	for (i=0; i<nGrids; i++) {
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}


//	for (i=round(RSolute/dr); i<nGrids; i++) {
	for (i=0; i<nGrids; i++) {
		r = rVec[i];

		rpMin = r - 0.5; 
		rpMax = r + 0.5;
		if (rpMin < RSolute) {
			rpMin = RSolute;
		}
		if (rpMax > R) {
			rpMax = R; 
		}

		rpMinIndex = i-(int)(0.5/dr);
		rpMaxIndex = i+(int)(0.5/dr);
		if (rpMinIndex < 0) {
			rpMinIndex = 0;
		}
		if (rpMaxIndex > nGrids) {
			rpMaxIndex = nGrids;
		}
		/*
		rpMinIndex = (int)((rpMin-RSolute)/dr);
		rpMaxIndex = (int)((rpMax-RSolute)/dr);
		*/
		N = (rpMaxIndex-rpMinIndex)+1;


		n2Func = dVector(N);
		n3Func = dVector(N);
		v2Func = dVector(N);

		for (j=0; j<N; j++) {
			rp = rVec[rpMinIndex+j];
			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];
			n2Func[j] = rp*fphi2(n2,n3,v2);
			n3Func[j] = rp*(0.25 - ((r-rp)*(r-rp)))*fphi3(n2,n3,v2);
			v2Func[j] = rp*(r*r - rp*rp + 0.25)*fphiv2(n2,n3,v2);
		}
		phi2 = trapz(n2Func,rpMin,rpMax,N);
		phi3 = trapz(n3Func,rpMin,rpMax,N);
		phiv2 = trapz(v2Func,rpMin,rpMax,N);

		free_dVector(n2Func);
		free_dVector(n3Func);
		free_dVector(v2Func);
			
		phi2Vec[i] += phi2*(PI/r);
		phi3Vec[i] += phi3*(PI/r);
		phiv2Vec[i] += phiv2*(PI/(r*r));
	}

	for (i=nGrids-(int)(rc*200);i<nGrids;i++) {
		phi2Vec[i] = phi2Vec[nGrids-(int)(rc*200+1)];
		phi3Vec[i] = phi3Vec[nGrids-(int)(rc*200+1)];
		phiv2Vec[i] = phiv2Vec[nGrids-(int)(rc*200+1)];
	}


	for (i=0; i<nGrids; i++) {
		c1HSVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}



	free_dVector(n2Vec);
	free_dVector(n3Vec);
	free_dVector(v2Vec);
	free_dVector(phi2Vec);
	free_dVector(phi3Vec);
	free_dVector(phiv2Vec);

}

/***************************************/
double FHSFlat(double *rVec, double *denVec, double R, double d, double rc, int FHSFlag) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *n2Func, *n3Func, *v2Func;
	double *FFunc;
	double F = 0;
	double dr;

	dr = rVec[1]-rVec[0];

	nGrids = (int)(R/dr)+1;

	n2Vec = dVector(nGrids);
	n3Vec = dVector(nGrids);
	v2Vec = dVector(nGrids);
	FFunc = dVector(nGrids);


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

		n2Func = dVector(N);
		n3Func = dVector(N);
		v2Func = dVector(N);
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

		free_dVector(n2Func);
		free_dVector(n3Func);
		free_dVector(v2Func);

		n2Vec[i] += n2*(PI*d);
		n3Vec[i] += n3*(PI);
		v2Vec[i] += v2*(2*PI);
	}
	

	// FFunc time
	for (i=((int)round(0.5/dr)); i<nGrids-(int)round(0.5/dr)-1; i++) {
	//for (i=0; i<nGrids; i++) {
		n2 = n2Vec[i];
		n3 = n3Vec[i];
		v2 = v2Vec[i];
		r = rVec[i];

		if (n3 >= 1E-8) {
			FFunc[i] = (-(n2/(PI*d*d)*log(1-n3))+((1/(2*PI*d))*((n2*n2-v2*v2)/(1-n3)))+(((n2*n2*n2-(3*n2*v2*v2))/(36*PI*n3*n3*n3))*((n3*log(1-n3))+(n3*n3/((1-n3)*(1-n3))))));
		} else {
			FFunc[i] = 0;
		}

	}


	/*
	for (i=0;i<(int)(rc*200);i++) {
		FFunc[i] = FFunc[(int)(rc*200)];
	}
	for (i=nGrids-(int)(rc*200);i<nGrids;i++) {
		FFunc[i] = FFunc[nGrids-(int)(rc*200)-1];
	}
	*/


	F = trapz(FFunc, 0, R, nGrids);

	free_dVector(n2Vec);
	free_dVector(n3Vec);
	free_dVector(v2Vec);
	free_dVector(FFunc);

	return F;
}
/*******************************************************/


void FHSSpherical(double *FHSVec, double *rVec, double *denVec, double R, double d, double rc, int FHSFlag) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *n2Func, *n3Func, *v2Func;
	double *FFunc;
	double F = 0;
	double dr;

	dr = rVec[1]-rVec[0];

	nGrids = (int)(R/dr)+1;

	n2Vec = dVector(nGrids);
	n3Vec = dVector(nGrids);
	v2Vec = dVector(nGrids);
	FFunc = dVector(nGrids);

	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		if (r >= 0.5) {
			rpMin = r-0.5; 
		} else {
			rpMin = 0.5-r; 
		}
		rpMax = r+0.5;
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

			n2Func[j] = den*rp;
			n3Func[j] = den*rp*((0.25)-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*((0.25)-((r-rp)*(r-rp)));
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] += n2*(PI/r);
		n3Vec[i] += n3*(PI/r);
		v2Vec[i] += v2*(PI/(r*r));

	}


	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		if (r >= 0.5) {
			rpMin = r-0.5; 
		} else {
			rpMin = 0.5-r; 
		}
		rpMax = r+0.5;
		if (rpMax > R) {
			rpMax = R; }

		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		v2Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			v2Func[j] = den*rp*(r-rp);
		}

		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(v2Func);

		v2Vec[i] += v2*(2*PI/(r));

	}

	// n3Vec calculate alone
	for (i=0; i<round(0.5/dr); i++) {
		r = rVec[i];
		rpMin = 0.;
		rpMax = 0.5 - r;
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n3Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n3Func[j] = 4*PI*(rp*rp)*den;
		}
		n3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		n3Vec[i] += n3;
	}

	n2Vec[0] = n2Vec[1];
	n3Vec[0] = n3Vec[1];
	v2Vec[0] = v2Vec[1];



	// FFunc time
	for (i=0; i<nGrids; i++) {
		n2 = n2Vec[i];
		n3 = n3Vec[i];
		v2 = v2Vec[i];
		r = rVec[i];

		if (n3 >= 1E-8) {
			FHSVec[i] = (-(n2/(PI*d*d)*log(1-n3))+((1/(2*PI*d))*((n2*n2-v2*v2)/(1-n3)))+(((n2*n2*n2-(3*n2*v2*v2))/(36*PI*n3*n3*n3))*((n3*log(1-n3))+(n3*n3/((1-n3)*(1-n3))))));
		} else {
			FHSVec[i] = 0;
		}


	}





	free_dVector(n2Vec);
	free_dVector(n3Vec);
	free_dVector(v2Vec);
	free_dVector(FFunc);

}


/**********************************************************/
double FHSSphericalSolute(double *rVec, double *denVec, double R, double RSolute, double d, double rc,int FHSFlag) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *n2Func, *n3Func, *v2Func;
	double *FFunc;
	double F = 0;
	double dr;

	dr = rVec[1]-rVec[0];

	nGrids = (int)((R-RSolute)/dr)+1;
	n2Vec = dVector(nGrids);
	n3Vec = dVector(nGrids);
	v2Vec = dVector(nGrids);
	FFunc = dVector(nGrids);

	for (i=0; i<nGrids; i++) {
	//for (i=1; i<nGrids; i++) {
		r = rVec[i];
		rpMin = r - 0.5; 
		rpMax = r + 0.5;
		if (rpMin < RSolute) {
			rpMin = RSolute;
		}
		if (rpMax > R) {
			rpMax = R; }

		rpMinIndex = i-(int)(0.5/dr);
		rpMaxIndex = i+(int)(0.5/dr);
		if (rpMinIndex < 0) {
			rpMinIndex = 0;
		}
		if (rpMaxIndex > nGrids) {
			rpMaxIndex = nGrids;
		}
		/*
		rpMinIndex = (floor)((rpMin-RSolute)/dr);

		rpMaxIndex = (floor)((rpMax-RSolute)/dr);
		*/
		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = dVector(N);
		n3Func = dVector(N);
		v2Func = dVector(N);
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n2Func[j] = den*rp;
			n3Func[j] = den*rp*(0.25-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*(r*r-rp*rp+0.25);
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free_dVector(n2Func);
		free_dVector(n3Func);
		free_dVector(v2Func);

		n2Vec[i] = n2*(PI/r);
		n3Vec[i] = n3*(PI/r);
		v2Vec[i] = v2*(PI/(r*r));
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
		r = rVec[i];

		if (n3 >= 1E-8) {
			FFunc[i] = 4*PI*r*r*(-(n2/(PI*d*d)*log(1-n3))+((1/(2*PI*d))*((n2*n2-v2*v2)/(1-n3)))+(((n2*n2*n2-(3*n2*v2*v2))/(36*PI*n3*n3*n3))*((n3*log(1-n3))+(n3*n3/((1-n3)*(1-n3))))));

		//	FFunc[i] = 4*PI*r*r*(-(n2/(PI)*log(1-n3))+((1/(2*PI))*((n2*n2-v2*v2)/(1-n3)))+(((n2*n2*n2-(3*n2*v2*v2))/(36*PI*n3*n3*n3))*((n3*log(1-n3))+(n3*n3/((1-n3)*(1-n3))))));
//			FFunc[i] = 4*PI*r*r*(-n2/PI*log(1-n3));
//			FFunc[i] += 4*PI*r*r*(1./2./PI*n2*n2/(1-n3));
//			FFunc[i] += 4*PI*r*r*n2*n2*n2*(log(1-n3)/36./PI/n3/n3 + 1/36./PI/n3/(1-n3)/(1-n3));
//			FFunc[i] += 4*PI*r*r*(-1./2./PI*v2*v2/(1-n3));
//			FFunc[i] += 4*PI*r*r*(-1)*(n2*v2*v2*(log(1-n3)/12./PI/n3/n3 + 1/12./PI/n3/(1-n3)/(1-n3)));
		} else {
			FFunc[i] = 0;
		}

	}

	/*
	for (i=0; i<nGrids; i++) {
		printf("%f\n",FFunc[i]);
	}
	exit(1);
	*/



	F = trapz(FFunc, RSolute, R, nGrids);

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(FFunc);

	return F;
}


