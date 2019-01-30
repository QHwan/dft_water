#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "constant.h"
#include "math2.h"


double fphi2(double n2,double n3,double nv2,double d, int FHSFlag) {
	double fphi2 = 0.0;
	//printf("%.3f \t %.3f \t %.3f \n", n2, n3, nv2);
	if (n3==0) {
		fphi2 = 0.;
	}
	else {
		if (FHSFlag == 1) {
			fphi2 = (-log(1-n3)/(PI*d*d)) + n2/(PI*d*(1-n3)) + ((n2*n2-nv2*nv2)/(12*PI*(pow(n3,3.))))*((n3*log(1-n3))+(n3*n3/((1-n3)*(1-n3))));
		}
		else if (FHSFlag == 2) {
			fphi2 = (-log(1-n3)/(PI*d*d)) + n2/(PI*d*(1-n3)) + (1/8/PI)*(1/(1-n3)/(1-n3))*(n2*n2-nv2*nv2);
		}
	}
	return fphi2;
}

double fphi3(double n2,double n3,double nv2,double d, int FHSFlag) {
	double fphi3 = 0.0;
	if (n3==0) {
		fphi3 = 0.;
	}
	else {
		if (FHSFlag == 1) {
			fphi3 = (n2)/((1-n3)*(PI*d*d)) + (n2*n2-nv2*nv2)/((2*PI*d)*((1-n3)*(1-n3))) + ((-2*log(1-n3)/(pow(n3,3.)))+((-2+5*n3-n3*n3)/(n3*n3*(pow((1-n3),3.)))))*((pow(n2,3.)-3*n2*nv2*nv2)/(36*PI));
		}
		else if (FHSFlag == 2) {
			fphi3 = n2/(1-n3*PI*d*d) + (n2*n2-nv2*nv2)/((2*PI*d)*((1-n3)*(1-n3))) + (1/(12*PI))*(1/(pow(1-n3,3)))*(pow(n2,3)-3*n2*nv2*nv2);
		}
	}
	return fphi3;
}

double fphiv2(double n2,double n3,double nv2,double d, int FHSFlag) {
	double fphiv2 = 0.0;
	if (n3==0) {
		fphiv2 = 0.;
	}
	else {
		if (FHSFlag == 1) {
			fphiv2 = (-nv2)/((PI*d)*(1-n3)) - (n2*nv2)/(6*PI*(pow(n3,3.)))*(n3*log(1-n3)+((n3*n3)/((1-n3)*(1-n3))));
		}
		else if (FHSFlag == 2) {
			fphiv2 = (-nv2)/((PI*d)*(1-n3)) - (1/(4*PI))*(1/((1-n3)*(1-n3)))*(nv2*nv2);
		}
	}
	return fphiv2;
}

/*********************************************/

void c1HSFlat(double *c1HSVec, double *rVec, double *denVec, double R, double d, double dr, int FHSFlag) {
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

	for (i=round(d*0.5/dr); i<nGrids-round(d*0.5/dr); i++) {
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
			n2Func[j] = fphi2(n2,n3,v2,d,FHSFlag);
			n3Func[j] = ((d*d*0.25)-((rp-r)*(rp-r)))*fphi3(n2,n3,v2,d,FHSFlag);
			v2Func[j] = (rp-r)*fphiv2(n2,n3,v2,d,FHSFlag);
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

	for (i=0; i<nGrids; i++) {
		c1HSVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}
/*
	for (i=0; i<nGrids; i++) {
		printf("%f\n",c1HSVec[i]);
	}
	exit(1);
	*/
	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(phi2Vec);
	free(phi3Vec);
	free(phiv2Vec);
}

/***************************************/

void c1HSSpherical(double *c1HSVec, double *rVec, double *denVec, double R, double d, double dr, int FHSFlag) {
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


	for (i=1; i<nGrids; i++) {
		r = rVec[i];
		if (r >= d*0.5) {
			rpMin = r-d*0.5; 
		} else {
			rpMin = d*0.5-r; 
		}
		rpMax = r+d*0.5;
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
			n3Func[j] = den*rp*((0.25*(d*d))-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*(r*r-rp*rp+(d*d*0.25));
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] += n2*(PI*d/r);
		n3Vec[i] += n3*(PI/r);
		v2Vec[i] += v2*(PI/(r*r));

	}

	// n3Vec calculate alone
	/*
	for (i=0; i<round(d*0.5/dr); i++) {
		r = rVec[i];
		rpMin = 0.;
		rpMax = d*0.5 - r;
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n3Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n3Func[j] = (rp*rp)*den;
		}
		n3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		n3Vec[i] += n3*4*PI;
	}
	*/

	n2Vec[0] = n2Vec[1];
	v2Vec[0] = v2Vec[1];



	// phi2Vec time
	for (i=0; i<nGrids; i++) {
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}

//	for (i=1; i<nGrids-round(d*0.5/dr); i++) {
	for (i=1; i<nGrids; i++) {
		r = rVec[i];

		if (r >= d*0.5) {
			rpMin = (r-d*0.5); 
		} else {
			rpMin = (d*0.5-r); 
		}
		rpMax = r+d*0.5;
		if (rpMax > R) {
			rpMax = R; }
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
			n2Func[j] = rp*fphi2(n2,n3,v2,d,FHSFlag);
			n3Func[j] = rp*((d*d*0.25)-((r-rp)*(r-rp)))*fphi3(n2,n3,v2,d,FHSFlag);
			v2Func[j] = (r*r-rp*rp+((d*d)*0.25))*fphiv2(n2,n3,v2,d,FHSFlag);
		}
		phi2 = trapz(n2Func,rpMin,rpMax,N);
		phi3 = trapz(n3Func,rpMin,rpMax,N);
		phiv2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);
			
		phi2Vec[i] += phi2*(PI*d/r);
		phi3Vec[i] += phi3*(PI/r);
		phiv2Vec[i] += phiv2*(PI/(r));
	}

	// phi3Vec calculate alone
	/*
	for (i=0; i<round(d*0.5/dr); i++) {
		r = rVec[i];
		
		rpMin = 0;
		rpMax = d*0.5 - r;
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n3Func = (double*)malloc(N*sizeof(double));

		for (j=0; j<N; j++) {
			rp = rVec[rpMinIndex+j];
			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];

			n3Func[j] = (rp*rp)*fphi3(n2,n3,v2,d,FHSFlag);
		}
		phi3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		phi3Vec[i] += phi3*4*PI;
	}
	*/
	phi2Vec[0] = phi2Vec[1];
	phiv2Vec[0] = phiv2Vec[1];

/*
	FILE *f;
		f = fopen("ph.xvg", "w");
		if(f==NULL) { printf("File error!"); }
		else {
			for (i=0; i<nGrids; i++) {
				fprintf(f, "%d \t %lf\t%lf\t%lf \n", i, phi2Vec[i], phi3Vec[i], phiv2Vec[i]);
			}
		}
		fclose(f);exit(1);
		*/

	for (i=0; i<nGrids; i++) {
		c1HSVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}
/*
	for (i=0;i<nGrids;i++) {
		printf("%f\n",c1HSVec[i]);
	}
	exit(1);
	*/

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(phi2Vec);
	free(phi3Vec);
	free(phiv2Vec);

}

/***************************************/

void c1HSSphericalSolute(double *c1HSVec, double *rVec, double *denVec, double R, double RSolute, double d, double dr, int FHSFlag) {
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

	for (i=round(RSolute/dr); i<nGrids; i++) {
		r = rVec[i];
		rpMin = r-d*0.5; 
		if (rpMin < RSolute) {
			rpMin = RSolute;
		}
		rpMax = r+d*0.5;
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
			n2Func[j] = den*rp;
			n3Func[j] = den*rp*((0.25*(d*d))-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*(r*r-rp*rp+(d*d*0.25));
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] += n2*(PI*d/r);
		n3Vec[i] += n3*(PI/r);
		v2Vec[i] += v2*(PI/(r*r));

	}

	// phi2Vec time
	for (i=0; i<nGrids; i++) {
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}

	for (i=round(RSolute/dr)+round(d*0.5/dr); i<nGrids-round(d*0.5/dr); i++) {
		r = rVec[i];

		rpMin = r-d*0.5; 
		if (rpMin < RSolute) {
			rpMin = RSolute;
		}
		rpMax = r+d*0.5;
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
			n2Func[j] = rp*fphi2(n2,n3,v2,d,FHSFlag);
			n3Func[j] = rp*((d*d*0.25)-((r-rp)*(r-rp)))*fphi3(n2,n3,v2,d,FHSFlag);
			v2Func[j] = (r*r-rp*rp+((d*d)*0.25))*fphiv2(n2,n3,v2,d,FHSFlag);
		}
		phi2 = trapz(n2Func,rpMin,rpMax,N);
		phi3 = trapz(n3Func,rpMin,rpMax,N);
		phiv2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);
			
		phi2Vec[i] += phi2*(PI*d/r);
		phi3Vec[i] += phi3*(PI/r);
		phiv2Vec[i] += phiv2*(PI/(r));
	}
	/*
	FILE *f;
		f = fopen("ph.xvg", "w");
		if(f==NULL) { printf("File error!"); }
		else {
			for (i=0; i<nGrids; i++) {
				fprintf(f, "%d \t %lf\t%lf\t%lf \n", i, phi2Vec[i], phi3Vec[i], phiv2Vec[i]);
			}
		}
		fclose(f);exit(1);
		*/

	for (i=0; i<nGrids; i++) {
		c1HSVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}
/*
	for (i=0;i<nGrids;i++) {
		printf("%f\n",c1HSVec[i]);
	}
	exit(1);
	*/

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(phi2Vec);
	free(phi3Vec);
	free(phiv2Vec);

}

/***************************************/

double FHSSpherical(double *rVec, double *denVec, double R, double d, double dr, int FHSFlag) {
	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *n2Func, *n3Func, *v2Func;
	double *FFunc;
	double F = 0;

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


	for (i=1; i<nGrids; i++) {
		r = rVec[i];
		if (r >= d*0.5) {
			rpMin = r-d*0.5; 
		} else {
			rpMin = d*0.5-r; 
		}
		rpMax = r+d*0.5;
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
			n3Func[j] = den*rp*((0.25*(d*d))-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*(r*r-rp*rp+(d*d*0.25));
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free(n2Func);
		free(n3Func);
		free(v2Func);

		n2Vec[i] += n2*(PI*d/r);
		n3Vec[i] += n3*(PI/r);
		v2Vec[i] += v2*(PI/(r*r));

	}

	// n3Vec calculate alone
	for (i=0; i<round(d*0.5/dr); i++) {
		r = rVec[i];
		rpMin = 0.;
		rpMax = d*0.5 - r;
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n3Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n3Func[j] = (rp*rp)*den;
		}
		n3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		n3Vec[i] += n3*4*PI;
	}
	n2Vec[0] = n2Vec[1];
	v2Vec[0] = v2Vec[1];


	// FFunc time
	for (i=0; i<nGrids-round(d*0.5/dr); i++) {
		n2 = n2Vec[i];
		n3 = n3Vec[i];
		v2 = v2Vec[i];
		r = rVec[i];

		if (n3 >= 1E-8) {
			FFunc[i] = 4*PI*r*r*(-(n2/(PI*d*d)*log(1-n3))+((1/(2*PI*d))*((n2*n2-v2*v2)/(1-n3)))+(((n2*n2*n2-(3*n2*v2*v2))/(36*PI*n3*n3*n3))*((n3*log(1-n3))+(n3*n3/((1-n3)*(1-n3))))));
		} else {
			FFunc[i] = 0;
		}

	}

	F = trapz(FFunc, 0, R, nGrids);

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(FFunc);

	return F;
}


