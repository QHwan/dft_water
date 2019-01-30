#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "constant.h"
#include "math2.h"
#include "memory.h"

double get_D(double n2,double n3,double v2,double K,double T_assoc,double d) {
	double ghs;
	ghs = get_ghs(n2,n3,v2,d);
	return 4*PI*K*ghs*(exp(1./T_assoc)-1);
}

double get_Dp2(double n2,double n3,double v2,double K,double T_assoc,double d) {
	double ghsp2;
	ghsp2 = get_ghsp2(n2,n3,v2,d);
	return 4*PI*K*ghsp2*(exp(1./T_assoc)-1);
}

double get_Dp3(double n2,double n3,double v2,double K,double T_assoc, double d) {
	double ghsp3;
	ghsp3 = get_ghsp3(n2,n3,v2,d);
	return 4*PI*K*ghsp3*(exp(1./T_assoc)-1);
}

double get_Dpv2(double n2,double n3,double v2,double K,double T_assoc,double d) {
	double ghspv2;
	ghspv2 = get_ghspv2(n2,n3,v2,d);
	return 4*PI*K*ghspv2*(exp(1./T_assoc)-1);
}

double get_X(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d) {
	double n0, zeta, D;
	n0 = n2/PI;
	zeta = 1-(v2*v2/n2/n2);
	D = get_D(n2,n3,v2,K,T_assoc,d);
//	return (-1+sqrt(1+2.*Ma*zeta*n0*D))/(Ma*zeta*n0*D);

	return (-1.+sqrt(1.+4.*Ma*zeta*n0*D))/(2.*Ma*zeta*n0*D);
}

double get_Xp2(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d) {
	double n0, zeta, D, Dp2;
	n0 = n2/PI;
	zeta = 1-(v2*v2/n2/n2);
	D = get_D(n2,n3,v2,K,T_assoc,d);
	Dp2 = get_Dp2(n2,n3,v2,K,T_assoc,d);

//	return ((1/(zeta*n0*D*sqrt(1+2.*Ma*zeta*n0*D))) + ((1-sqrt(1+2.*Ma*zeta*n0*D))/(Ma*zeta*zeta*n0*n0*D*D))) * \
		(D*(1/PI/d/d)*(1+v2*v2/n2/n2) + zeta*n0*Dp2);

	return ((1./(zeta*n0*D*sqrt(1.+4.*Ma*zeta*n0*D))) + ((1.-sqrt(1.+4.*Ma*zeta*n0*D))/(2.*Ma*zeta*zeta*n0*n0*D*D))) * \
		(D*(1./PI)*(1+v2*v2/n2/n2) + zeta*n0*Dp2);
}

double get_Xp3(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d) {
	double n0, zeta, D, Dp3;
	n0 = n2/PI;
	zeta = 1-(v2*v2/n2/n2);
	D = get_D(n2,n3,v2,K,T_assoc,d);
	Dp3 = get_Dp3(n2,n3,v2,K,T_assoc,d);
//	return ((1./(zeta*n0*D*sqrt(1.+2.*Ma*zeta*n0*D))) + ((1-sqrt(1+2.*Ma*zeta*n0*D))/(Ma*zeta*zeta*n0*n0*D*D))) * \
		(D*0 + zeta*n0*Dp3);

	return ((1/(zeta*n0*D*sqrt(1+4*Ma*zeta*n0*D))) + ((1-sqrt(1+4*Ma*zeta*n0*D))/(2*Ma*zeta*zeta*n0*n0*D*D))) * \
		(D*0 + zeta*n0*Dp3);
}

double get_Xpv2(double n2, double n3, double v2, double Ma,double K,double T_assoc, double d) {
	double n0, zeta, D, Dpv2;
	n0 = n2/PI;
	zeta = 1-(v2*v2/n2/n2);
	D = get_D(n2,n3,v2,K,T_assoc,d);
	Dpv2 = get_Dpv2(n2,n3,v2,K,T_assoc,d);
//	return ((1/(zeta*n0*D*sqrt(1+2.*Ma*zeta*n0*D))) + ((1-sqrt(1+2.*Ma*zeta*n0*D))/(Ma*zeta*zeta*n0*n0*D*D))) * \
		(D*(-2/PI/d/d)*(v2/n2) + zeta*n0*Dpv2);

	return ((1/(zeta*n0*D*sqrt(1+4*Ma*zeta*n0*D))) + ((1-sqrt(1+4*Ma*zeta*n0*D))/(2*Ma*zeta*zeta*n0*n0*D*D))) * \
		(D*(-2/PI/d/d)*(v2/n2) + zeta*n0*Dpv2);
}



double phip_assoc_n2(double n2, double n3, double v2, double Ma, double K, double T_assoc, double d, int FHSFlag) {
	double fphi2 = 0.0;
	double n0, zeta, X, Xp2;

	n0 = n2/PI;
	zeta = 1-(v2*v2/n2/n2);

	/*
	if (zeta < 1e-8) {
		fphi2 = 0;
	} else {
	*/
		X = get_X(n2,n3,v2,Ma,K,T_assoc,d);
		Xp2 = get_Xp2(n2,n3,v2,Ma,K,T_assoc,d);

		fphi2 = Ma*(log(X) - 0.5*X + 0.5)*(1./PI)*(1+v2*v2/n2/n2) + \
			Ma*zeta*n0*(1./X - 0.5)*Xp2;
	//}
	return fphi2;
}

double phip_assoc_n3(double n2, double n3, double v2, double Ma, double K, double T_assoc, double d, int FHSFlag) {
	double fphi3 = 0.0;
	double n0, zeta, X, Xp3;

	n0 = n2/PI;
	zeta = 1-v2*v2/n2/n2;

	/*
	if (zeta < 1e-8) {
		fphi3 = 0;
	} else {
	*/
		X = get_X(n2,n3,v2,Ma,K,T_assoc,d);
		Xp3 = get_Xp3(n2,n3,v2,Ma,K,T_assoc,d);
		fphi3 = Ma*(log(X) - 0.5*X + 0.5)*(1/PI/d/d)*(0) + \
			Ma*zeta*n0*(1./X - 0.5)*Xp3;
	//}
	return fphi3;
}

double phip_assoc_v2(double n2, double n3, double v2, double Ma, double K, double T_assoc, double d, int FHSFlag) {
	double fphiv2 = 0.0;
	double n0, zeta, X, Xpv2;

	n0 = n2/PI;
	zeta = 1-v2*v2/n2/n2;
	
	if (zeta < 1e-8) {
		fphiv2 = 0;
	} else {
		X = get_X(n2,n3,v2,Ma,K,T_assoc,d);
		Xpv2 = get_Xpv2(n2,n3,v2,Ma,K,T_assoc,d);
		fphiv2 = Ma*(log(X) - 0.5*X + 0.5)*(-2/PI/d/d)*(v2/n2) + \
			Ma*zeta*n0*(1./X - 0.5)*Xpv2;
	}
	return fphiv2;
}

double phi_assoc(double n2, double n3, double v2, double Ma, double K, double T_assoc, double d, int FHSFlag) {
	double n0, zeta, X;
	n0 = n2/PI;
	zeta = 1-v2*v2/n2/n2;
	
	if (zeta==0) {
		return 0;
	} else {
		X = get_X(n2,n3,v2,Ma,K,T_assoc,d);
		return Ma*n0*zeta*(log(X)-X*0.5+0.5);
	}
}

/*********************************************/
/***************************************/
void c1AssocSpherical(double *c1AssocVec, double *rVec, double *denVec, double R, double Ma,double K,double T_assoc,double d,double dr, double rc, int FHSFlag, int wallFlag) {
	int i, j, range;
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
			n3Func[j] = (rp*rp)*den;
		}
		n3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		n3Vec[i] += n3*4*PI;
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

		/*
		if (r >= 0.5) {
			rpMin = (r-0.5); 
		} else {
			rpMin = (0.5-r); 
		}
		*/
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
			rp = rVec[rpMinIndex+j];
			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];
			n2Func[j] = rp*phip_assoc_n2(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);
			n3Func[j] = rp*((0.25)-((rp-r)*(rp-r)))*phip_assoc_n3(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);
			v2Func[j] = rp*(0.25 - (r-rp)*(r-rp))*phip_assoc_v2(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);

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

			v2Func[j] = rp*(r-rp)*phip_assoc_v2(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);
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

			n3Func[j] = rp*rp*phip_assoc_n3(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);

		}
		phi3 = trapz(n3Func,rpMin,rpMax,N);

		free(n3Func);

		phi3Vec[i] += phi3*4*PI;
	}
	phi2Vec[0] = phi2Vec[1];
	phi3Vec[0] = phi3Vec[1];
	phiv2Vec[0] = phiv2Vec[1];


	for (i=0; i<nGrids; i++) {
		c1AssocVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
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



void c1AssocSphericalSolute(double *c1AssocVec, double *rVec, double *denVec, double R,double RSolute, double Ma,double K,double T_assoc,double d,double dr, double rc, int FHSFlag, int wallFlag) {
	int i, j, range;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *phi2Vec, *phi3Vec, *phiv2Vec;
	double *n2Func, *n3Func, *v2Func;



	nGrids = (round)((R-RSolute)/dr)+1;
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
		rpMin = r - 0.5; 
		rpMax = r + 0.5;
		if (rpMin < RSolute) {
			rpMin = RSolute;
		}
		if (rpMax > R) {
			rpMax = R; 
		}

		rpMinIndex = i-(round)(0.5/dr);
		rpMaxIndex = i+(round)(0.5/dr);
		if (rpMinIndex < 0) {
			rpMinIndex = 0;
		}
		if (rpMaxIndex > nGrids) {
			rpMaxIndex = nGrids;
		}

		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n2Func[j] = den*rp;
			n3Func[j] = den*rp*((0.25)-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*((0.25)+r*r-rp*rp);
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
	n2Vec[0] = n2Vec[1];
	n3Vec[0] = n3Vec[1];
	v2Vec[0] = v2Vec[1];




	// phi2Vec time
	for (i=0; i<nGrids; i++) {
		phi2Vec[i] = 0;
		phi3Vec[i] = 0;
		phiv2Vec[i] = 0;
	}

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

		rpMinIndex = i-(round)(0.5/dr);
		rpMaxIndex = i+(round)(0.5/dr);
		if (rpMinIndex < 0) {
			rpMinIndex = 0;
		}
		if (rpMaxIndex > nGrids) {
			rpMaxIndex = nGrids;
		}

		N = (rpMaxIndex-rpMinIndex)+1;

		n2Func = (double*)malloc(N*sizeof(double));
		n3Func = (double*)malloc(N*sizeof(double));
		v2Func = (double*)malloc(N*sizeof(double));

		for (j=0; j<N; j++) {
			rp = rVec[rpMinIndex+j];
			n2 = n2Vec[rpMinIndex+j];
			n3 = n3Vec[rpMinIndex+j];
			v2 = v2Vec[rpMinIndex+j];
			n2Func[j] = rp*phip_assoc_n2(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);
			n3Func[j] = rp*((0.25)-((rp-r)*(rp-r)))*phip_assoc_n3(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);
			v2Func[j] = rp*(r*r-rp*rp+(0.25))*phip_assoc_v2(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);

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
	
	phi2Vec[0] = phi2Vec[1];
	phi3Vec[0] = phi3Vec[1];
	phiv2Vec[0] = phiv2Vec[1];


	// n2Vec, n3Vec, v2Vec calibration for infinite system
	for (i=nGrids-(int)(rc*200);i<nGrids;i++) {
		phi2Vec[i] = phi2Vec[nGrids-(int)(rc*200+1)];
		phi3Vec[i] = phi3Vec[nGrids-(int)(rc*200+1)];
		phiv2Vec[i] = phiv2Vec[nGrids-(int)(rc*200+1)];
	}

	for (i=0; i<nGrids; i++) { 
		c1AssocVec[i] = phi2Vec[i]+phi3Vec[i]-phiv2Vec[i];
	}

	free(n2Vec);
	free(n3Vec);
	free(v2Vec);
	free(phi2Vec);
	free(phi3Vec);
	free(phiv2Vec);
}



/************************/
void FAssocSpherical(double *FAssocVec, double *rVec, double *denVec, double R, double Ma, double K, double T_assoc, double d, double dr, double rc, int FHSFlag) {

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

	n2Vec = dVector(nGrids);
	n3Vec = dVector(nGrids);
	v2Vec = dVector(nGrids);
	FFunc = dVector(nGrids);


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

		n2Func = dVector(N);
		n3Func = dVector(N);
		v2Func = dVector(N);
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n2Func[j] = den*rp;
			n3Func[j] = den*rp*((0.25*(d*d))-((r-rp)*(r-rp)));
			v2Func[j] = den*rp*(0.25 - (r-rp)*(r-rp));
		}

		n2 = trapz(n2Func,rpMin,rpMax,N);
		n3 = trapz(n3Func,rpMin,rpMax,N);
		v2 = trapz(v2Func,rpMin,rpMax,N);

		free_dVector(n2Func);
		free_dVector(n3Func);
		free_dVector(v2Func);

		n2Vec[i] += n2*(PI*d/r);
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
	for (i=0; i<round(d*0.5/dr); i++) {
		r = rVec[i];
		rpMin = 0.;
		rpMax = d*0.5 - r;
		rpMinIndex = round(rpMin/dr);
		rpMaxIndex = round(rpMax/dr);
		N = (rpMaxIndex-rpMinIndex)+1;

		n3Func = dVector(N);
		for (j=0; j<N; j++) {
			den = denVec[rpMinIndex+j];
			rp = rVec[rpMinIndex+j];
			n3Func[j] = (rp*rp)*den;
		}
		n3 = trapz(n3Func,rpMin,rpMax,N);

		free_dVector(n3Func);

		n3Vec[i] += n3*4*PI;
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
			FAssocVec[i] = phi_assoc(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);
		} else {
			FAssocVec[i] = 0;
		}


	}




}
/////////////////////////////


void FAssocSphericalSolute(double *FAssocVec, double *rVec, double *denVec, double R, double RSolute,double Ma, double K, double T_assoc, double d, double dr, double rc, int FHSFlag) {

	int i, j;
	int rpMaxIndex, rpMinIndex, N, nGrids;
	double r, rp, rpMax, rpMin;
	double n2, n3, v2, phi2, phi3, phiv2;
	double den;
	double *n2Vec, *n3Vec, *v2Vec;
	double *n2Func, *n3Func, *v2Func;


	nGrids = (round)((R-RSolute)/dr)+1;

	n2Vec = dVector(nGrids);
	n3Vec = dVector(nGrids);
	v2Vec = dVector(nGrids);


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

		rpMinIndex = i-(round)(0.5/dr);
		rpMaxIndex = i+(round)(0.5/dr);
		if (rpMinIndex < 0) {
			rpMinIndex = 0;
		}
		if (rpMaxIndex > nGrids) {
			rpMaxIndex = nGrids;
		}

		N = (rpMaxIndex-rpMinIndex)+1;


		n2Func = dVector(N);
		n3Func = dVector(N);
		v2Func = dVector(N);
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

		free_dVector(n2Func);
		free_dVector(n3Func);
		free_dVector(v2Func);

		n2Vec[i] += n2*(PI*d/r);
		n3Vec[i] += n3*(PI/r);
		v2Vec[i] += v2*(PI/(r*r));
	}

	n2Vec[0] = n2Vec[1];
	n3Vec[0] = n3Vec[1];
	v2Vec[0] = v2Vec[1];

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
			FAssocVec[i] = phi_assoc(n2,n3,v2,Ma,K,T_assoc,d,FHSFlag);
		} else {
			FAssocVec[i] = 0;
		}
	}

	free_dVector(n2Vec);
	free_dVector(n3Vec);
	free_dVector(v2Vec);
}


