#include <stdlib.h>
#include <math.h>

#include "constant.h"
#include "math2.h"
#include "memory.h"

#define P4(x) ((x)*(x)*(x)*(x))
#define P6(x) ((x)*(x)*(x)*(x)*(x)*(x))

double uAtt(double r, double sig, double T, double rc) {
	double rpMin;
	//rpMin = pow(2,(1./6.))*sig;
	rpMin = 1.122462;
	double u = 0.;
	if (r<rpMin) {
		u = -1.*(1/T);
	}
	else if (r>=rpMin && r<rc) {
		u = 4*(1/T)*(pow(sig/r,12)-pow(sig/r,6));
		//u = 4*(1/T)*(pow(sig/r,12)-P6((sig/r)));
	}
	else {
		u = 0.;
	}
	return u;
}



void c1DispSpherical(double *c1DispVec, double *rVec, double *denVec, double R,  double T, double rc) {
	int i, j, k;
	int jMin, jMax;
	int rpMaxIndex, rpMinIndex, nGrids;
	double r, rp, rpMax;
	double cAtt, dist;
	double den, denp, u, xi, c;
	double *cAttVec;
	double dr;
	double rpMin = pow(2,(1./6.));
	double ddr;
	int range;
	double inv_rc4,inv_rc10,inv_rpMin4,inv_rpMin10;


	dr = rVec[1]-rVec[0];

	range = round(rc/dr);
	inv_rc4 = 1/(rc*rc*rc*rc);
	inv_rc10 = 1/(rc*rc*rc*rc*rc*rc*rc*rc*rc*rc);
	inv_rpMin4 = 1/(rpMin*rpMin*rpMin*rpMin);
	inv_rpMin10 = 1/(rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin);

	double d = 1;

	nGrids = (int)(R/dr)+1;


	// Why error?
	cAttVec = dVector(nGrids);

	
	for (i=1; i<nGrids; i++) {
		r = rVec[i];
		cAttVec[0] = 0;

		//for (j=0; j<nGrids; j++) {
		//	cAttVec[j] = 0;
		//}
	//	for (j=1; j<nGrids-round(sig*0.5/dr); j++) {
		jMin = i-range;
		if (jMin < 1) {
			jMin = 1;
		}
		jMax = i+range;
		if (jMax > nGrids) {
			jMax = nGrids;
		}

		//for (j=1; j<nGrids; j++) {
		for (j=jMin; j<jMax; j++) {
			rp = rVec[j];
			denp = denVec[j];
			dist = abs(j-i)*dr;
			
			ddr = (i+j)*dr;

			if (dist > rc) {
				cAttVec[j] = 0;
			} else if (dist > rpMin) {
				if (ddr > rc) {
					cAttVec[j] = inv_rc4 - 0.4*inv_rc10 - (1)/(dist*dist*dist*dist) + 0.4/(dist*dist*dist*dist*dist*dist*dist*dist*dist*dist);
				} else {
					cAttVec[j] = (1)/(ddr*ddr*ddr*ddr) - 0.4/(ddr*ddr*ddr*ddr*ddr*ddr*ddr*ddr*ddr*ddr) - (1)/(dist*dist*dist*dist) + 0.4/(dist*dist*dist*dist*dist*dist*dist*dist*dist*dist);
				}
			} else {
				if (ddr > rc) {
					cAttVec[j] = inv_rc4 - 0.4*inv_rc10 - inv_rpMin4 + 0.4*inv_rpMin10 +0.5*(dist*dist-rpMin*rpMin);
				} else if (ddr > rpMin) {
					cAttVec[j] = (1)/(ddr*ddr*ddr*ddr) - 0.4/pow(ddr,10.) - inv_rpMin4 + 0.4*inv_rpMin10 +0.5*(dist*dist-rpMin*rpMin);
				} else {
					cAttVec[j] = 0.5*(dist*dist-ddr*ddr);
				}
			}
			
			cAttVec[j] *= rp*denp;

		}

		c1DispVec[i] = trapz(cAttVec,0.,R,nGrids)*2*PI/(r*T);
		//c1DispVec[i] *= (2*PI/(r*T));
	}
	c1DispVec[0] = c1DispVec[1];

	/*
	for (i=0; i<(int)(rc*200); i++) {
		c1DispVec[i] = c1DispVec[(int)(rc*200)+1];
	}
	for (i=nGrids-(int)(rc*200);i<nGrids;i++) {
		c1DispVec[i] = c1DispVec[nGrids-(int)(rc*200)-1];
	}
	*/

	free_dVector(cAttVec);


}



/*************************************************************/
void c1DispSphericalSolute(double *c1DispVec, double *rVec, double *denVec, double R,  double RSolute,double T, double rc) {
	int i, j, k;
	int jMin, jMax;
	int rpMaxIndex, rpMinIndex, nGrids;
	int N;
	double r, rp, rpMax;
	double cAtt, dist;
	double den, denp, u, xi, c;
	double *cAttVec;
	double dr;
	double rpMin = pow(2,(1./6.));
	double ddr;
	int range;
	double inv_rc4,inv_rc10,inv_rpMin4,inv_rpMin10;


	dr = rVec[1]-rVec[0];

	range = round(rc/dr);
	inv_rc4 = 1/(rc*rc*rc*rc);
	inv_rc10 = 1/(rc*rc*rc*rc*rc*rc*rc*rc*rc*rc);
	inv_rpMin4 = 1/(rpMin*rpMin*rpMin*rpMin);
	inv_rpMin10 = 1/(rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin*rpMin);

	double d = 1;

	nGrids = round((R-RSolute)/dr)+1;
	for (i=0; i<nGrids; i++) { 

		r = rVec[i];

		jMin = i-range;
		if (jMin < 1) {
			jMin = 1;
		}

		jMax = i+range;
		if (jMax > nGrids) {
			jMax = nGrids;
		}

		//for (j=1; j<nGrids; j++) {
		N = (jMax-jMin)+1;

		cAttVec = dVector(N);

		for (j=0; j<N; j++) {
			rp = RSolute + dr*(jMin+j);
			if (rp<RSolute+nGrids*dr) {
				denp = denVec[jMin+j];
			} else {
				denp = denVec[nGrids-1];
			}
			dist = abs(i-(jMin+j))*dr;
			
			ddr = RSolute + (i+(jMin+j))*dr;

			if (dist > rc) {
				cAttVec[j] = 0;
			} else if (dist > rpMin) {
				if (ddr > rc) {
					cAttVec[j] = inv_rc4 - 0.4*inv_rc10 - (1)/(dist*dist*dist*dist) + 0.4/(dist*dist*dist*dist*dist*dist*dist*dist*dist*dist);
				} else {
					cAttVec[j] = (1)/(ddr*ddr*ddr*ddr) - 0.4/(ddr*ddr*ddr*ddr*ddr*ddr*ddr*ddr*ddr*ddr) - (1)/(dist*dist*dist*dist) + 0.4/(dist*dist*dist*dist*dist*dist*dist*dist*dist*dist);
				}
			} else {
				if (ddr > rc) {
					cAttVec[j] = inv_rc4 - 0.4*inv_rc10 - inv_rpMin4 + 0.4*inv_rpMin10 +0.5*(dist*dist-rpMin*rpMin);
				} else if (ddr > rpMin) {
					cAttVec[j] = (1)/(ddr*ddr*ddr*ddr) - 0.4/pow(ddr,10.) - inv_rpMin4 + 0.4*inv_rpMin10 +0.5*(dist*dist-rpMin*rpMin);
				} else {
					cAttVec[j] = 0.5*(dist*dist-ddr*ddr);
				}
			}
			
			cAttVec[j] *= rp*denp;
		}

		c1DispVec[i] = trapz(cAttVec,jMin*dr,jMax*dr,N)*2*PI/(r*T);


		free_dVector(cAttVec);
	}


	for (i=nGrids-(int)(rc*200);i<nGrids;i++) {
		c1DispVec[i] = c1DispVec[nGrids-(int)(rc*200)-1];
	}

}




/**************************************************************/
void c1ChainDispFlat(double **c1DispMat,double *rVec,double **denMat,double Z, double M, double T,double rc) {
	int i, j, k, l;
	double cAtt;
	double den;
	double z1, z2, r12;
	int nGrids,n1Grids, n2Grids;
	int denIndex;
	double theta;
	double dTheta = PI/200;
	double dr = 0.01;
	double rMin = pow(2,1./6.);
	double *c1Vec, *c2Vec;
	double *lut_UrLJ;
	double *lut_Ure;
	double *lut_sin, *lut_cos;
	double dz;

	for (i=0; i<M; i++) {
		printf("M = %d\n",i);
		dz = rVec[1]-rVec[0];

		nGrids = (int)(Z/dz)+1;
		n1Grids = (int)(rc/dr)+1;
		n2Grids = (int)(PI/dTheta)+2;

		c1Vec = (double*)malloc(n1Grids*sizeof(double));
		c2Vec = (double*)malloc(n2Grids*sizeof(double));
		for (j=0; j<nGrids; j++) {
			c1DispMat[i][j] = 0.;
		}
		for (j=0; j<n1Grids; j++) {
			c1Vec[j] = 0.;
		}
		for (j=0; j<n2Grids; j++) {
			c2Vec[j] = 0.;
		}

		lut_UrLJ = (double*)malloc(n1Grids*sizeof(double));
		lut_Ure = (double*)malloc(n1Grids*sizeof(double));
		for (j=0; j<n1Grids; j++) {
			r12 = j*dr;
			lut_UrLJ[j] = (4/T)*(pow(1/r12,12)-pow(1/r12,6));
			lut_Ure[j] = (-1/T);
		}

		lut_sin = (double*)malloc(n2Grids*sizeof(double));
		lut_cos = (double*)malloc(n2Grids*sizeof(double));
		for (j=0; j<n2Grids; j++) {
			theta = dTheta*j;
			lut_sin[j] = sin(theta);
			lut_cos[j] = cos(theta);
		}




		for (j=0; j<nGrids; j++) {

			z1 = rVec[j];

			for (k=0; k<n1Grids; k++) {
				r12 = k*dr;

				for (l=0; l<n2Grids; l++) {
					denIndex = (int)(floor)(((z1+r12*lut_cos[l])/dz)+0.5);


					if (denIndex<0) {
						den = 0;
					} else if (denIndex>=nGrids) {
						den = 0;
					} else {
						den = denMat[i][denIndex];
					}

					if (r12 > rc) {
						c2Vec[l] = 0;
					} else if (r12 <= rc && r12 > rMin) {
						c2Vec[l] = 2*PI*r12*r12*lut_sin[l]*den*lut_UrLJ[k];
					} else {
						c2Vec[l] = 2*PI*r12*r12*lut_sin[l]*den*lut_Ure[k];
					}
				}
				
				c1Vec[k] = trapz(c2Vec,0,(n2Grids-1)*dTheta,n2Grids);
			}
			c1DispMat[i][j] = trapz(c1Vec,0,(n1Grids-1)*dr,n1Grids);
		}


		free(c1Vec);
		free(c2Vec);
		free(lut_UrLJ);
		free(lut_Ure);
		free(lut_sin);
		free(lut_cos);
	}
}




void FDispSpherical(double *FDispVec, double *c1DispVec, double *rVec, double *denVec, double R, double dr, double sig, double T, double rc) {
	int i;
	int nGrids = (int)round(R/dr)+1;
	double r, den;
	
	for (i=0; i<nGrids; i++) {
		den = denVec[i];
		FDispVec[i] = 0.5*den*c1DispVec[i];
	}


}

/*************************************************************/
/*******************************/
void FDispSphericalSolute(double *FDispVec, double *c1DispVec, double *rVec, double *denVec, double R, double RSolute,double dr, double sig, double T, double rc) {
	int i;
	int nGrids = round((R-RSolute)/dr)+1;
	double r, den;
	

	nGrids = round((R-RSolute)/dr)+1;

	for (i=0; i<nGrids; i++) { 
		den = denVec[i];
		FDispVec[i] = 0.5*den*c1DispVec[i];
	}
}

