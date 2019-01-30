#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <math.h>

#include "constant.h"
#include "program.h"
#include "input.h"
#include "ideal.h"
#include "hard.h"
#include "disp.h"
#include "math2.h"
#include "analysis.h"

double fp(double eta, double T, double rMin, double rc ) {
	double p;
	p = (6*eta/PI)*((1+eta+eta*eta-eta*eta*eta)/(pow((1-eta),3)) + (48*eta)/(T)*(-1/(9*pow(rc,9))+1/(9*pow(rMin,9))+1/(3*pow(rc,3))-1/(3*pow(rMin,3))) - 4/T*eta*(pow(rMin,3)));
	return p;
}

double fmu(double eta,double T,double rMin,double rc) {
	double mu;
	mu = log(eta*pow(1/T,1.5)) + (eta*(8-9*eta+3*eta*eta))/(pow(1-eta,3)) + (96*eta)/(T)*(-1/(9*pow(rc,9))+1/(9*pow(rMin,9))+1/(3*pow(rc,3))-1/(3*pow(rMin,3))) - 8/T*eta*(pow(rMin,3));
	return mu;
}

double fmuHS(double eta,double T,double rMin,double rc) {
	double mu;
	mu = log(eta*pow(1/T,1.5)) + (eta*(8-9*eta+3*eta*eta))/(pow(1-eta,3));
	return mu;
}

int main(int argc, char *argv[])
{
	double rhob, rhob_v, eta, mu, muHS, muDisp, muBulk;
	double R, Ri, RSolute, T, d, r, dr, q, tol, tol2;
	unsigned int nGrids, nGrids2, nRGrids, again, step, s;
	unsigned int i, j, k, N;
	unsigned int Ntot;
	int tolFlag;
	double ff;
	double den, rp;
	double sig, eps, rc;
	double ref1 = 0, ref2 = 0, ref3 = 0;
	double radius;
	double d1 = 0, d2 = 0, d3 = 0;

	double F, F1, F2, F3, F4, Fs, FBulk;

	FILE *f;

	clock_t before;
	before = clock();

	char *fileName = argv[1];
	char *oFileName = argv[2];
	INPUT input;

	inputInitialize(&input);
	readInput(fileName,&input);
	
	R = input.size;
	Ri = input.Ri;
	Ntot = input.N;
	rhob = input.den;
	rhob_v = input.den_v;

	step = input.step;


	sig = input.sig;
	T = input.T;
	rc = input.rCut;

	eta = (PI*rhob)/6;
	muHS = eta*(8-9*eta+3*eta*eta)/(pow(1-eta,3.));
	double rMin = pow(2,1./6.)*sig;
	if (input.FDispFlag==0) {
		muDisp = 0;
		d = 1;
	} else {
		//d = (1+0.2977*T)*sig/(1+0.33163*T+0.00104771*T*T);
		d = 1;
		eta = (PI*rhob/6)*(pow(d,3.));
	//muDisp = (16*PI*rhob*(1/T)/9)*(-(pow(sig,12.)/pow(rc,9.))+(3*pow(sig,6.)/pow(rc,3.))-(2*pow(sig,3.)));

		muDisp = (96*eta/T)*(-(pow(sig,12.)/(9*pow(rc,9.)))+(pow(sig,6.)/(3*pow(rc,3.)))+(pow(sig,12.)/(9*pow(rMin,9.)))-(pow(sig,6.)/(3*pow(rMin,3.)))) - 8*eta*pow(rMin,3);

	}

	mu = muHS + muDisp;
	mu = input.mu;
	printf("Mu = %f\n",mu);

	// Initial setting , all length unit is dma
	R = input.size;		// R is size of slab 
	RSolute = 0;
	if (input.systemFlag == 3) {
		RSolute = R;
		R += 20.;
	}
	dr = input.dr;	// dr is grid step, here, 0.01*d
	//nRGrids = round(R/dr)+1;
	//nGrids = nRGrids+round(5*d/dr); // external vapor region length = 5*d
	nGrids = round(R/dr)+1;
	nGrids2 = round(1/0.0001)+1;

	q = input.q;

	tol = input.tol;
	tol2 = 0;
	again = 0;
	s = 0;

	double *rVec;
	double *denOutVec, *denInVec, *denVec, *dVec;
	double *c1HSVec;
	double *c1DispVec;
	double *cAttVec;
	double *diffVec;
	double *tolVec;
	double *funcVec;

	rVec = (double *)malloc(nGrids*sizeof(double));
	denVec = (double*)malloc(nGrids*sizeof(double));
	denInVec = (double*)malloc(nGrids*sizeof(double));
	denOutVec = (double*)malloc(nGrids*sizeof(double));
	c1HSVec = (double*)malloc(nGrids*sizeof(double));
	c1DispVec = (double*)malloc(nGrids*sizeof(double));
	tolVec = (double*)malloc(nGrids*sizeof(double));
	funcVec = (double*)malloc(nGrids*sizeof(double));

	dVec = (double*)malloc(nGrids2*sizeof(double));

	for (i=0; i<nGrids; i++) { 
		rVec[i] = 0;
		denVec[i] = 0;
		denInVec[i] = 0; 
		denOutVec[i] = 0;
		c1HSVec[i] = 0;
		c1DispVec[i] = 0;
		tolVec[i] = 0;
	}

	for (i=0; i<nGrids2; i++) {
		dVec[i] = 0.0001*i;
	}

	// DFT start
	//for (i=0; i<(nGrids-(int)(d/2/dr)); i++) {

	if (input.systemFlag == 2) {
		if (input.ensembleFlag == 1) {
			radius = pow(((3*Ntot/4/PI)-(pow(R-d*0.5,3)*rhob_v))/(rhob-rhob_v),1/3);
			for (i=0; i<(int)(radius/dr); i++) {
				denInVec[i] += rhob;
			}
			for (i=(int)(radius/dr)+1; i<nGrids; i++) {
				denInVec[i] += rhob_v;
			}
		}
		else if (input.ensembleFlag == 2) {
			for (i=0; i<(int)round(Ri/dr); i++) {
				denInVec[i] += rhob;
			} for (i=(int)round(Ri/dr); i<nGrids;i++) {
				denInVec[i] += rhob_v;
			}
		}
	}


	for (i=0; i<nGrids; i++) {
		rVec[i] += dr*i;
	}

	// Main DFT loop
	char *omFileName = argv[3];
	FILE *fm;
	fm = fopen(omFileName,"w");

	while (again < 1) {
		s += 1;

		for (i=0; i<nGrids; i++) {
			denVec[i] = denInVec[i];
			denOutVec[i] = denInVec[i];
			denInVec[i] = 0;
		}

		c1DispSpherical(c1DispVec,rVec,denVec,R,dr,sig,T,rc);

		// From Here
		for (i=0; i<nGrids; i++) {
			r = rVec[i];


			ref1 = 0;
			ref2 = 0;
			while (ref1*ref2 >= 0) {
				d1 = (double)rand()/RAND_MAX;
				d2 = (double)rand()/RAND_MAX;

				if (d1 == d2) {
					d1 = (double)rand()/RAND_MAX;
					d2 = (double)rand()/RAND_MAX;
				} else if (d1>d2) {
					d3 = d2;
					d2 = d1;
				    	d1 = d3;	
				}
				ref1 = fmuHS(PI/6*d1,T,rMin,rc) + c1DispVec[i] - mu;
				ref2 = fmuHS(PI/6*d2,T,rMin,rc) + c1DispVec[i] - mu;
				tol2 = fabs(ref1-ref2);

			}


			while (tol2>1E-5) {
				d3 = (d1+d2)/2;
				ref3 = fmuHS(PI/6*d3,T,rMin,rc) + c1DispVec[i] - mu;

				if (ref1*ref3 < 0) {
					d2 = d3;
					ref2 = ref3;
				} else if (ref3*ref2 < 0) {
					d1 = d3;
					ref1 = ref3;
				}
				tol2 = fabs(ref1-ref2);
			}

			denOutVec[i] = (d1+d2)/2;
			printf("%f\n",denOutVec[i]);

		}

		tolFlag = checkTol(tolVec, denVec, denOutVec, nGrids, tol);

		if (tolFlag == 1 || step == s){
			char *oFileName = argv[2];
			f = fopen(oFileName, "w");
			if(f==NULL) { printf("File error!"); }
			else {
				for (i=0; i<nGrids; i++) {
					fprintf(f, "%d \t %f \n", i, denOutVec[i]);
				}
			}
			fclose(f);

			again = 1;
		}
		else {
			updateDen(denVec, denInVec, denOutVec, nGrids, q);

			F1 = FIdealSpherical(denOutVec, R, d, dr, 6/PI*pow(T,1.5));
			F2 = FHSSpherical(rVec, denOutVec, R, d, dr, input.FHSFlag);
			F3 = FDispSpherical(c1DispVec, rVec, denOutVec, R, dr, sig, T, rc);
			double *FFunc;
			FFunc = (double*)malloc(nGrids*sizeof(double));
			for (i=0; i<nGrids; i++) {
				FFunc[i] = 0;
			}
			for (i=0; i<nGrids; i++) {
				r = rVec[i];
				FFunc[i] = 4*PI*r*r*denOutVec[i]*(-mu);
			}
			F4 = trapz(FFunc,0,R,nGrids);
			free(FFunc);

			printf("%f\t%f\t%f\t%f\t%f\t%f\n",denOutVec[0], F1, F2, F3, F4, F1+F2+F3+F4);

			fprintf(fm, "%d\t%lf\n",s,F1+F2+F3+F4);
		}
	}
	fclose(fm);
	free(dVec);
	free(denVec);


	/************** Analysis section *******************/
	double rhol, rhog, Re, Rs, Rb = 0;
	double tension, tension2;
	double  V, A;
	double Vl, Vv;
	double mu_id, mu_ext;
	double dv, dl, p, pv, pl, eta_l, eta_v, eta1, eta2, mu1, mu2;
	double ohm;
	int len, index_l, index_v;

	nGrids = nGrids - (int)round(10*d/dr);
	R = R-10;
	denVec = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		denVec[i] = denOutVec[i];
	}

	len = (int)((1.0-0.001)/0.0005);
	dVec = (double*)malloc(len*sizeof(double));
	for (i=0; i<len; i++) {
		dVec[i] = 0.001+0.0005*i;
	}

//	mu_id = log(eta);
//	mu_ext = log((4*PI/3*R*R*R)/(4*PI*ff));
	//mu = log((Ntot*PI/6*pow(1/T,1.5))/(4*PI*ff));
	printf("mu = %f\n", mu);

	for (i=0; i<len; i++) {
		d1 = dVec[i];
		d2 = dVec[i+1];

		eta1 = PI/6*d1;
		eta2 = PI/6*d2;

		mu1 = fmu(eta1,T,rMin,rc);
		mu2 = fmu(eta2,T,rMin,rc);

		if (mu2-mu1 < 0) { 
			index_v = i;
			break;
		}
	}


	for (i=len-1; i>0; i--) {
		d1 = dVec[i];
		d2 = dVec[i-1];

		eta1 = PI/6*d1;
		eta2 = PI/6*d2;

		mu1 = fmu(eta1,T,rMin,rc);
		mu2 = fmu(eta2,T,rMin,rc);

		if (mu1-mu2 < 0) {
			index_l = i;
			break;
		}
	}

	for (i=0; i<index_v; i++) {
		d1 = dVec[i];
		d2 = dVec[i+1];
		eta1 = PI/6*d1;
		eta2 = PI/6*d2;
		mu1 = fmu(eta1,T,rMin,rc);
		mu2 = fmu(eta2,T,rMin,rc);
		if (mu1<mu && mu2>mu) {
			dv = (d1+d2)/2;
			break;
		}
	}

	for (i=len-1; i>index_l; i--) {
		d1 = dVec[i];
		d2 = dVec[i-1];
		eta1 = PI/6*d1;
		eta2 = PI/6*d2;
		mu1 = fmu(eta1,T,rMin,rc);
		mu2 = fmu(eta2,T,rMin,rc);
		if (mu1>mu && mu2<mu) {
			dl = (d1+d2)/2;
			break;
		}
	}
	free(dVec);
	eta_l = PI/6*dl;
	eta_v = PI/6*dv;
	pl = fp(eta_l,T,rMin,rc);
	pv = fp(eta_v,T,rMin,rc); // Now we know dl, dv, pl, pv


	double *fVec;
	fVec = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		den = denVec[i];
		fVec[i] = 3/(dl-dv)*r*r*(den-dv);
	}
	free(fVec);

	Re = trapz(fVec,0,rVec[nGrids-1],nGrids);
	Re = pow(Re,1./3.); // Now we know Re

/*
	for (i=0; i<nGrids; i++) {
		denInVec[i] = 0;
	}
	for (i=0; i<nGrids; i++) {
		denInVec[i] = Ntot/(4*PI/3*pow(R,3));
	}
	*/
	F1 = FIdealSpherical(denVec, R, d, dr, 6/PI*pow(T,1.5));
	F2 = FHSSpherical(rVec, denVec, R, d, dr, input.FHSFlag);
	F3 = FDispSpherical(c1DispVec, rVec, denVec, R, dr, sig, T, rc);
	double *FFunc;
	FFunc = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		FFunc[i] = 0;
	}
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		FFunc[i] = 4*PI*r*r*denVec[i]*(-mu);
	}
	F4 = trapz(FFunc,0,R,nGrids);
	free(FFunc);

	F = F1+F2+F3+F4;
/*
	double *FFunc;
	FFunc = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		FFunc[i] = 0;
	}
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		FFunc[i] = 4*PI*r*r*denVec[i]*(-mu);
	}
	F += trapz(FFunc,0,R,nGrids);
	*/
/*
	F1 = FIdealSpherical(denInVec, R, d, dr, 6/PI*pow(T,1.5));
	F2 = FHSSpherical(rVec, denInVec, R, d, dr, input.FHSFlag);
	c1DispSpherical(c1DispVec,rVec,denInVec,R,dr,sig,T,rc);
	F3 = FDispSpherical(c1DispVec, rVec, denInVec, R, dr, sig, T, rc);
	FBulk = F1+F2+F3;
	*/
	
	//printf("%f\t%f\t%f\t%f\n",F1, F2, F3, F1+F2+F3);
	A = 4*PI*Re*Re;
	V = 4*PI/3*R*R*R;

	//Rs = 3*(F+pv*V)/(2*PI*(pl-pv));
	//Rs = pow(Rs,1./3.);
	tension2 = (3*(F+pv*V)*(pl-pv)*(pl-pv))/(16*PI);
	tension2 = pow(tension2,1./3.);

	tension = (F+pv*V)/A+((pl-pv)*Re)/3;

	char *opFileName = strcat(oFileName,".energy");
	f = fopen(opFileName,"w");
	fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Re,dl,dv,pl,pv,tension);
	fclose(f);



/*
	p = (6*eta/PI)*((1+eta+eta*eta-eta*eta*eta)/(pow((1-eta),3)) + (48*eta)/(T)*(-1/(9*pow(rc,9))+1/(9*pow(rMin,9))+1/(3*pow(rc,3))-1/(3*pow(rMin,3))) - 4*eta*(pow(rMin,3)));
	//p = (6*eta/PI)*((1+eta+eta*eta-eta*eta*eta)/pow((1-eta),3));
	V = 4*PI/3*pow((R-0.5),3);
	A = 4*PI*pow(R-0.5,2);


	double *FFunc;
	FFunc = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		FFunc[i] = 0;
	}
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		FFunc[i] = 4*PI*r*r*denVec[i]*(-mu);
	}
	F4 = trapz(FFunc,0,R,nGrids);

	free(FFunc);
*/
	//F = F1+F2+F3+F4;
	/*
	F = F1+F2+F3+F4;
	tension = (F+p*V)/A;
	char *opFileName = strcat(oFileName,".energy");
	f = fopen(opFileName,"w");
	fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",F1,F2,F3,F4,F,p,V,A,tension);
	fclose(f);

	exit(1);
	*/

	/*
	double *FFunc;
	FFunc = (double*)malloc(nGrids*sizeof(double));
	for (i=0; i<nGrids; i++) {
		FFunc[i] = 0;
	}
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		FFunc[i] = 4*PI*r*r*denVec[i]*(-mu);
	}
	F += trapz(FFunc,0,R,nGrids);
	for (i=0; i<nGrids; i++) {
		r = rVec[i];
		FFunc[i] = 4*PI*r*r*denBulkVec[i]*(-mu);
	}
	FBulk += trapz(FFunc,0,R,nGrids);
	free(FFunc);

	printf("%lf\t%lf\n",F,FBulk);
	*/


	/*************************************************/
	free(rVec);
	free(denVec);
	free(denInVec);
	free(denOutVec);
	free(c1HSVec);
	free(c1DispVec);
	free(tolVec);
	free(funcVec);

	printf("Total time is %f s\n", ((double)clock()-before)/CLOCKS_PER_SEC);

}

