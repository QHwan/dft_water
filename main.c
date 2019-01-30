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
#include "assoc.h"
#include "chain.h"
#include "math2.h"
#include "analysis.h"
#include "memory.h"


int main(int argc, char *argv[])
{
	double min;
	double rhob, eta, mu, mu_ext, muIdeal, muHS, muAssoc, muDisp, muChain;
	double pIdeal, pHS, pAssoc, pDisp;
	double dl, dv, pl, pv;
	double R, RSolute, Ri, dR, T, d, r, dr, rp, q, tol, tolval, tolval2;
	double rpMax, rpMin;
	int rpMaxIndex, rpMinIndex;
	double A, V, p;
	unsigned int nGrids, again, step, step_DFT, s;
	unsigned int nImages;
	unsigned int i, j, k, N, len;
	double Ntot;
	int tolFlag;
	double den;
	double sig, eps, rc, rMin;
	double tension, tension_Re;
	double buf;
	double dist;

	double F, F_id, F_hs, F_assoc, F_disp, F_chain, F_mu;
	double Re, Rs;
	double R1, R2, R3;
	double N1, N2, N3, dN, NMax, NMin, NRef;
	double M, Ma, K, T_assoc;

	double a;

	double wallT;

	FILE *f;

	clock_t before;
	before = clock();

	char *fileName = argv[1];
	char *oFileName = argv[2];
	char *o2FileName = argv[3];
	char *testFileName;


	int evol1, evol2, evol3;
	int try;

	// Read inputfile
	INPUT input;

	inputInitialize(&input);
	readInput(fileName,&input);
	
	Ntot = input.N;
	M = input.M;
	Ma = input.Ma;
	dl = input.den;
	dv = input.den_v;

	printf("%f\n",Ma);

	pl = input.pl;
	pv = input.pv;


	step = input.step;


	sig = input.sig;
	T = input.T;
	T_assoc = input.T_assoc;
	//T_assoc = T_assoc*T;
	rc = input.rCut;

	rMin = pow(2,1./6.)*sig;

	R = input.R;		// R is size of slab 
	Ri = input.Ri;
	RSolute = input.RSolute;
	dr = input.dr;	// dr is grid step, here, 0.01*d

	d = 1.0;
	K = 1.4849*1e-4*d*d*d;

	mu = input.mu;


	nGrids = round((R-RSolute)/dr)+1;


	q = input.q;

	tol = input.tol;
	again = 0;
	s = 0;

	Ntot += 4*PI*R*R*R*dv/3;

	// mu and p
	a = get_a(T,rc);	
	muIdeal = fmu_id(dl);
	muHS = fmu_hs(dl);
	if (input.FDispFlag!=0) {
		muDisp = fmu_disp(dl,a);
	} else {
		muDisp = 0;
	}

	if (Ma != 0) {
		muAssoc = fmu_assoc(dl,Ma,K,T_assoc);
	} else {
		muAssoc = 0;
	}
	mu_ext = muHS + muDisp + muAssoc;
	mu = muIdeal + muHS + muDisp + muAssoc;

	printf("%f\t%f\t%f\t%f\n",muIdeal,muHS,muDisp,muAssoc);
	printf("mu = %f\n",mu);


	// allocate vector and matrix
	double *rVec, *rEVec;
	double *denOutVec, *denInVec, *denVec, *denEVec;
	double *c1HSVec, *c1AssocVec;
	double *c1DispVec, *c1DispEVec;
	double *FIdealVec, *FHSVec, *FDispVec, *FAssocVec, *FVec;
	double *cAttVec;
	double *diffVec;
	double *tolVec;
	double *funcVec;
	double *wallVec;
	double *sVec, *snewVec;

	double **denOutMat, **denInMat, **denMat;
	double **c1DispMat;
	double **c1HSMat;
	double **c1ChainMat;
	double **GMat;

	double func;



	rVec = dVector(nGrids);
	denVec = dVector(nGrids);
	denInVec = dVector(nGrids);
	denOutVec = dVector(nGrids);
	c1HSVec = dVector(nGrids);
	c1AssocVec = dVector(nGrids);
	c1DispVec = dVector(nGrids);
	tolVec = dVector(nGrids);
	wallVec = dVector(nGrids);
	funcVec = dVector(nGrids);

	FIdealVec = dVector(nGrids);
	FHSVec = dVector(nGrids);
	FDispVec = dVector(nGrids);
	FAssocVec = dVector(nGrids);
	FVec = dVector(nGrids);

	if (M > 1) {
		denOutMat = dMatrix(M,nGrids);
		denInMat = dMatrix(M,nGrids);
		denMat = dMatrix(M,nGrids);
		c1HSMat = dMatrix(M,nGrids);
		c1DispMat = dMatrix(M,nGrids);
		c1ChainMat = dMatrix(M,nGrids);
		GMat = dMatrix(M,nGrids);
	}




	for (i=0; i<nGrids; i++) {
		rVec[i] = dr*i;
	}

	for (i=0; i<(int)round(Ri/dr); i++) {
		denInVec[i] = dl;
	}
	for (i=(int)round(Ri/dr); i<nGrids; i++) {
		denInVec[i] = dv;
	}

	Ntot = (4.*PI/3.)*Ri*Ri*Ri*dl + (4.*PI/3.)*(R*R*R-Ri*Ri*Ri)*dv;

	
	again = 0;
	s = 0;

	printf("dl, dv = %f\t%f\n",dl,dv);





	// Iterate using picard iteration
	while (again < 1) {
		s++;

		// update density profile
		for (i=0; i<nGrids; i++) {
			denVec[i] = denInVec[i];
			denOutVec[i] = denInVec[i];
			denInVec[i] = 0;
		}

		if (M > 1) {
			for (i=0; i<M; i++) {
				for (j=0; j<nGrids; j++) {
					denMat[i][j] = denInMat[i][j];
					denOutMat[i][j] = denInMat[i][j];
					denInMat[i][j] = 0;
				}
			}
		}

		// calculate c1 from density profile
		if (M == 1) {
			c1HSSpherical(c1HSVec,rVec,denVec,R,rc); 

			if (input.FDispFlag!=0) {
				c1DispSpherical(c1DispVec, rVec, denVec,R,T, rc);
			}

			if (Ma != 0) {
				c1AssocSpherical(c1AssocVec,rVec,denVec,R,Ma,K,T_assoc, d,dr, rc, input.FHSFlag, input.wallFlag);
			}

			
			/*
			for (i=nGrids-(int)(rc*200); i<nGrids; i++) {
				c1HSVec[i] = c1HSVec[nGrids-(int)(rc*200)-1];
				c1DispVec[i] = c1DispVec[nGrids-(int)(rc*200)-1];
				c1AssocVec[i] = c1AssocVec[nGrids-(int)(rc*200)-1];
			}
			*/
			


			for (i=0; i<nGrids; i++) {
				funcVec[i] = 4.*PI*rVec[i]*rVec[i]*exp(-c1HSVec[i]-c1DispVec[i]-c1AssocVec[i]);
			}
			func = trapz(funcVec,0,R,nGrids);

			for (i=0; i<nGrids; i++) {
				denOutVec[i] = (Ntot/func)*exp(-c1HSVec[i]-c1DispVec[i]-c1AssocVec[i]);
			}

			/*
			for (i=0; i<(int)(rc*200); i++) {
				denOutVec[i] = dl;
			}

			for (i=nGrids-(int)(rc*200); i<nGrids;i++) {
				denOutVec[i] = denOutVec[nGrids-(int)(rc*200)-1];
			}
			*/

			
			tolFlag = checkTol(tolVec, denVec, denOutVec, nGrids, tol);


			if (tolFlag == 1 || step == s){

				/*
				for (i=0; i<(int)(rc*200);i++) {
					denOutVec[i] = denOutVec[(int)(rc*200)+1];
				}
*/

				/*
				for (i=nGrids-(int)(rc*200); i<nGrids;i++) {
					denOutVec[i] = denOutVec[nGrids-(int)(rc*200)-1];
				}
				*/

				f = fopen(oFileName, "w");
				if(f==NULL) { printf("File error!"); }
				else {
					for (i=0; i<nGrids; i++) {
						fprintf(f, "%f \t %f \n", i*dr, denOutVec[i]);
					}
				}
				fclose(f);

				again = 1;
			}
			else {
				for (i=0; i<nGrids; i++) {
					//funcVec[i] = 4*PI*rVec[i]*rVec[i]*(denOutVec[i]-dv);
					funcVec[i] = 4*PI*rVec[i]*rVec[i]*(denOutVec[i]);
				}
				Ntot = trapz(funcVec,0,R,nGrids);
				printf("%d\t%f\n",s,Ntot);

				

				updateDen(denVec, denInVec, denOutVec, nGrids, q);
			}
		} else {
			continue;
		}
	}



	/************** Analysis section *******************/
	/*************************************************/
	printf("Analysis...\n");


	if (M == 1) {

		FIdealSpherical(FIdealVec, denOutVec, R, d, dr);

		c1HSSpherical(c1HSVec,rVec,denOutVec,R,rc); 
		FHSSpherical(FHSVec, rVec, denOutVec, R, d,rc, input.FHSFlag); 

		if (Ma != 0) {
			c1AssocSpherical(c1AssocVec,rVec,denOutVec,R,Ma,K,T_assoc, d,dr, rc, input.FHSFlag, input.wallFlag);
			FAssocSpherical(FAssocVec, rVec, denOutVec, R,Ma, K, T_assoc, d, dr, rc, input.FHSFlag); 
		} 

		if (input.FDispFlag!=0) {
			c1DispSpherical(c1DispVec, rVec, denOutVec,R,T, rc);
			/*
			for (i=nGrids-(int)(rc*200); i<nGrids; i++) {
				c1DispVec[i] = c1DispVec[nGrids-(int)(rc*200)-1];
			}
			*/
			FDispSpherical(FDispVec, c1DispVec,rVec,denOutVec,R,dr,sig,T,rc); 
		} 


		/*
		for (i=nGrids-(int)(rc*200); i<nGrids; i++) {
			c1HSVec[i] = c1HSVec[nGrids-(int)(rc*200)-1];
			c1DispVec[i] = c1DispVec[nGrids-(int)(rc*200)-1];
			c1AssocVec[i] = c1AssocVec[nGrids-(int)(rc*200)-1];
			FHSVec[i] = FHSVec[nGrids-(int)(rc*200)-1];
			FDispVec[i] = FDispVec[nGrids-(int)(rc*200)-1];
			FAssocVec[i] = FAssocVec[nGrids-(int)(rc*200)-1];
		}
		*/

		/* func */
		for (i=0; i<nGrids; i++) {
			funcVec[i] = 4*PI*rVec[i]*rVec[i]*exp(-c1HSVec[i]-c1DispVec[i]-c1AssocVec[i]);
		}
		func = trapz(funcVec,0,R,nGrids);

		/* Ntot */
		for (i=0; i<nGrids; i++) {
			funcVec[i] = 4*PI*rVec[i]*rVec[i]*denOutVec[i];
		}
		Ntot = trapz(funcVec,0,R,nGrids);

		dv = denOutVec[nGrids-(int)(rc*200)-1];
		mu = fmu_id(dv) + fmu_hs(dv) + fmu_disp(dv,get_a(T,rc)) + fmu_assoc(dv,Ma,K,T_assoc);

		//mu = log(Ntot/func);
		printf("Ntot = %f\n",Ntot);
		printf("mu = %f\n",mu);


		for (i=0; i<nGrids; i++) {
			FVec[i] = 4*PI*rVec[i]*rVec[i]*(FIdealVec[i] + FHSVec[i] + FDispVec[i] + FAssocVec[i] + denOutVec[i]*(-mu));
		}
		F = trapz(FVec,0,R,nGrids);


		// mu and p
		// find dl, dv from mu
		double den_before, den_after;
		double mu_before, mu_after;
		double den1, den2, den3;
		double mu1, mu2, mu3;
		double tol, buf;
		tol = 1.e-8;
		buf = 1.;

		for (i=0; i<(int)(1.2/0.000001); i++) {
			den_before = 1.2 - (i+1)*0.000001;
			den_after = 1.2 - (i)*0.000001;

			mu_before = fmu_id(den_before) + fmu_hs(den_before) + fmu_disp(den_before,get_a(T,rc)) + fmu_assoc(den_before,Ma,K,T_assoc);
			mu_after = fmu_id(den_after) + fmu_hs(den_after) + fmu_disp(den_after,get_a(T,rc)) + fmu_assoc(den_after,Ma,K,T_assoc);

			if ((mu-mu_before)*(mu-mu_after)<0) {
				den1 = den_before;
				den2 = den_after;
				den3 = (den_before+den_after)*0.5;

				while (buf > tol) {
					mu1 = fmu_id(den1) + fmu_hs(den1) + fmu_disp(den1,get_a(T,rc)) + fmu_assoc(den1,Ma,K,T_assoc);
					mu2 = fmu_id(den2) + fmu_hs(den2) + fmu_disp(den2,get_a(T,rc)) + fmu_assoc(den2,Ma,K,T_assoc);
					mu3 = fmu_id(den3) + fmu_hs(den3) + fmu_disp(den3,get_a(T,rc)) + fmu_assoc(den3,Ma,K,T_assoc);

					if ((mu-mu1)*(mu-mu3) < 0) {
						den2 = den3;
						den3 = (den1+den2)*0.5;
					} else if ((mu-mu3)*(mu-mu2) < 0) {
						den1 = den3;
						den3 = (den1+den2)*0.5;
					}
					buf = fabs(mu1-mu2);
				}

				dl = den1;
				break;
			}
		}

		printf("dl, dv = %f\t%f\n",dl,dv);
		printf("%f\n",fmu_id(dv)+fmu_hs(dv)+fmu_disp(dv,get_a(T,rc))+fmu_assoc(dv,Ma,K,T_assoc));
		printf("%f\n",fmu_id(dl)+fmu_hs(dl)+fmu_disp(dl,get_a(T,rc))+fmu_assoc(dl,Ma,K,T_assoc));

		pl = fp_id(dl) + fp_hs(dl) + fp_disp(dl,get_a(T,rc)) + fp_assoc(dl,Ma,K,T_assoc);
		pv = fp_id(dv) + fp_hs(dv) + fp_disp(dv,get_a(T,rc)) + fp_assoc(dv,Ma,K,T_assoc);




		/*


		FIdealSpherical(FIdealVec, denOutVec, R, d, dr);

		c1HSSpherical(c1HSVec,rVec,denOutVec,R,rc); 
		FHSSpherical(FHSVec, rVec, denOutVec, R, d,rc, input.FHSFlag); 

		if (Ma != 0) {
			c1AssocSpherical(c1AssocVec,rVec,denOutVec,R,Ma,K,T_assoc, d,dr, rc, input.FHSFlag, input.wallFlag);
			FAssocSpherical(FAssocVec, rVec, denOutVec, R,Ma, K, T_assoc, d, dr, rc, input.FHSFlag); 
		} 

		if (input.FDispFlag!=0) {
			c1DispSpherical(c1DispVec, rVec, denOutVec,R,T, rc);
			for (i=nGrids-(int)(rc*200); i<nGrids; i++) {
				c1DispVec[i] = c1DispVec[nGrids-(int)(rc*200)-1];
			}
			FDispSpherical(FDispVec, c1DispVec,rVec,denOutVec,R,dr,sig,T,rc); 
		} 


		for (i=0; i<nGrids; i++) {
			FVec[i] = 0;
			FVec[i] = 4*PI*rVec[i]*rVec[i]*(FIdealVec[i] + FHSVec[i] + FDispVec[i] + FAssocVec[i] + denOutVec[i]*(-mu));
		}
		F = trapz(FVec,0,R,nGrids);
		*/





		
		
		//printf("%f\t%f\n",FIdealVec[nGrids-2],FHSVec[nGrids-2]);
		printf("pl, pv = %f\t%f\n",pl,pv);

		V = (4.*PI/3.)*pow(R,3.);
		printf("R = %f\n",R);


		Rs = 3*(F+pv*V)/(2*PI*(pl-pv));
		Rs = pow(Rs,1./3.);

		for (i=0; i<nGrids; i++) {
			funcVec[i] = rVec[i]*rVec[i]*(denOutVec[i]-dv);
		}
		Re = (3/(dl-dv))*trapz(funcVec,0,R,nGrids);
		Re = pow(Re,1./3.);

		A = (4.*PI)*Re*Re;
		tension_Re = (F+pv*V)/A + ((pl-pv)*Re)/3.;

		tension = 3*(F+pv*V)*pow(pl-pv,2.)/(16.*PI);
		tension = pow(tension,1./3.);




		printf("%f\t%f\n",F,pv*V);


	} else {
/*
		for (i=0; i<M; i++) {
			for (j=0; j<nGrids; j++) {
				if (i==0) {
					denVec[j] = denOutMat[i][j];
				} else {
					denVec[j] += denOutMat[i][j];
				}
			}
		}

		double *fVec;
		fVec = (double*)malloc(nGrids*sizeof(double));
		for (i=0; i<nGrids; i++) {
			r = rVec[i];
			den = denVec[i];
			fVec[i] = 3/(dl-dv)*r*r*(den-dv);
		}

		Re = trapz(fVec,0,rVec[nGrids-1],nGrids);
		Re = pow(Re,1./3.); // Now we know Re
		free(fVec);

		//F_id = FIdealFlat(denVec, R, d, dr);

		F_hs = FHSFlat(rVec, denVec, R, d, rc, input.FHSFlag); 

		if (Ma != 0) {
			F_assoc = FAssocFlat(rVec,denVec,R,Ma,K,T_assoc,d,dr, rc,input.FHSFlag);
		} else {
			F_assoc = 0;
		}

		F_chain = FChainFlat(rVec,denVec,R,M,d,dr,rc,input.FHSFlag); 

		c1DispFlat(c1DispVec, rVec, denVec, R, T, rc); 
		F_disp = FDispFlat(c1DispVec, rVec, denVec, R, dr, sig, T, rc);

		for (i=0; i<nGrids; i++) {
			funcVec[i] = 0;
		}
		for (i=0; i<nGrids; i++) {
			funcVec[i] = - denMat[0][i] - denVec[i]*(c1HSMat[0][i]+c1ChainMat[0][i]+c1DispVec[i]);
		}
		F_mu = trapz(funcVec,0,R,nGrids);



		F = (F_hs + F_assoc + F_chain + F_disp) + F_mu;

		printf("%f\t%f\t%f\t%f\t%f\t%f\n",F_hs,F_assoc,F_chain,F_disp,F_mu,F);

		tension = F + pv*R;
		*/

	}




	char *opFileName = strcat(oFileName,".energy");
	f = fopen(opFileName,"w");
	fprintf(f, "%lf\t %lf\t %lf\t %lf \t %lf\t %lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",F,Re,Rs,Re-Rs,mu,dl,dv,pl,pv,tension_Re, tension);
	fclose(f);

	free_dVector(rVec);
	free_dVector(denVec);
	free_dVector(denInVec);
	free_dVector(denOutVec);
	free_dVector(c1HSVec);
	free_dVector(c1AssocVec);
	free_dVector(c1DispVec);
	free_dVector(tolVec);
	free_dVector(wallVec);
	free_dVector(funcVec);

	free_dVector(FIdealVec);
	free_dVector(FHSVec);
	free_dVector(FDispVec);
	free_dVector(FAssocVec);
	free_dVector(FVec);


	if (M > 1) {
		free_dMatrix(denMat,M);
		free_dMatrix(denInMat,M);
		free_dMatrix(denOutMat,M);
		free_dMatrix(c1HSMat,M);
		free_dMatrix(c1DispMat,M);
		free_dMatrix(c1ChainMat,M);
		free_dMatrix(GMat,M);
	}

	printf("Total time is %f s\n", ((double)clock()-before)/CLOCKS_PER_SEC);


}


