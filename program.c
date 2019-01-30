#include <stdio.h>
#include <math.h>
#include "math2.h"
#include "memory.h"
#include "constant.h"



double get_a(double T,double rc) {
	double rMin, a;
	rMin = pow(2,1./6.);
	a = 2*(PI*pow(rMin,3))/(3*T);
	a += (8.)*PI*(1/T)*(1./(9*pow(rc,9.))-1./(3*pow(rc,3))- 1./(9*pow(rMin,9.)) + 1./(3*pow(rMin,3)));
	return a;
	}

double fmu_id(double rho) {
	return log(rho);
}

double fp_id(double rho) {
	return rho;
}

double fmu_hs(double rho) {
	double eta;
	eta = rho*PI/6.;
	return eta*(8.-9.*eta+3.*eta*eta)/(pow(1.-eta,3.));
}

double fp_hs(double rho) {
	double eta;
	eta = rho*PI/6.;
	return rho*(4.*eta-2.*eta*eta)/(pow(1.-eta,3.));
}

double fmu_disp(double rho, double a) {
	return -2.*a*rho;
}

double fp_disp(double rho, double a) {
	return -a*rho*rho;
}

///////
double get_ghs_bulk(double rho) {
	return 1./(1.-PI*rho/6.) + (1./4.)*((PI*rho)/(pow(1.-PI*rho/6.,2.))) + (1./72.)*(pow(PI*rho,2.))/(pow(1.-PI*rho/6.,3.));
}

double get_gphs_bulk(double rho) {
	return (5.*PI/12.)/(pow(1.-PI*rho/6.,2.)) + (PI*PI/9.)*(rho)/(pow(1.-PI*rho/6.,3.)) + (PI*PI*PI/144.)*(pow(rho,2.))/(pow(1.-PI*rho/6.,4.));
}

double get_D_bulk(double rho,double K,double T_assoc) {
	double ghs;
	ghs = get_ghs_bulk(rho);
	return 4.*PI*K*ghs*(exp(1./T_assoc)-1.);
}

double get_Dp_bulk(double rho,double K,double T_assoc) {
	double gphs;
	gphs = get_gphs_bulk(rho);
	return 4.*PI*K*gphs*(exp(1./T_assoc)-1.);
}

double get_X_bulk(double rho,double Ma,double K,double T_assoc) {
	double D;
	D = get_D_bulk(rho,K,T_assoc);
	//return (-1.+sqrt(1.+2.*Ma*rho*D))/(Ma*rho*D);
	return (-1.+sqrt(1.+4.*Ma*rho*D))/(2.*Ma*rho*D);
}

double get_Xp_bulk(double rho,double Ma,double K,double T_assoc) {
	double D;
	double Dp;
	D = get_D_bulk(rho,K,T_assoc);
	Dp = get_Dp_bulk(rho,K,T_assoc);
	//return (1./(rho*D*sqrt((1.+2.*Ma*rho*D))))*(D+rho*Dp) + ((1.-sqrt(1.+2.*Ma*rho*D))/(Ma*rho*rho*D*D))*(D+rho*Dp);
	return (1./(rho*D*sqrt(1.+4.*Ma*rho*D)))*(D+rho*Dp) + ((1.-sqrt(1.+4.*Ma*rho*D))/(2.*Ma*rho*rho*D*D))*(D+rho*Dp);
}


double fmu_assoc(double rho,double Ma,double K,double T_assoc) {
	double X, Xp;
	if (Ma == 0) {
		return 0;
	} else {
		X = get_X_bulk(rho,Ma,K,T_assoc);
		Xp = get_Xp_bulk(rho,Ma,K,T_assoc);

		return Ma*(log(X)-0.5*X+0.5) + Ma*rho*Xp*((1./X) - 0.5);
	}
}

double fp_assoc(double rho,double Ma,double K,double T_assoc) {
	double X, Xp;
	if (Ma == 0) {
		return 0;
	} else {
		X = get_X_bulk(rho,Ma,K,T_assoc);
		Xp = get_Xp_bulk(rho,Ma,K,T_assoc);

		return  -Ma*rho*(log(X) - 0.5*X + 0.5) + rho*(Ma*(log(X)-0.5*X+0.5) + Ma*rho*Xp*((1/X) - 0.5));
	}
}



/*
void find_peaks_mu(double *muVec, double T,double a) {
	int i, count;
	double f, f_before, f_after;
	double mu_max, mu0, mu_min;
	double *denVec, *test_solVec;
	denVec = dVector(10000);
	for (i=0; i<10000; i++) {
		denVec[i] = 0.0001*i;
	}

	test_solVec = dVector(2);
	count = 0;
	for (i=1; i<10000-1; i++) {
		f = fmu_id(denVec[i]) + fmu_hs(denVec[i]) + fmu_disp(denVec[i],a);
		f_before = fmu_id(denVec[i-1]) + fmu_hs(denVec[i-1]) + fmu_disp(denVec[i-1],a);
		f_after = fmu_id(denVec[i+1]) + fmu_hs(denVec[i+1]) + fmu_disp(denVec[i+1],a);

		if ((f-f_before)*(f_after-f)<0) {
			test_solVec[count] += denVec[i];
			count++;
			if (count > 1) {
				printf("Too much dl and dv\n");
				exit(1);
			}
		}
	}

	if (test_solVec[1] == 0) {
		printf("Too small dl and dv\n");
		exit(1);
	}
	

	mu_max = fmu_id(test_solVec[0]) + fmu_hs(test_solVec[0]) + fmu_disp(test_solVec[0],a);
	mu0 = fmu_id(denVec[0]) + fmu_hs(denVec[0]) + fmu_disp(denVec[0],a);
	mu_min = fmu_id(test_solVec[1]) + fmu_hs(test_solVec[1]) + fmu_disp(test_solVec[1],a);
	if (mu_min<mu0) {
		mu_min = mu0;
	}

	free_dVector(denVec);
	free_dVector(test_solVec);

	muVec[0] = mu_max;
	muVec[1] = mu_min;
}
*/


void find_roots_mu(double *dVec, double T, double a, double mu_test) {
	int i, step, count;
	double tol, tol_ref;
	double sol1, sol2, sol3;
	double f, f_before;
	double f1, f2, f3;
	double *denVec, *test_solVec;
	double rho_v, rho_l;

	tol_ref = 1e-6;
	denVec = dVector(10000);
	for (i=0; i<10000; i++) {
		denVec[i] = i*0.0001;
	}

	test_solVec = dVector(4);

	count = 0;
	for (i=1; i<10000; i++) {
		f = fmu_id(denVec[i]) + fmu_hs(denVec[i]) + fmu_disp(denVec[i],a) - mu_test;
		f_before = fmu_id(denVec[i-1]) + fmu_hs(denVec[i-1]) + fmu_disp(denVec[i-1],a) - mu_test;
		if (f>0 && f_before<0) {
			test_solVec[count] = denVec[i-1];
			test_solVec[count+1] = denVec[i];
			if (count > 3) {
				printf("Too much solution\n");
				exit(1);
			}
			count += 2;
		}
	}
	
	tol = 1.;
	step = 0;
	while (tol>tol_ref) {
		if (step == 0) {
			sol1 = test_solVec[0];
			sol2 = test_solVec[1];
		}

		sol3 = (sol1+sol2)/2;

		f1 = fmu_id(sol1) + fmu_hs(sol1) + fmu_disp(sol1,a) - mu_test;
		f2 = fmu_id(sol2) + fmu_hs(sol2) + fmu_disp(sol2,a) - mu_test;
		f3 = fmu_id(sol3) + fmu_hs(sol3) + fmu_disp(sol3,a) - mu_test;


		if (f3>0 && f1<0) {
			sol2 = sol3;
		} else if (f2>0 && f3<0) {
			sol1 = sol3;
		} else {
			printf("Something is wrong..\n");
			exit(1);
		}
		tol = fabs(sol2-sol1);
		step++;
	}

	rho_v = (sol1+sol2)*0.5;

	tol = 1.;
	step = 0;
	while (tol>tol_ref) {
		if (step == 0) {
			sol1 = test_solVec[2];
			sol2 = test_solVec[3];
		}

		sol3 = (sol1+sol2)/2;

		f1 = fmu_id(sol1) + fmu_hs(sol1) + fmu_disp(sol1,a) - mu_test;
		f2 = fmu_id(sol2) + fmu_hs(sol2) + fmu_disp(sol2,a) - mu_test;
		f3 = fmu_id(sol3) + fmu_hs(sol3) + fmu_disp(sol3,a) - mu_test;



		if (f3>0 && f1<0) {
			sol2 = sol3;
		} else if (f2>0 && f3<0) {
			sol1 = sol3;
		} else {
			printf("Something is wrong..\n");
			exit(1);
		}
		tol = abs(sol2-sol1);
		step++;
	rho_l = (sol1+sol2)*0.5;
	}

	dVec[0] = rho_l;
	dVec[1] = rho_v;

	free_dVector(denVec);
	free_dVector(test_solVec);

}





int checkTol(double *tolVec, double *denVec, double *denOutVec, int nGrids, double tol) {
	int i, tolFlag;
	double r;
	double max;
	for (i=0; i<nGrids; i++) { tolVec[i] = 0; }
	for (i=0; i<nGrids; i++) {
		tolVec[i] = fabs(denOutVec[i]-denVec[i]);
	}
	max = tolVec[0];

	for (i=1; i<nGrids; i++) {
		if (tolVec[i]>max) {
			max = tolVec[i];
		}
	}

	if (max>tol) {
		tolFlag = 0;
	} else {
		tolFlag = 1;
	}

	printf("Max is %f.\n", max);


	return tolFlag;

}

double returnTol(double *tolVec, double *denVec, double *denOutVec, int nGrids, double tol) {
	int i, tolFlag;
	double r;
	double max;
	for (i=0; i<nGrids; i++) { tolVec[i] = 0; }
	for (i=0; i<nGrids; i++) {
		tolVec[i] = fabs(denOutVec[i]-denVec[i]);
	}
	max = tolVec[0];

	for (i=1; i<nGrids; i++) {
		if (tolVec[i]>max) {
			max = tolVec[i];
		}
	}

	if (max>tol) {
		tolFlag = 0;
	} else {
		tolFlag = 1;
	}

	printf("Max is %f.\n", max);

	return max;

}

void updateDen(double *denVec, double *denInVec, double *denOutVec, int nGrids, double q) {
	int i;
	for (i=0; i<nGrids; i++) {
		denInVec[i] = denVec[i]*(1.0-q) + q*denOutVec[i];
	}
}
