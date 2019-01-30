#include <stdlib.h>

int * iVector(int n) {
	int i;
	int *r;

	r = (int *)malloc(n*sizeof(int));
	for (i=0; i<n; i++) {
		r[i] = 0.;
	}
	return r;
}

double * dVector(int n) {
	int i;
	double *r;

	r = (double *)malloc(n*sizeof(double));
	for (i=0; i<n; i++) { 
		r[i] = 0.;
	}

	return r;
}

double ** dMatrix(int row, int col) {
	int i, j;
	double **r;

	r = (double **)malloc(row*sizeof(double*));
	for (i=0; i<row; i++) {
		r[i] = (double *)malloc(col*sizeof(double));
	}

	for (i=0; i<row; i++) {
		for (j=0; j<col; j++) {
			r[i][j] = 0;
		}
	}

	return r;
}

void free_iVector(int *r) {
	free(r);
}

void free_dVector(double *r) {
	free(r);
}

void free_dMatrix(double **r, int row) {
	int i, j;

	for (i=0; i<row; i++) {
		free(r[i]);
	}
	free(r);
}

