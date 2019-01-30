#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structure.h"
#include "constant.h"

void inputInitialize(INPUT *input) {
	input->size = 0;
	input->dr = 0;
	input->ensembleFlag = -1;
	input->N = 0;
	input->den = 0;
	input->den_v = 0;

	input->mu = 0;
	input->Ri = 0;

	input->step = 0;

	input->HSDiam = 0;
	input->sig = 0;
	input->T = 0;

	input->rCut = 0;
	
	input->FHSFlag = -1;
	input->FDispFlag = -1;

	input->q = 0;
	input->tol = 0;

	input->pressure = 0;
}

void readInput(char *fileName, INPUT *input) {
	int i;
	char line[MAXSTRING], word[MAXSTRING], null[MAXSTRING];
	char *w = NULL;
	FILE *file;

	file = fopen(fileName,"r");
	if (file==NULL) {
		printf("File pointer is NULL\n");
		exit(1);
	}

	while(fgets(line,MAXSTRING,file) != NULL) {
		strcpy(word,line);
		w = strtok(word,"\t");

		if (strcmp(word,"System") == 0) {
			sscanf(line,"%s %d",null,&input->systemFlag);
		}
		else if (strcmp(word,"Size") == 0) {
			sscanf(line,"%s %lf",null,&input->size);
		}
		else if (strcmp(word,"dr") == 0) {
			sscanf(line,"%s %lf",null,&input->dr);
		}
		else if (strcmp(word,"Ri") == 0) {
			sscanf(line,"%s %lf",null,&input->Ri);
		}
		else if (strcmp(word,"Step") == 0) {
			sscanf(line,"%s %lf",null,&input->step);
		}
		else if (strcmp(word,"mu") == 0) {
			sscanf(line,"%s %lf",null,&input->mu);
		}
		else if (strcmp(word,"Ensemble") == 0) {
			sscanf(line,"%s %s",null,w);
			if (strcmp(w,"NVT") == 0) {
				input->ensembleFlag = 1;
			}
			else if (strcmp(w,"UVT") == 0) {
				input->ensembleFlag = 2;
			}
			else {
				printf("You choose wrong ensemble!\n");
				exit(1);
			}
		}
		else if (strcmp(word,"N") == 0) {
			sscanf(line,"%s %d",null,&input->N);
		}
		else if (strcmp(word,"den") == 0) {
			if (input->ensembleFlag == 1) {
				sscanf(line,"%s %lf %lf",null,&input->den,&input->den_v);
			} else {
				sscanf(line,"%s %lf %lf",null,&input->den,&input->den_v);
			}
		}
		else if (strcmp(word,"HardSphere") == 0) {
			sscanf(line,"%s %s",null,w);
			if (strcmp(w,"No") == 0) {
				input->FHSFlag = 0;
			}
			else if (strcmp(w,"MFMT") == 0) {
				input->FHSFlag = 1;
			}
			else if (strcmp(w,"FMT") == 0) {
				input->FHSFlag = 2;
			}
			else {
				printf("You choose wrong HS energy functional!\n");
				exit(1);
			}
		}
		else if (strcmp(word,"HSDiam") == 0) {
			sscanf(line,"%s %lf",null,&input->HSDiam);
		}

		else if (strcmp(word,"Dispersion") == 0) {
			sscanf(line,"%s %s",null,w);
			if (strcmp(w,"No") == 0) {
				input->FDispFlag = 0;
			}
			else if (strcmp(w,"LJ") == 0) {
				input->FDispFlag = 1;
			}
			else {
				printf("You choose wrong Disp energy functional!\n");
				exit(1);
			}
		}
		else if (strcmp(word,"Sig") == 0) {
			sscanf(line,"%s %lf",null,&input->sig);
		}
		else if (strcmp(word,"Temperature") == 0) {
			sscanf(line,"%s %lf",null,&input->T);
		}
		else if (strcmp(word,"RCut") == 0) {
			sscanf(line,"%s %lf",null,&input->rCut);
		}
		else if (strcmp(word,"q") == 0) {
			sscanf(line,"%s %lf",null,&input->q);
		}	
		else if (strcmp(word,"Tolerance") == 0) {
			sscanf(line,"%s %lf",null,&input->tol);
		}
		else if (strcmp(word,"Pressure") == 0) {
			sscanf(line,"%s %d",null,&input->pressure);
		}
	}

	fclose(file);
}


