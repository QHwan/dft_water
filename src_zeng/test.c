#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXSTRING 1024

int main() {
	int i;
	char line[MAXSTRING], word[MAXSTRING], null[MAXSTRING];
	char *w = NULL;
	int system_flag;
	FILE *file;

	file = fopen("test.txt","r");
	if (file==NULL) {
		printf("File pointer is NULL\n");
		exit(1);
	}

	while(fgets(line,MAXSTRING,file) != NULL) {
		strcpy(word,line);
		w = strtok(word,"\t");

		if (strcmp(word,"System") == 0) {
			sscanf(line,"%s %d",null,&system_flag);
		}
	}

	printf("System flag is %d.\n", system_flag);
	
}	
