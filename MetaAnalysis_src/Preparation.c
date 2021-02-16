#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Preparation.h"

// Hash function to carry out over the base pair position of a SNP
long int hashFunc(long int i) {
	return i%10000;
}

// Function to push a SNP into the hash table
void hashSNPpush(Data SNP) {
	long int ind = hashFunc(SNP.Pos);
	hashTable[SNP.CHR-1][ind][hashLen[SNP.CHR-1][ind]] = SNP;
	hashLen[SNP.CHR-1][ind]++;
}

// Function to search for a SNP by position and rsID from the hash table and return the index off its position or 0
Data *hashSNPsearch(long int CHR, long int Pos, char *SNPID) {
	Data *tmp;
	tmp = malloc(sizeof(Data));
	long int ind = hashFunc(Pos);
	long int i;
	for (i = 0; i < hashLen[CHR-1][ind]; i++) {
		tmp = &hashTable[CHR-1][ind][i];
		if (strcmp(SNPID, (*tmp).SNP) == 0)
			return tmp;
	}
	return NULL;
}

// Function to read the parameters from the input text file
void ReadParam(char *ParIn) {
	FILE *ParFile;
	ParFile = fopen(ParIn,"r");
	int nInfile, nInCaCa, nInCaCon, nInConCon;
	if (ParFile == NULL) {
        printf("Cannot open parameter file.\n");
        exit(0);
    }
    else {
    	while (fgets(buffer, sizeof(buffer), ParFile) != NULL) {
	    	char *tok = strtok(buffer," \t\n");
	    	char *toktmp;
	    	size = 2; //default analyzing 2 datasets, using 1 indicator vector
	    	if (tok != NULL) {
	    		if (strcmp(tok, "Input") == 0) {
	    			strcpy(bufftmp, strtok(NULL, " \t\n"));
	    			toktmp = strtok(bufftmp,", \t\n");
	    			nInfile = 0;
	    			while( toktmp != NULL ) {
				    	strcpy(Input[nInfile], toktmp);
				  		toktmp = strtok(NULL, ", \t\n");
				  		nInfile += 1;
				  	}
				  }
	    		else if (strcmp(tok, "Output") == 0)
	    			strcpy(Output, strtok(NULL, " \t\n"));
	    		else if (strcmp(tok, "Firth") == 0)
	    			FirthThres = atof(strtok(NULL, " \t\n"));
	    		else if (strcmp(tok, "Zthres") == 0)
	    			Zthres = atof(strtok(NULL, " \t\n"));
	    		else if (strcmp(tok, "CaseInCase") == 0) {
	    			strcpy(bufftmp, strtok(NULL, " \t\n"));
	    			toktmp = strtok(bufftmp,", \t\n");
	    			nInCaCa = 0;
	    			while( toktmp != NULL ) {
	    				nCaCa[nInCaCa] = atoi(toktmp);
				  		toktmp = strtok(NULL, ", \t\n");
				  		nInCaCa += 1;
				  	}
	    		}
	    		else if (strcmp(tok, "ControlInControl") == 0) {
	    			strcpy(bufftmp, strtok(NULL, " \t\n"));
	    			toktmp = strtok(bufftmp,", \t\n");
	    			nInConCon = 0;
	    			while( toktmp != NULL ) {
	    				nConCon[nInConCon] = atoi(toktmp);
				  		toktmp = strtok(NULL, ", \t\n");
				  		nInConCon += 1;
				  	}
	    		}
	    		else if (strcmp(tok, "CaseInControl") == 0) {
	    			strcpy(bufftmp, strtok(NULL, " \t\n"));
	    			toktmp = strtok(bufftmp,", \t\n");
	    			nInCaCon = 0;
	    			while( toktmp != NULL ) {
	    				nCaCon[nInCaCon] = atoi(toktmp);
				  		toktmp = strtok(NULL, ", \t\n");
				  		nInCaCon += 1;
				  	}
	    		}
	    		/*else if (strcmp(tok, "nFiles") == 0)
	    			size = atoi(strtok(NULL, " \t\n"));*/
	    	}
    	}
    	if ((nInCaCa != pow(nInfile,2)) || (nInConCon != pow(nInfile,2)) || (nInCaCon != pow(nInfile,2))) {
    		printf("Missing parameter!\n");
    		exit(0);
    	}
    }
    size = nInfile;
    int i = 0; int j = 0;
    for (i = 0; i < size; i++) {
    	for (j = 0; j < size; j++) {
    		CaCa[i][j] = nCaCa[i*size+j];
    		ConCon[i][j] = nConCon[i*size+j];
    		CaCon[i][j] = nCaCon[i*size+j];
    	}
    }

}









