#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Prep.h"

// Function to read the parameters from the input text file
void ReadParam(char *ParIn) {
	int nCa, nCon;
	FILE *ParFile;
	ParFile = fopen(ParIn,"r");
	if (ParFile == NULL) {
        printf("Cannot open parameter file.\n");
        exit(0);
    }
    else {
    	thresP = 1.0; //defaul taking all SNPs, if no threshold specified
    	while (fgets(buffer, sizeof(buffer), ParFile) != NULL) {
	    	char *tok = strtok(buffer," \t\n");
	    	char *toktmp;
	    	if (tok != NULL) {
	    		if (strcmp(tok, "Target") == 0) {
	    			strcpy(bufftmp, strtok(NULL, " \t\n"));
	    			toktmp = strtok(bufftmp,", \t\n");
	    			nInfile = 0;
	    			while( toktmp != NULL ) {
				    	strcpy(TargetIn[nInfile], toktmp);
				  		toktmp = strtok(NULL, ", \t\n");
				  		nInfile++;
				  	}
				}
	    		else if (strcmp(tok, "Base") == 0)
	    			strcpy(Base, strtok(NULL, " \t\n"));
	    		else if (strcmp(tok, "Output") == 0)
	    			strcpy(Output, strtok(NULL, " \t\n"));
	    		else if (strcmp(tok, "Pthres") == 0)
	    			thresP = atof(strtok(NULL, " \t\n"));
	    		else if (strcmp(tok, "nCase") == 0) {
	    			strcpy(bufftmp, strtok(NULL, " \t\n"));
	    			toktmp = strtok(bufftmp,", \t\n");
	    			nCa = 0;
	    			while( toktmp != NULL ) {
	    				nCase[nCa++] = atoi(toktmp);
				  		toktmp = strtok(NULL, ", \t\n");
				  	}
	    		}
	    		else if (strcmp(tok, "nControl") == 0) {
	    			strcpy(bufftmp, strtok(NULL, " \t\n"));
	    			toktmp = strtok(bufftmp,", \t\n");
	    			nCon = 0;
	    			while( toktmp != NULL ) {
	    				nControl[nCon++] = atoi(toktmp);
				  		toktmp = strtok(NULL, ", \t\n");
				  	}
	    		}
	    	}
    	}

    	if ((nCa != nInfile) || (nCon != nInfile)) {
    		printf("Missing parameter!\n");
    		exit(0);
    	}
    }
}


// Function to read the Base data into Hash table
void ReadBase(char *InputBase, FILE *LogFile) {
	long int i = 0;
	Frq Basefrq;
	FILE *Base;
	char snp[50], a1[100], a2[100];
	Base = fopen(InputBase, "r");
	if (Base == NULL) {
	    printf("Cannot open PRS input file.\n");
	    exit(0);
	}
	else {
		int SNPc = 0, CHRc = 0, Posc = 0, ORc = 0, SEc = 0, Affc = 0, Unaffc = 0, BETAc = 0, Pc = 0;
		char *tok;
		fgets(buffer, sizeof(buffer), Base);
		tok = strtok(buffer," \t\n");
		while (tok != NULL) {
			if (strcmp(tok, "SNP") == 0)
				SNPc = i+1;
			else if (strcmp(tok, "CHR") == 0)
				CHRc = i+1;
			else if (strcmp(tok, "BP") == 0)
				Posc = i+1;
			else if (strcmp(tok, "OR") == 0)
				ORc = i+1;
			else if (strcmp(tok, "SE") == 0)
				SEc = i+1;
			else if (strcmp(tok, "P") == 0)
				Pc = i+1;
			else if (strcmp(tok, "A1") == 0)
                		Affc = i+1;
           		 else if (strcmp(tok, "A2") == 0)
                		Unaffc = i+1;
           		 else if (strcmp(tok, "Beta") == 0)
                		BETAc = i+1;
			tok = strtok(NULL, " \t\n");
			i++;
		}
		if (!(SNPc && CHRc && Posc && ( ORc || BETAc) && SEc && Pc)) {
			printf("Missing token\n");
			exit(0);
		}

		i = 0;
		while (fgets(buffer, sizeof(buffer), Base) != NULL) {
			int k = 0;
			char *token[20];
			char *p = strtok(buffer," \t\n");
			while (p != NULL) {
				token[k++] = p;
				p = strtok(NULL," \t\n");
			}
			Data tmp;
		    	strcpy(snp,token[SNPc-1]);
		    	strcpy(a1,token[Affc-1]);
		    	strcpy(a2,token[Unaffc-1]);
		    	convertToUpperCase(a1);
		    	convertToUpperCase(a2);
		    	strcpy(tmp.SNP,snp);
		    	strcpy(tmp.Aff,a1);
		    	strcpy(tmp.Unaff,a2);
		    	tmp.CHR = atoi(token[CHRc-1]);
		    	tmp.Pos = atoi(token[Posc-1]);
		    	if (BETAc)
				tmp.Beta = atof(token[BETAc-1]);
		    	else
				tmp.Beta = log(atof(token[ORc-1]));
			    // tmp.Beta = log(atof(token[ORc-1]))/atof(token[SEc-1]); //possible shrinkage?
			if (atof(token[Pc-1]) <= thresP) {
				hashSNPpush(tmp);
				i++;
			}
		}
	}
	fclose(Base);
	fprintf(LogFile, "Analysis Starts.\n");
	fprintf(LogFile, "P value threshold for base SNPs : %.2e.\n", thresP);
	fprintf(LogFile, "%ld SNPs below P threshold read from base.\n", i);
}




