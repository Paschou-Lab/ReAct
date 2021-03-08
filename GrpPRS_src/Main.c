#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cdflib.h"
#include "Prep.h"
#include "CountConstruct.h"
#include "SupportFunc.h"


//arguments: 1-Experiment index, 2-Input Base file, 3-Input Map file directory, 4-Output, 5-seed, 6-maxGen
int main(int argc, char* argv[]) {

    int m; long int n;
	ReadParam(argv[1]);
    for (m = 0; m < 23; m++) {
        for (n = 0; n < 10000; n++) {
            hashTable[m][n] = malloc(500 * sizeof(struct data));
            hashLen[m][n] = 0;
        }
    }

	FILE *OutFile;
	FILE *LogFile;
	OutFile = fopen(Output,"w");
	LogFile = fopen(strcat(Output,".log"),"w");
	ReadBase(Base, LogFile);

	if (OutFile == NULL) {
        printf("Cannot open output file.\n");
        exit(0);
    } //check first wether the output file can be opened
	fprintf(OutFile, "InFile\tPthres\tnSNPs\tCasePRS\tControlPRS\tStudyPRS\tCasePRS_SE\tControlPRS_SE\tStudyPRS_SE\tPval\n");

	char snp[50], a1[50], a2[50];
	int chr, k, j;
	long int pos, ncase, ncontrol;;
	double or, se, freq;
	Frq DFfrq;
    double VarCaS, VarConS, VarPopS, Tstat, df;
	Data *tmpSNP;
	tmpSNP = malloc(sizeof(Data));
	FILE *InFile;

	for (k = 0; k < nInfile; k++) {
        InFile = fopen(TargetIn[k], "r");
        if (InFile == NULL) {
            printf("Cannot open the %d-th input file.\n", k+1);
            exit(0);
        }
        else {
        	ScoreCase[k] = 0;
        	ScoreControl[k] = 0;
        	ScoreDF[k] = 0;
        	nSNPDf[k] = 0;
        	Pval[k] = 1;
            VarCaS = VarConS = VarPopS = 0;

            int SNPc = 0, Affc = 0, Unaffc = 0, CHRc = 0, Posc = 0, ORc = 0, BETAc = 0, SEc = 0, nCasec = 0, nControlc = 0, Frqc = 0;
            char *tok;
            fgets(buffer, sizeof(buffer), InFile);
            int i = 0;
            tok = strtok(buffer," \t\n");
            while (tok != NULL) {
                if (strcmp(tok, "SNP") == 0)
                    SNPc = i+1;
                else if (strcmp(tok, "CHR") == 0)
                    CHRc = i+1;
                else if (strcmp(tok, "BP") == 0)
                    Posc = i+1;
                else if (strcmp(tok, "A1") == 0)
                    Affc = i+1;
                else if (strcmp(tok, "A2") == 0)
                    Unaffc = i+1;
                else if (strcmp(tok, "OR") == 0)
                    ORc = i+1;
                else if (strcmp(tok, "Beta") == 0)
                    BETAc = i+1;
                else if (strcmp(tok, "SE") == 0)
                    SEc = i+1;
                else if (strcmp(tok, "nCase") == 0)
                    nCasec = i+1;
                else if (strcmp(tok, "nControl") == 0)
                    nControlc = i+1;
                else if (strcmp(tok, "Frq") == 0)
                    Frqc = i+1;
                tok = strtok(NULL, " \t\n");
                i++;
            } // read header 

            if (!(SNPc && CHRc && Posc && Affc && Unaffc && ORc && SEc)) {
                printf("Missing token\n");
                exit(0);
            } // input format sanity check

            while (fgets(buffer, sizeof(buffer), InFile) != NULL) {
                char *token[i];
                char *p = strtok(buffer," \t\n");
                j = 0;
                while (p != NULL) {
                    token[j++] = p;
                    p = strtok(NULL," \t\n");;
                }
                strcpy(snp,token[SNPc-1]);
                strcpy(a1,token[Affc-1]);
                strcpy(a2,token[Unaffc-1]);
                convertToUpperCase(a1);
                convertToUpperCase(a2);
                chr = atoi(token[CHRc-1]);
                pos = atoi(token[Posc-1]);
                if (ORc)
                    or = atof(token[ORc-1]);
                else if (BETAc)
                    or = exp(atof(token[BETAc-1]));
                se = atof(token[SEc-1]);

                if (!Frqc)
                    freq = 0.0f;
                else
                    freq = atof(token[Frqc-1]);
                if (!nCasec)
                    ncase = nCase[k];
                else
                    ncase = atoi(token[nCasec-1]);
                if (!nControlc)
                    ncontrol = nControl[k];
                else
                    ncontrol = atoi(token[nControlc-1]);

                tmpSNP = hashSNPsearch(chr, pos, snp);
                if (tmpSNP) {
                	DFfrq = GroupFreq(se, ncase, ncontrol, or, freq);
                	if ((strcmp(a1, (*tmpSNP).Aff) == 0) && (strcmp(a2, (*tmpSNP).Unaff) == 0)) {
                        ScoreCase[k] += (*tmpSNP).Beta * DFfrq.pCa;
                        VarCaS += pow((*tmpSNP).Beta,2)*DFfrq.pCa*(1-DFfrq.pCa);
                		ScoreControl[k] += (*tmpSNP).Beta * DFfrq.pCon;
                        VarConS += pow((*tmpSNP).Beta,2)*DFfrq.pCon*(1-DFfrq.pCon);
                		ScoreDF[k] += (*tmpSNP).Beta * DFfrq.pPop;
                        VarPopS += pow((*tmpSNP).Beta,2)*DFfrq.pPop*(1-DFfrq.pPop);
                		nSNPDf[k]++;
                	}
                	else if ((strcmp(a2, (*tmpSNP).Aff) == 0) && (strcmp(a1, (*tmpSNP).Unaff) == 0)) {
                		ScoreCase[k] += (*tmpSNP).Beta * (1-DFfrq.pCa);
                        VarCaS += pow((*tmpSNP).Beta,2)*DFfrq.pCa*(1-DFfrq.pCa);
                		ScoreControl[k] += (*tmpSNP).Beta * (1-DFfrq.pCon);
                        VarConS += pow((*tmpSNP).Beta,2)*DFfrq.pCon*(1-DFfrq.pCon);
                		ScoreDF[k] += (*tmpSNP).Beta * (1-DFfrq.pPop);
                        VarPopS += pow((*tmpSNP).Beta,2)*DFfrq.pPop*(1-DFfrq.pPop);
                		nSNPDf[k]++;
                	}
                }
            }
            fclose(InFile);
            ScoreCase[k] = ScoreCase[k]/nSNPDf[k];
            SDCase[k] = sqrt(2*VarCaS/pow(2*nSNPDf[k],2));
            ScoreControl[k] = ScoreControl[k]/nSNPDf[k];
            SDControl[k] = sqrt(2*VarConS/pow(2*nSNPDf[k],2));
            ScoreDF[k] = ScoreDF[k]/nSNPDf[k];
            SDDF[k] = sqrt(2*VarPopS/pow(2*nSNPDf[k],2));
            df = nCase[k]+nControl[k]-2;
            Tstat = fabs(ScoreCase[k] - ScoreControl[k])/(SDDF[k]*sqrt(1/(double)nCase[k] + 1/(double)nControl[k]));
            // printf("df = %lf, Tstat = %lf\n", df, Tstat);
            Pval[k] = GetPval(Tstat, df);
            fprintf(OutFile, "%s\t%lf\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%.4e\n", TargetIn[k], thresP, nSNPDf[k], ScoreCase[k], ScoreControl[k], ScoreDF[k], SDCase[k], SDControl[k], SDDF[k], Pval[k]);
            fprintf(LogFile, "Study %s Finished, %ld SNPs taken for PRS computation.\n", TargetIn[k], nSNPDf[k]);
        }
    }
    fclose(OutFile);
    fclose(LogFile);
	return (0);
}
