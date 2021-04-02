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
    memset(nOverlapCa, 0, sizeof(nOverlapCa[0]) * maxDF);
    memset(nOverlapCon, 0, sizeof(nOverlapCon[0]) * maxDF);
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
	fprintf(OutFile, "InFile\tPthres\tnSNPs\tCasePRS\tControlPRS\tCasePRS_SE\tControlPRS_SE\tR2\tPval\n");

	char snp[50], a1[50], a2[50];
	int chr, k, j, flag;
	long int pos, ncase, ncontrol;;
	double or, se, freq, Zk;
	Frq DFfrq;
	Data *tmpSNP;
	tmpSNP = malloc(sizeof(Data));
	FILE *InFile;

////
    memset(SumZik, 0, sizeof(SumZik[0]) * maxDF);
    memset(SqrSumZk, 0, sizeof(SqrSumZk[0][0]) * maxDF * 2);
    memset(Nik, 0, sizeof(Nik[0]) * maxDF);
    memset(CorrR, 0, sizeof(CorrR[0]) * maxDF);
////

	for (k = 0; k < nInfile; k++) {
        InFile = fopen(TargetIn[k], "r");
        if (InFile == NULL) {
            printf("Cannot open the %d-th input file.\n", k+1);
            exit(0);
        }
        else {
        	ScoreCase[k] = ScoreControl[k] = ScoreDF[k] = 0;
            VarCaS[k] = VarConS[k] = VarPopS[k] = 0;
            BaseScoreCase[k] = BaseScoreControl[k] = BaseScore[k] = 0; 
            BaseVarCaS[k] = BaseVarConS[k] = BaseVarPopS[k] = 0;
        	nSNPDf[k] = 0;
        	Pval[k] = 1;

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

            if (!(SNPc && CHRc && Posc && Affc && Unaffc && ( ORc || BETAc) && SEc)) {
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
                    flag = -1;
                	DFfrq = GroupFreq(se, ncase, ncontrol, or, freq);
                	if ((strcmp(a1, (*tmpSNP).Aff) == 0) && (strcmp(a2, (*tmpSNP).Unaff) == 0)) {
                        flag = 1;
                        if (Zthres > 0.0)
                            UpdateZsum((*tmpSNP).Zbase, or, se, k);
                	}
                	else if ((strcmp(a2, (*tmpSNP).Aff) == 0) && (strcmp(a1, (*tmpSNP).Unaff) == 0)) {
                        flag = 0;
                        if (Zthres > 0.0)
                            UpdateZsum((*tmpSNP).Zbase, 1/or, se, k);
                	}
                    if (((*tmpSNP).Pbase <= thresP) && (flag != -1)) {
                        UpdateScore((*tmpSNP), DFfrq, k, flag);
                        nSNPDf[k]++;
                    }
                }
            }
            fclose(InFile);
        }
    }

    if (Zthres > 0.0) {
        GetCorrR(nInfile);
        for (k = 0; k < nInfile; k++) {
            nOverlapCa[k] = sqrt(nCase[k]*nCaBase)*CorrR[k];
            nOverlapCon[k] = sqrt(nControl[k]*nConBase)*CorrR[k];
        }
    }
    for (k = 0; k < nInfile; k++) {
        qCa[k] = nOverlapCa[k]/(double)nCase[k];
        qCon[k] = nOverlapCon[k]/(double)nControl[k];
        q[k] = (nOverlapCa[k]+nOverlapCon[k])/(double)(nCase[k] + nControl[k]);
        if (Zthres > 0.0) 
            fprintf(LogFile, "Est. sample overlap for %d-th target study and the base: %lf\n", k+1, q[k]);
    }

    for (k = 0; k < nInfile; k++) {
        UpdateObsStat(k);
        UpdateDisStat(k);
        ShiftBaseScore(k);
        if ( q[k] > 0.0) {
            UpdateRealStat(k);
            Pval[k] = GetPval(TstatReal[k], df[k]);
            fprintf(OutFile, "%s\t%lf\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%.4e\n", TargetIn[k], thresP, nSNPDf[k], ScoreCase[k], ScoreControl[k], SDCase[k], SDControl[k], R2Real[k], Pval[k]);
            fprintf(LogFile, "Study %s Finished, %ld SNPs taken for PRS computation.\n", TargetIn[k], nSNPDf[k]);
        }
        else {
            Pval[k] = GetPval(TstatObs[k], df[k]);
            fprintf(OutFile, "%s\t%lf\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%.4e\n", TargetIn[k], thresP, nSNPDf[k], ScoreCase[k], ScoreControl[k], SDCase[k], SDControl[k], R2Obs[k], Pval[k]);
            fprintf(LogFile, "Study %s Finished, %ld SNPs taken for PRS computation.\n", TargetIn[k], nSNPDf[k]);
        }
    }
    fclose(OutFile);
    fclose(LogFile);
	return (0);
}
