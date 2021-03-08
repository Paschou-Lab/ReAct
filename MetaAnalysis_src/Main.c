#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cdflib.h"
#include "Preparation.h"
#include "CountConstruct.h"
#include "PvalCompute.h"
#include "SupportFunc.h"


// Main function
int main(int argc,char* argv[]) {
	long int i=0; // SNP counter
	long int j=0; // others
	long int k=0; // input counter
	long int m=0, n=0; //paired item counter
    double EntTmpCa[2][3], EntTmpCon[2][3];
    double EntTmp[20];
    int monoFlag, FirthFlag;
    FirthThres = 5;
    Zthres = 0.0;
	ReadParam(argv[1]);
	for (m = 0; m < 23; m++) {
		for (n = 0; n < 10000; n++) {
			hashTable[m][n] = malloc(500 * sizeof(struct data));
			hashLen[m][n] = 0;
        }
	}

    double Ybar, Beta0;
    double tmpSumCa, tmpSumCon;

    for (m = 0; m < size; m++) {
        tmpSumCa = 0;
        tmpSumCon = 0;
        for (n = 0; n < size; n++) {
            tmpSumCa += CaCa[m][n];
            tmpSumCon += ConCon[m][n];
        }
        DefFacDF[0][m] = CaCa[m][m]/tmpSumCa;
        DefFacDF[1][m] = ConCon[m][m]/tmpSumCon;
    }


	FILE *OutFile;
	FILE *LogFile;
	OutFile = fopen(Output,"w");
	LogFile = fopen(strcat(Output,".log"),"w");
	if (OutFile == NULL) {
        printf("Cannot open output file.\n");
        exit(0);
    } //check first wether the output file can be opened
	fprintf(OutFile, "SNP\tCHR\tBP\tA1\tA2\tnCase\tnControl\tOR\tSE\tPval\n");

	char snp[50], a1[100], a2[100];
	int chr;
    char plus = '+', minus = '-', non = '?';
	long int pos, ncase, ncontrol;
	double or, se, freq;
	Data *tmpSNP;
	tmpSNP = malloc(sizeof(Data));
	FILE *InFile;
    int lenCount = 0;

    //read all inputs
    memset(SumZik, 0, sizeof(SumZik[0][0]) * 20 * 20);
    memset(SqrSumZk, 0, sizeof(SqrSumZk[0][0]) * 20 * 20);
    memset(Nik, 0, sizeof(Nik[0][0]) * 20 * 20);
    memset(CorrR, 0, sizeof(CorrR[0][0]) * 20 * 20);

    for (k = 0; k < size; k++) {
        InFile = fopen(Input[k], "r");
        if (InFile == NULL) {
            printf("Cannot open the %ld-th input file.\n", k);
            exit(0);
        }
        else {
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
                    ncase = CaCa[k][k];
                else
                    ncase = atoi(token[nCasec-1]);
                if (!nControlc)
                    ncontrol = ConCon[k][k];
                else
                    ncontrol = atoi(token[nControlc-1]);
                
                tmpSNP = hashSNPsearch(chr, pos, snp);
                if (tmpSNP) {
                    if ((strcmp(a1, (*tmpSNP).Aff) == 0) && (strcmp(a2, (*tmpSNP).Unaff) == 0)) {
                        (*tmpSNP).OR[k] = or;
                        (*tmpSNP).SE[k] = se;
                        (*tmpSNP).nCase[k] = ncase;
                        (*tmpSNP).totalNCase += ncase;
                        (*tmpSNP).nControl[k] = ncontrol;
                        (*tmpSNP).totalNControl += ncontrol;

                        if (Zthres > 0.0)
                            UpdateSumZ((*tmpSNP).OR, (*tmpSNP).SE, k);
                    }
                    else if ((strcmp(a2, (*tmpSNP).Aff) == 0) && (strcmp(a1, (*tmpSNP).Unaff) == 0)) {
                        (*tmpSNP).OR[k] = 1/or;
                        (*tmpSNP).SE[k] = se;
                        (*tmpSNP).nCase[k] = ncase;
                        (*tmpSNP).totalNCase += ncase;
                        (*tmpSNP).nControl[k] = ncontrol;  
                        (*tmpSNP).totalNControl += ncontrol;

                        if (Zthres > 0.0)
                            UpdateSumZ((*tmpSNP).OR, (*tmpSNP).SE, k);
                    }
                    else
                        fprintf(LogFile, "Allele mismatch for SNP %s, a1 = %s, A1 = %s, a2 = %s, A2 = %s.\n", snp, a1, (*tmpSNP).Aff, a2, (*tmpSNP).Unaff);
                }
                else {
                    Data tmp;
                    strcpy(tmp.SNP,snp);
                    strcpy(tmp.Aff,a1);
                    strcpy(tmp.Unaff,a2);
                    tmp.CHR = chr;
                    tmp.Pos = pos;
                    memset(tmp.OR, 0, sizeof(tmp.OR[0]) * 20);
                    memset(tmp.SE, 0, sizeof(tmp.SE[0]) * 20);
                    memset(tmp.nCase, 0, sizeof(tmp.nCase[0]) * 20);
                    memset(tmp.nControl, 0, sizeof(tmp.nControl[0]) * 20);
                    tmp.OR[k] = or;
                    tmp.SE[k] = se;
                    tmp.nCase[k] = ncase;
                    tmp.totalNCase = ncase;
                    tmp.nControl[k] = ncontrol;
                    tmp.totalNControl = ncontrol;
                    tmp.Freq = freq;
                    hashSNPpush(tmp);
                }
            }
        }
        fclose(InFile);
    }

    if (Zthres > 0.0) {
        GetCorrR(size);
        for (m = 0; m < size; m++) {
            tmpSumCa = 0;
            tmpSumCon = 0;
            for (n = 0; n < size; n++) {
                CaCa[m][n] = sqrt(CaCa[m][m]*CaCa[n][n])*CorrR[m][n];
                ConCon[m][n] = sqrt(ConCon[m][m]*ConCon[n][n])*CorrR[m][n];
                tmpSumCa += CaCa[m][n];
                tmpSumCon += ConCon[m][n];
                if (n>m)
                    fprintf(LogFile, "Est. sample overlap for %ld-th input study and %ld-th input study: %lf\n", m+1, n+1, CorrR[m][n]);
            }
            DefFacDF[0][m] = CaCa[m][m]/tmpSumCa;
            DefFacDF[1][m] = ConCon[m][m]/tmpSumCon;
        }
    }


    Data SNP;
    Frq FreqMat[20];
    double tmpShrCa, tmpShrCon;
    Stat SNPout;
    for (m = 0; m < 23; m++) {
        for (n = 0; n < 10000; n++) {
            for (i = 0; i < hashLen[m][n]; i++) {
                SNP = hashTable[m][n][i];
                monoFlag = 1;
                FirthFlag = 1;
                memset(GenoMatOut, 0, sizeof(GenoMatOut[0][0]) * 20 * 6);
                memset(EnCode, 0, sizeof(EnCode[0]) * 20);
                memset(FreqMat, 0, sizeof(FreqMat[0]) * 20);
                for (j = 0; j < size; j++) {
                    FreqMat[j] = GroupFreq(SNP.SE[j], SNP.nCase[j], SNP.nControl[j], SNP.OR[j], SNP.Freq);
                }

                ObsCount(size, FreqMat, SNP.nCase, SNP.nControl);
                IsMonoCode(&monoFlag, size, EnCode);
                IsFirthCorr(&FirthFlag, &Ybar, size, SNP.nCase, SNP.nControl);
                Beta0 = log(Ybar/(1-Ybar));
                if (monoFlag == 0) {
                    SNPout = CountLogit(Beta0, GenoMatOut, EnCode, size, FirthFlag);
                }
                else {
                   SNPout = MonoCodeP(GenoMatOut, size);
                }
                fprintf(OutFile, "%s\t%d\t%ld\t%s\t%s\t%ld\t%ld\t%lf\t%lf\t%.4e\n", SNP.SNP, SNP.CHR, SNP.Pos, SNP.Aff, SNP.Unaff, SNP.totalNCase, SNP.totalNControl, SNPout.OR, SNPout.SE, SNPout.P);
            }
        }
    }
    fclose(OutFile);
    fclose(LogFile);
	return (0);
}
