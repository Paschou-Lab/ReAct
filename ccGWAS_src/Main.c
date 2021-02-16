#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "Preparation.h"
#include "CountConstruct.h"
#include "SupportFunc.h"


// Main function
int main(int argc,char* argv[]) {
	long int i=0; // SNP counter
	long int j=0; // others
	long int k=0; // input counter
	long int m=0, n=0; //paired item counter
    Zthres = 0.0;
    
	ReadParam(argv[1]);
	for (m = 0; m < 23; m++) {
		for (n = 0; n < 10000; n++) {
			hashTable[m][n] = malloc(500 * sizeof(struct data));
			hashLen[m][n] = 0;
        }
	}


    double N[3], alpha[3];
    double VarInd, VarReal;

	FILE *OutFile;
	FILE *LogFile;
	OutFile = fopen(Output,"w");
	LogFile = fopen(strcat(Output,".log"),"w");
	if (OutFile == NULL) {
        printf("Cannot open output file.\n");
        exit(0);
    } //check first wether the output file can be opened
	fprintf(OutFile, "SNP\tCHR\tBP\tA1\tA2\tOR\tSE\tPval\tControlOR\tControlSE\n");

	char snp[50], a1[50], a2[50];
	int chr;
	long int pos, ncase, ncontrol, counter;
	double or, se, freq;
	Data *tmpSNP;
	tmpSNP = malloc(sizeof(Data));
	FILE *InFile;

    //read all inputs
    for (k = 0; k < 2; k++) {
        InFile = fopen(Input[k], "r");
        if (InFile == NULL) {
            printf("Cannot open the %ld-th input file.\n", k);
            exit(0);
        }
        else {
            int SNPc = 0, Affc = 0, Unaffc = 0, CHRc = 0, Posc = 0, ORc = 0, SEc = 0, nCasec = 0, nControlc = 0, Frqc = 0;
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
            counter = 0;
            while (fgets(buffer, sizeof(buffer), InFile) != NULL) {
                counter++;
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
                chr = atoi(token[CHRc-1]);
                pos = atoi(token[Posc-1]);
                or = atof(token[ORc-1]);
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

                if (k == 0) {
                    Data tmp;
                    strcpy(tmp.SNP,token[SNPc-1]);
                    strcpy(tmp.Aff,token[Affc-1]);
                    strcpy(tmp.Unaff,token[Unaffc-1]);
                    tmp.CHR = chr;
                    tmp.Pos = pos;
                    tmp.OR[0] = or;
                    tmp.SE[0] = se;
                    tmp.nCase[0] = ncase;
                    tmp.nControl[0] = ncontrol;
                    tmp.Freq = freq;
                    hashSNPpush(tmp);
                }

                else {
                    tmpSNP = hashSNPsearch(chr, pos, snp);    
                    if (tmpSNP) {
                        if ((strcmp(a1, (*tmpSNP).Aff) == 0) && (strcmp(a2, (*tmpSNP).Unaff) == 0)) {
                            (*tmpSNP).OR[k] = or;
                            (*tmpSNP).SE[k] = se;
                            (*tmpSNP).nCase[k] = ncase;
                            (*tmpSNP).nControl[k] = ncontrol;
                            if (Zthres > 0.0)
                                UpdateSumZ((*tmpSNP).OR, (*tmpSNP).SE);
                        }
                        else if ((strcmp(a2, (*tmpSNP).Aff) == 0) && (strcmp(a1, (*tmpSNP).Unaff) == 0)) {
                            (*tmpSNP).OR[k] = 1/or;
                            (*tmpSNP).SE[k] = se;
                            (*tmpSNP).nCase[k] = ncase;
                            (*tmpSNP).nControl[k] = ncontrol;
                            if (Zthres > 0.0)
                                UpdateSumZ((*tmpSNP).OR, (*tmpSNP).SE);
                        }
                        else
                            fprintf(LogFile, "Allele mismatch for SNP %s, a1 = %s, (*tmpSNP).Aff = %s, a2 = %s, (*tmpSNP).Unaff = %s.\n", snp, a1, (*tmpSNP).Aff, a2, (*tmpSNP).Unaff);
                    }
                }
            }
        }
        fclose(InFile);
    }

    if (Zthres > 0.0) {
        CorrR = GetCorrR();
        fprintf(LogFile, "Est. sample overlap: %lf\n", CorrR);
        CaCa[0][1] = CaCa[1][0] = sqrt(CaCa[0][0]*CaCa[1][1])*CorrR;
        ConCon[0][1] = ConCon[1][0] = sqrt(ConCon[0][0]*ConCon[1][1])*CorrR;
    }
    DefCase[0] = CaCa[0][0]/(CaCa[0][0]+CaCa[0][1]);
    DefCase[1] = CaCa[1][1]/(CaCa[1][1]+CaCa[1][0]);
    DefControl[0] = ConCon[0][0]/(ConCon[0][0]+ConCon[0][1]);
    DefControl[1] = ConCon[1][1]/(ConCon[1][1]+ConCon[1][0]);

    Data SNP;
    Frq FreqMat[20];
    double tmpShrCa, tmpShrCon;
    Stat SNPout;
    Stat SNPoutCa, SNPoutCon;
    for (m = 0; m < 23; m++) {
        for (n = 0; n < 10000; n++) {
            for (i = 0; i < hashLen[m][n]; i++) {
                SNP = hashTable[m][n][i];
                memset(FreqMat, 0, sizeof(FreqMat[0]) * 2);
                for (j = 0; j < 2; j++) {
                    FreqMat[j] = GroupFreq(SNP.SE[j], SNP.nCase[j], SNP.nControl[j], SNP.OR[j], SNP.Freq);
                }
                if (FreqNotNan(FreqMat[0]) && FreqNotNan(FreqMat[1])) {
                    SNPoutCon = ORstat(FreqMat[0].pCon, SNP.nControl[0], FreqMat[1].pCon, SNP.nControl[1], DefControl);
                    SNPoutCa = ORstat(FreqMat[0].pCa, SNP.nCase[0], FreqMat[1].pCa, SNP.nCase[1], DefCase);
                    SNPout.OR = exp(log(SNPoutCa.OR)-log(SNPoutCon.OR));
                    SNPout.SE = SNPoutCa.SE;
                    SNPout.P = WaldP(log(SNPout.OR), SNPout.SE);
                    fprintf(OutFile, "%s\t%d\t%ld\t%s\t%s\t%lf\t%lf\t%.4e\t%lf\t%lf\n", SNP.SNP, SNP.CHR, SNP.Pos, SNP.Aff, SNP.Unaff, SNPout.OR, SNPout.SE, SNPout.P, SNPoutCon.OR, SNPoutCon.SE);
                }
            }
        }
    }
    fclose(OutFile);
    fclose(LogFile);
	return (0);
}
