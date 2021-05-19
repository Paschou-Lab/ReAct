#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "SupportFunc.h"

void fOverlap() {
	double sum = 0.0;
	int i;
	for (i = 0; i < maxDF; i++)
		sum += (nOverlapCa[i] + nOverlapCon[i]);
	if ((sum > 0.0) || (Zthres > 0.0))
		OverlapFlag = 1;
} 

void convertToUpperCase(char *Str) {
  while(*Str) {
    *Str = toupper((unsigned char)*Str);
    Str ++;
    }
}

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

void UpdateZsum(double Zbase, double or, double se, int k) {
	double Zk = (or*se != 0) ? log(or)/se : 0.0f;
    if ((fabs(Zk) < Zthres)&&(fabs(Zbase) < Zthres)) {
        Nik[k] += 1;
        SumZik[k] += Zbase*Zk;
        SqrSumZk[k][0] += pow(Zbase,2);
        SqrSumZk[k][1] += pow(Zk,2);
    }
}

void UpdateScore(Data SNP, Frq freq, int k, int flag) {
	if (flag == 1) {
		ScoreCase[k] += SNP.Beta * freq.pCa;
	    VarCaS[k] += pow(SNP.Beta,2)*freq.pCa*(1-freq.pCa);
		ScoreControl[k] += SNP.Beta * freq.pCon;
	    VarConS[k] += pow(SNP.Beta,2)*freq.pCon*(1-freq.pCon);
		ScoreDF[k] += SNP.Beta * freq.pPop;
	    VarPopS[k] += pow(SNP.Beta,2)*freq.pPop*(1-freq.pPop);
	}
	else if (flag == 0) {
		ScoreCase[k] += SNP.Beta * (1-freq.pCa);
	    VarCaS[k] += pow(SNP.Beta,2)*freq.pCa*(1-freq.pCa);
		ScoreControl[k] += SNP.Beta * (1-freq.pCon);
	    VarConS[k] += pow(SNP.Beta,2)*freq.pCon*(1-freq.pCon);
		ScoreDF[k] += SNP.Beta * (1-freq.pPop);
	    VarPopS[k] += pow(SNP.Beta,2)*freq.pPop*(1-freq.pPop);
	}

    BaseScoreCase[k] += SNP.Beta*SNP.Basefrq.pCa;
    BaseScoreControl[k] += SNP.Beta*SNP.Basefrq.pCon;
    BaseScore[k] += SNP.Beta*SNP.Basefrq.pPop;
    BaseVarCaS[k] += pow(SNP.Beta,2)*SNP.Basefrq.pCa*(1-SNP.Basefrq.pCa);
    BaseVarConS[k] += pow(SNP.Beta,2)*SNP.Basefrq.pCon*(1-SNP.Basefrq.pCon);
    BaseVarPopS[k] += pow(SNP.Beta,2)*SNP.Basefrq.pPop*(1-SNP.Basefrq.pPop);
}

double LogLikelihood(double x2, double y2, double xy, double rho, double N) {
	double l0 = - N/2 * log(1-pow(rho,2)) - (x2 - 2*xy*rho + y2)/(2 * (1 - pow(rho,2)));
	double l1 = N*log(bivnor(-Zthres,-Zthres,rho) + bivnor(Zthres,Zthres,rho) - bivnor(Zthres,-Zthres,rho) - bivnor(-Zthres,Zthres,rho));
	return(l0-l1);
}

void GetCorrR(int size) {
	double rho, tmpLogL, maxLogL, bestRho;
	int k;
	for (k = 0; k < size; k++) {
		bestRho = 0.0;
		maxLogL = LogLikelihood(SqrSumZk[k][0], SqrSumZk[k][1], SumZik[k], bestRho, Nik[k]);
		for (rho = 0.001; rho < 1.0; rho = rho + 0.001) {
			tmpLogL = LogLikelihood(SqrSumZk[k][0], SqrSumZk[k][1], SumZik[k], rho, Nik[k]);
			bestRho = (tmpLogL > maxLogL) ? rho : bestRho;
			maxLogL = (tmpLogL > maxLogL) ? tmpLogL : maxLogL;
		}
		CorrR[k] = bestRho;
	}
}

double TtoR2(double Tstat, double df) {
	double a = pow(Tstat,2)/df;
	double R2 = a/(1+a);
	return(R2);
}

void UpdateObsStat(int k) {
	ScoreCase[k] = ScoreCase[k]/nSNPDf[k];
    SDCase[k] = sqrt(2*VarCaS[k]/pow(2*nSNPDf[k],2));
    ScoreControl[k] = ScoreControl[k]/nSNPDf[k];
    SDControl[k] = sqrt(2*VarConS[k]/pow(2*nSNPDf[k],2));
    ScoreDF[k] = ScoreDF[k]/nSNPDf[k];
    SDDF[k] = sqrt(2*VarPopS[k]/pow(2*nSNPDf[k],2));
    df[k] = nCase[k]+nControl[k]-2;
    TstatObs[k] = fabs(ScoreCase[k] - ScoreControl[k])/(SDDF[k]*sqrt(1/(double)nCase[k] + 1/(double)nControl[k]));
	R2Obs[k] = TtoR2(TstatObs[k], df[k]);
}

void UpdateDisStat(int k) {
	BaseScoreCase[k] = BaseScoreCase[k]/nSNPDf[k];
    BaseSDCase[k] = sqrt(2*BaseVarCaS[k]/pow(2*nSNPDf[k],2));
    BaseScoreControl[k] = BaseScoreControl[k]/nSNPDf[k];
    BaseSDControl[k] = sqrt(2*BaseVarConS[k]/pow(2*nSNPDf[k],2));
	BaseScore[k] = BaseScore[k]/nSNPDf[k];
    BaseSDPopS[k] = sqrt(2*BaseVarPopS[k]/pow(2*nSNPDf[k],2));
}

void ShiftBaseScore(int k) {
	double shift = BaseScore[k]-ScoreDF[k];
	BaseScoreCase[k] -= shift;
	BaseScoreControl[k] -= shift;
}

void UpdateRealStat(int k) {
	double RealScoreCase, RealScoreControl;
	double RealSDDF, RealSDCase, RealSDControl;
	double Realdf;

	double ScoreShareCa, ScoreShareCon;
	RealScoreCase = (ScoreCase[k] - qCa[k]*BaseScoreCase[k])/(1-qCa[k]);
	RealScoreControl = (ScoreControl[k] - qCon[k]*BaseScoreControl[k])/(1-qCon[k]);

	RealSDDF = (pow(SDDF[k],2) - pow(q[k],2)*pow(BaseSDPopS[k],2))/pow((1-q[k]),2);
	RealSDDF = sqrt(RealSDDF);
	RealSDCase = (pow(SDCase[k],2) - pow(qCa[k],2)*pow(BaseSDCase[k],2))/pow((1-qCa[k]),2);
	RealSDCase = sqrt(RealSDCase);
	RealSDControl = (pow(SDControl[k],2) - pow(qCon[k],2)*pow(BaseSDControl[k],2))/pow((1-qCon[k]),2);
	RealSDControl = sqrt(RealSDControl);
	Realdf = (double)nCase[k] * (1-qCa[k]) + (double)nControl[k] * (1-qCon[k]) - 2;
	TstatReal[k] = fabs(RealScoreCase - RealScoreControl)/(RealSDDF*sqrt(1/(double)nCase[k] + 1/(double)nControl[k]));
	R2Real[k] = TtoR2(TstatReal[k], Realdf);

	ScoreCase[k] = RealScoreCase;
	ScoreControl[k] = RealScoreControl;
	SDDF[k] = RealSDDF;
	SDCase[k] = RealSDCase;
	SDControl[k] = RealSDControl;
	df[k] = Realdf;
}

double GetPval(double tStat, double df) {
    double p, q;
    int st = 0; // error variable
    int w = 1; // function variable
    double bnd = 1;
    cdft( &w,&p,&q,&tStat,&df, &st, &bnd);
    return(2*q);
}
