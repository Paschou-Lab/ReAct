#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "SupportFunc.h"


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

double R2toT(double R2, double df) {
	double t = sqrt( R2/(1-R2) * df);
	return(t);
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
    BaseSDPopS[k] = sqrt(2*BaseVarPopS[k]/pow(2*nSNPDf[k],2));
    TstatDis[k] = fabs(BaseScoreCase[k] - BaseScoreControl[k])/(BaseSDPopS[k]*sqrt(1/(double)nCaBase + 1/(double)nConBase));
    R2Dis[k] = TtoR2(TstatDis[k], dfBase);
}

void UpdateRealStat(int k) {
	R2Real[k] = (R2Obs[k] - q[k]*R2Dis[k])/(1-q[k]);
	TstatReal[k] = R2toT(R2Real[k], df[k]);
}

double GetPval(double tStat, double df) {
    double p, q;
    int st = 0; // error variable
    int w = 1; // function variable
    double bnd = 1;
    cdft( &w,&p,&q,&tStat,&df, &st, &bnd);
    return(2*q);
}











