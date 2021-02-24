#ifndef SupportFunc_H
#define SupportFunc_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "toms462.h"

#define maxDF 20


// Prelim
char TargetIn[maxDF][5000];
long int nCase[maxDF];
long int nControl[maxDF];
double nOverlap[maxDF];
char Output[5000];
char Base[5000];
char buffer[10000];
char bufftmp[10000];
int nInfile;
double thresP;
long int nCaBase, nConBase;
double dfBase;
double Zthres;


// Struct
struct ResTable {
	double res11, res12, res21, res22;
};
typedef struct ResTable Res;

struct GrpFreq {
	double pCa, pCon, pPop;
};
typedef struct GrpFreq Frq;

struct data {
	char SNP[50], Aff[50], Unaff[50];
	int CHR;
	long int Pos;
	double Beta;
	double Pbase, Zbase;
	Frq Basefrq;
};
typedef struct data Data;

Data *hashTable[23][100000];
int hashLen[23][100000];


// Score related
double SumZik[maxDF];
double SqrSumZk[maxDF][2];
double Nik[maxDF];
double CorrR[maxDF];
double q[maxDF];

double BaseScoreCase[maxDF];
double BaseScoreControl[maxDF];
double BaseSDCase[maxDF];
double BaseSDControl[maxDF];
double BaseSDPopS[maxDF];
double BaseVarCaS[maxDF];
double BaseVarConS[maxDF];
double BaseVarPopS[maxDF];

double ScoreCase[maxDF];
double ScoreControl[maxDF];
double ScoreDF[maxDF];
double VarCaS[maxDF];
double VarConS[maxDF];
double VarPopS[maxDF];
double SDCase[maxDF];
double SDControl[maxDF];
double SDDF[maxDF];
double Pval[maxDF];
long int nSNPDf[maxDF];


// Stat related
double df[maxDF];
double TstatObs[maxDF];
double R2Obs[maxDF];
double TstatDis[maxDF];
double R2Dis[maxDF];
double R2Real[maxDF];
double TstatReal[maxDF];


// Support Func
long int hashFunc(long int i);

void hashSNPpush(Data SNP);

Data *hashSNPsearch(long int CHR, long int Pos, char *SNPID);

void UpdateZsum(double Zbase, double or, double se, int k);

void UpdateScore(Data SNP, Frq freq, int k, int flag);

double LogLikelihood(double x2, double y2, double xy, double rho, double N);

void GetCorrR(int size);

double TtoR2(double Tstat, double df);

double R2toT(double R2, double df);

void UpdateObsStat(int k);

void UpdateDisStat(int k);

void UpdateRealStat(int k);

double GetPval(double tStat, double df);

#endif