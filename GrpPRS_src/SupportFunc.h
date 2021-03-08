#ifndef SupportFunc_H
#define SupportFunc_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cdflib.h"

#define maxDF 20

struct data {
	char SNP[50], Aff[50], Unaff[50];
	int CHR;
	long int Pos;
	double Beta;
};
typedef struct data Data;

struct ResTable {
	double res11, res12, res21, res22;
};
typedef struct ResTable Res;

struct GrpFreq {
	double pCa, pCon, pPop;
};
typedef struct GrpFreq Frq;

Data *hashTable[23][10000];
int hashLen[23][10000];

double ScoreCase[maxDF];
double ScoreControl[maxDF];
double ScoreDF[maxDF];
double SDCase[maxDF];
double SDControl[maxDF];
double SDDF[maxDF];
double Pval[maxDF];
long int nSNPDf[maxDF];

void convertToUpperCase(char *Str);

long int hashFunc(long int i);

void hashSNPpush(Data SNP);

Data *hashSNPsearch(long int CHR, long int Pos, char *SNPID);

double GetPval(double tStat, double df);

#endif
