#ifndef CountConstruct_H
#define CountConstruct_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SupportFunc.h"
#include "Preparation.h"


// Statistics to be computed from a Chisq test
struct Output {
	double OR, SE;
	double P;
};
typedef struct Output Stat;

// Structure of a 2 by 2 able
struct ResTable {
	double res11, res12, res21, res22;
};
typedef struct ResTable Res;

struct GrpFreq {
	double pCa, pCon, pPop;
};
typedef struct GrpFreq Frq;


Frq GroupFreq(double se,long int nCase,long int nControl,double or,double freq);

void GCount(double Gcount[6], double dfCase, double dfControl, Frq dfFrq);

int isEmpty(double Gcount[6]);

void GCountSub(double Gcount[6], double Gshr[6], double scaler);

void GCountPrint();

void ObsCount(int size, Frq FreqMat[20], long int nCase[20], long int nControl[20]);

#endif