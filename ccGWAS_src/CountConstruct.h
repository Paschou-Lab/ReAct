#ifndef CountConstruct_H
#define CountConstruct_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SupportFunc.h"


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

int FreqNotNan(Frq GroupFreq);

Stat ORstat(double Frq1, double N1, double Frq2, double N2, double DefFac[2]);


#endif