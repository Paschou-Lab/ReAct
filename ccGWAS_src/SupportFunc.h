#ifndef SupportFunc_H
#define SupportFunc_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cdflib.h"
#include "Preparation.h"
#include "CountConstruct.h"
#include "toms462.h"

#define iterMax 15 //can change
#define	MonoThres 1E-15

double nCa, nCon, nTotal;
double GenoMatOut[2][3]; // max (20+(20*19)/2)*6
double EnCode[2][3];
double DefCase[2];
double DefControl[2];

double SumZ12;
double SqrSumZ1, SqrSumZ2;
double N12;
double CorrR;

void convertToUpperCase(char *Str);
double LogLikelihood(double x2, double y2, double xy, double rho, double N);
void UpdateSumZ(double OR[2], double SE[2]);
double GetCorrR();
double WaldP(double beta, double SE);

#endif











