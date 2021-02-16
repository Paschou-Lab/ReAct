#ifndef SupportFunc_H
#define SupportFunc_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "Preparation.h"
#include "CountConstruct.h"
#include "toms462.h"

#define iterMax 15 //can change
#define	MonoThres 1E-15

double nCa, nCon, nTotal;
double GenoMatOut[20][6]; // max (20+(20*19)/2)*6
double EnCode[20];
double DefFacDF[2][20];

double SumZik[20][20];
double SqrSumZk[20][20];
double Nik[20][20];
double CorrR[20][20];


double WaldP(double beta, double SE);
// double StandBiNormCDF(double rho);
double LogLikelihood(double x2, double y2, double xy, double rho, double N);
void UpdateSumZ(double OR[20], double SE[20], int k);
void GetCorrR(int size);
double LogitPredY(double Eta);
double LogitEta(double *coeff, double *X);
double Ymean(int *Y, long int nObs);

void minor_def(double b[100][100], double a[100][100], int i, int n);
double det(double a[100][100], int n);
void transpose(double c[100][100], double d[100][100], double n, double det);
void cofactor(double a[100][100], double d[100][100], double n, double determinte);
void inverse(double a[100][100], double d[100][100], int n, double det);
void cout(double a[100][100],int n,int show);
void IsMonoCode(int *flag, int size, double EnCode[20]);
void IsFirthCorr(int *flag, double *Ybar, int size, long int nCase[20], long int nControl[20]);

#endif











