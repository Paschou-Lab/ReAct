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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define iterMax 100 //can change
#define	epsilon 1E-15

double nCa, nCon, nTotal;
double GenoMatOut[20][6]; // max (20+(20*19)/2)*6
double EnCode[20];
double GenoEnCode[3];
double DefFacDF[2][20];
double Const;

double SumZik[20][20];
double SqrSumZk[20][20];
double Nik[20][20];
double CorrR[20][20];

void convertToUpperCase(char *Str);
double WaldP(double beta, double SE);
// double StandBiNormCDF(double rho);
double LogLikelihood(double x2, double y2, double xy, double rho, double N);
void UpdateSumZ(double OR[20], double SE[20], int k);
void GetCorrR(int size);
double LogitPredY(double Eta);
double LogitEta(double *coeff, double *X);
double Ymean(int *Y, long int nObs);
void PrintMat(double a[3][3]);
double MaxNum(double a, double b, double c);
void MatInverse(double b[3][3], double a[3][3]);
void IsMonoCode(int *flag, int size, double EnCode[20]);
void IsFirthCorr(int *flag, double *Ybar, int size, long int nCase[20], long int nControl[20]);

#endif











