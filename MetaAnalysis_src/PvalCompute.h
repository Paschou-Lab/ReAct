#ifndef PvalCompute_H
#define PvalCompute_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "CountConstruct.h"
#include "SupportFunc.h"


// Hessian matrix using count 
void CountHessian(double H[100][100], double Beta[3], double Count[20][3], double EnCode[20], int size);

// XSz vector using count
void CountXSz(double G[3], double Beta[3], double Count[20][6], double EnCode[20], int size);
// double *CountXSz(double Beta0, double BetaX, double BetaZ,double n[2][3][2]);

// XSz vector using count, penalized (Firth's corretion)
void CountXSzPenalize(double G[3], double Beta[3], double Count[20][6], double EnCode[20], double InvH[100][100], int size);

// Function to compute adjusted standard error from hessian matrix
double adjustSE(double Beta[3], double Count[20][3], double EnCode[20], int size);

// Fake logistic regression with counts (fake IRLS??)
Stat CountLogit(double Beta0, double GenoMatOut[20][6], double EnCode[20], int size, int flag);

// If no population stratification, then use odds ratio test
Stat MonoCodeP(double GenoMatOut[20][6], int size);

#endif