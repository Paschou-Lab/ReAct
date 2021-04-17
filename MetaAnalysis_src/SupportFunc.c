#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SupportFunc.h"

void convertToUpperCase(char *Str) {
  while(*Str) {
    *Str = toupper((unsigned char)*Str);
    Str ++;
    }
}

double WaldP(double beta, double SE) {
	double Z = fabs(beta/SE);
	double p, q;
	int st = 0; // error variable
	int w = 1; // function variable
	double bnd = 1;
	double mean = 0;
	double sd = 1;
	cdfnor(&w,&p,&q,&Z,&mean,&sd,&st,&bnd);
	return(2*q);
}


double LogLikelihood(double x2, double y2, double xy, double rho, double N) {
	double l0 = - N/2 * log(1-pow(rho,2)) - (x2 - 2*xy*rho + y2)/(2 * (1 - pow(rho,2)));
	double l1 = N*log(bivnor(-Zthres,-Zthres,rho) + bivnor(Zthres,Zthres,rho) - bivnor(Zthres,-Zthres,rho) - bivnor(-Zthres,Zthres,rho));
	return(l0-l1);
}

//Check outside the function that OR[k] != 0.
void UpdateSumZ(double OR[20], double SE[20], int k) {
	int i;
	double Zk = log(OR[k]);
	double Zi;
	Zk = (OR[k]*SE[k] != 0) ? log(OR[k])/SE[k] : 0.0f;
	Zk = (fabs(Zk) <= Zthres) ? Zk : 0.0f;
	for (i = 0; i < k; i++) {
		Zi = (OR[i]*SE[i] != 0) ? log(OR[i])/SE[i] : 0.0f;
		Zi = (fabs(Zi) <= Zthres) ? Zi : 0.0f;
		if (Zi*Zk != 0) {
			Nik[i][k] += 1.0;
			SumZik[i][k] += Zi*Zk;
			SqrSumZk[k][i] += pow(Zk,2);
			SqrSumZk[i][k] += pow(Zi,2);
		}
	}
}

void GetCorrR(int size) {
	long int m,n;
	double rho, tmpLogL, maxLogL, bestRho;
	int i, k;
	for (i = 0; i < size; i++) {
		CorrR[i][i] = 1.0;
		for (k = i+1; k < size; k++) {
			bestRho = 0.0;
			maxLogL = LogLikelihood(SqrSumZk[k][i], SqrSumZk[i][k], SumZik[i][k], bestRho, Nik[i][k]);
			for (rho = 0.001; rho < 1.0; rho = rho + 0.001) {
				tmpLogL = LogLikelihood(SqrSumZk[k][i], SqrSumZk[i][k], SumZik[i][k], rho, Nik[i][k]);
				bestRho = (tmpLogL > maxLogL) ? rho : bestRho;
				maxLogL = (tmpLogL > maxLogL) ? tmpLogL : maxLogL;
			}
			CorrR[i][k] = CorrR[k][i] = bestRho;
		}
	}
}


double LogitPredY(double Eta) {
	int i;
	double Yhat;
	Yhat = 1/(1 + exp(-Eta));
	return (Yhat);
}

double LogitEta(double *coeff, double *X) {
	int i;
	double Eta;
	double tmp = coeff[0];
	for (i = 0; i < size; i++) {
		tmp += coeff[i+1]*X[i];
	}
	Eta = tmp;
	return (Eta);
}

double Ymean(int *Y, long int nObs) {
	long int i;
	double Sum;
	for (i = 0; i < nObs; i++)
		Sum += Y[i];
	double Ybar = (double)Sum/nObs;
	return(Ybar);
}

// The regular or Pseudo Inverse through SVD of a three by three matix -- GLS
void PrintMat(double a[3][3]) {
	int i,j;
	for(i=0;i < 3;i++){
		for(j=0;j < 3;j++)
			printf(" %.17g\t",a[i][j]);
		printf("\n");
	}
}

void MatInverse(double b[3][3], double a[3][3]) {
	int n = 3; int i, j;
	double tmp;
	gsl_matrix * gA = gsl_matrix_alloc(n, n);
	gsl_matrix * inv = gsl_matrix_alloc(n, n);
	gsl_matrix * U = gsl_matrix_alloc(n, n);
	gsl_matrix * V = gsl_matrix_alloc(n, n);
	gsl_vector * S = gsl_vector_alloc(n);
	gsl_matrix * PinvS = gsl_matrix_alloc(n, n);
	gsl_matrix_set_zero(PinvS);
	gsl_matrix * VtSinv = gsl_matrix_alloc(n, n);
	gsl_matrix_set_zero(VtSinv);
	gsl_vector * work = gsl_vector_alloc(n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
   			gsl_matrix_set(gA, i, j, a[i][j]);
   		}
   	}
   	
	gsl_linalg_SV_decomp (gA, V, S, work);
	gsl_vector_free(work);
	gsl_matrix_memcpy (U, gA);

	for (i = 0; i < n; i++) {
		if (fabs(gsl_vector_get(S, 0)) > epsilon)
			tmp = (gsl_vector_get(S, i) > (epsilon * gsl_vector_get(S, 0) * 3) ? 1.0/gsl_vector_get(S, i) : 0.0);
		else
			tmp = 0.0;
		gsl_matrix_set(PinvS, i, i, (fabs(tmp) < epsilon ? 0.0 : tmp));
	}

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, PinvS, 0., VtSinv);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., VtSinv, U, 0., inv);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			b[i][j] = ((fabs(gsl_matrix_get(inv, i, j)) < epsilon) ? 0.0f : gsl_matrix_get(inv, i, j));
	}
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(PinvS);
	gsl_matrix_free(VtSinv);
	gsl_vector_free(S);
	gsl_matrix_free(gA);
}

void IsMonoCode(int *flag, int size, double EnCode[20]) {
	int i;
	double tmp;
	for (i = 0; i < size; i++) {
		if (EnCode[i] > epsilon) {
			tmp = EnCode[i];
			break;
		}
	}
	if (i == (size-1)) 
		*flag = 1;
	else {
		for (i = i+1; i < size; i++) {
			if ((EnCode[i] > epsilon) &&  (fabs(EnCode[i]-tmp) > epsilon))  {
				*flag = 0;
			}
		}
	}
}

void IsFirthCorr(int *flag, double *Ybar, int size, long int nCase[20], long int nControl[20]) {
	long int nCa, nCon, minN, maxN, tmpN;
	int i;
	double r;
	nCa = nCon = 0;
	for (i = 0; i < size; i++) {
		nCa += nCase[i];
		nCon += nControl[i];
		tmpN = nCase[i] + nControl[i];
		if (tmpN > epsilon) {
			minN = maxN = tmpN;
			break;
		}
	}
	if (i == (size-1)) 
		*flag = 0;
	else {
		for (i = i+1; i < size; i++) {
			nCa += nCase[i];
			nCon += nControl[i];
			tmpN = nCase[i] + nControl[i];
			if (tmpN > maxN)
				maxN = tmpN;
			if ((tmpN > epsilon) && (tmpN < minN))
				minN = tmpN;
		}
		r = (nCa < nCon) ? ((double)nCon/(double)nCa) : ((double)nCa/(double)nCon);
		r = (r > ((double)maxN/(double)minN)) ? r : ((double)maxN/(double)minN);
		if (r >= FirthThres)
			*flag = 1;
		else
			*flag = 0;
	}
	*Ybar = (double)nCa/((double)nCa+(double)nCon);
}









