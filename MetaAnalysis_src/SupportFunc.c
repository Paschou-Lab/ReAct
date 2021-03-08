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



// Mat inverse bundle, scooped from https://github.com/md-akhi/Inverse-matrix/blob/master/Inverse-matrix.c
// Please be correct
// what the hell why does he have to name it as minor
void minor_def(double b[100][100], double a[100][100], int i, int n) {
	int j, l, h = 0, k = 0;
	for (l = 1; l < n;l++) {
		for (j = 0; j < n; j++){
			if (j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if (k == (n-1)){
				h++;
				k = 0;
			}
		}
	}
}

double det(double a[100][100], int n) {
	int i;
	double b[100][100], sum=0;
	if (n == 1)
		return a[0][0];
	else if (n == 2)
		return (a[0][0]*a[1][1] - a[0][1]*a[1][0]);
	else
		for (i = 0; i < n; i++) {
			minor_def(b, a, i, n);
			sum = (double) (sum + a[0][i]*pow(-1,i)*det(b,(n-1)));
		}
	return sum;
}


void transpose(double c[100][100], double d[100][100], double n, double det) {
	int i, j;
	double b[100][100];
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			b[i][j] = c[j][i];
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n;j++)
			d[i][j] = b[i][j]/det;
	}
}

void cofactor(double a[100][100], double d[100][100], double n, double determinte) {
	double b[100][100], c[100][100];
	int l, h, m, k, i, j;
	for (h = 0; h < n; h++)
		for (l = 0; l < n; l++){
			m = 0;
			k = 0;
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					if (i != h && j != l){
						b[m][k] = a[i][j];
						if (k < (n-2))
							k++;
						else{
							k = 0;
							m++;
						}
					}
			c[h][l] = pow(-1,(h+l)) * det(b,(n-1));
		}
	transpose(c, d, n, determinte);
}

void inverse(double a[100][100], double d[100][100], int n, double det) {
	if (det == 0)
		printf("\nInverse of Hessian Matrix is not possible\n");
	else if (n == 1)
		d[0][0] = 1;
	else
		cofactor(a, d, n, det);
}

void cout(double a[100][100],int n,int show) {
	int i,j;
	if(show == 1)
		for(i=0;i < n;i++){
			for(j=0;j < n;j++)
				printf(" %.2f \t",a[i][j]);
			printf("\n");
		}
	else if(show == 2){
		printf("\n\n The Inverse Of Matrix Is : \n\n");
		for (i=0;i<n;i++){
			for (j=0;j<n;j++)
				printf(" %.4f \t",a[i][j]);
			printf("\n");
		}
	}
}

void IsMonoCode(int *flag, int size, double EnCode[20]) {
	int i;
	double tmp;
	for (i = 0; i < size; i++) {
		if (EnCode[i] > MonoThres) {
			tmp = EnCode[i];
			break;
		}
	}
	if (i == (size-1)) 
		*flag = 1;
	else {
		for (i = i+1; i < size; i++) {
			if ((EnCode[i] > MonoThres) &&  (fabs(EnCode[i]-tmp) > MonoThres))  {
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
		if (tmpN > MonoThres) {
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
			if ((tmpN > MonoThres) && (tmpN < minN))
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









