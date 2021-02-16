#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "PvalCompute.h"


// Hessian size: 2+k(k-1)/2
void CountHessian(double H[100][100], double Beta[3], double Count[20][3], double EnCode[20], int size) {
	// printf("elements: %lf,%lf,%lf,%lf,%lf, %lf\n", H[0][0], H[0][1], H[0][2], H[1][1], H[1][2], H[3][3]);
	double Eta[20][3], Yhat, d[20][3];
	int i,k;
	for (i = 0; i < size; i++){ // indicator enumerate
		for (k = 0; k < 3; k++) {
			Eta[i][k] = Beta[0] + k*Beta[1] + EnCode[i]*Beta[2];
			Yhat = LogitPredY(Eta[i][k]);
			d[i][k] = Yhat*(1 - Yhat);
			// printf("Eta = %lf, d = %lf\n", Eta[i][k], d[i][k]);
		}
	}
	// printf("Check-CountHessian\n");
	for (i = 0; i < size; i++) {
		for (k = 0; k < 3; k++) { 
			H[0][0] += Count[i][k]*d[i][k];
			H[0][1] += k*Count[i][k]*d[i][k];
			H[1][0] += k*Count[i][k]*d[i][k];
			H[1][1] += pow(k,2)*Count[i][k]*d[i][k];
			H[0][2] += EnCode[i]*Count[i][k]*d[i][k];
			H[2][0] += EnCode[i]*Count[i][k]*d[i][k];
			H[1][2] += k*EnCode[i]*Count[i][k]*d[i][k];
			H[2][1] += k*EnCode[i]*Count[i][k]*d[i][k];
			H[2][2] += pow(EnCode[i],2)*Count[i][k]*d[i][k];
		}
		// printf("i = %d, encode = %lf\n", i, EnCode[i]);
	}
}


// XSz vector using count
void CountXSz(double G[3], double Beta[3], double Count[20][6], double EnCode[20], int size) {
	double Eta[20][3], Yhat, d[20][3], z[20][6]; // because the WeightedSum function only accept thi dim
	int i,k,y;
	// int len = size+(size*(size-1))/2;
	// printf("BP XSz: G1 = %lf, G2 = %lf, G3 = %lf\n", G[0], G[1], G[2]);
	for (i = 0; i < size; i++) {
		for (k = 0; k < 3; k++)  {// genotype enumerator --case come first in table 
			for (y = 0; y < 2; y++) { //trait enumerator (0,1), 1-case; 0-control
				Eta[i][k] = Beta[0] + k*Beta[1] + EnCode[i]*Beta[2];
				Yhat = LogitPredY(Eta[i][k]);
				d[i][k] = Yhat*(1 - Yhat);
				z[i][k+y*3] = (d[i][k]!=0.0) ? (Eta[i][k] + ((1-y)-Yhat)/d[i][k]) : 0.0;
				// printf("%d-%dth:Eta = %lf, d = %lf, z = %lf\n", i,k,Eta[i][k], d[i][k], z[i][k+y*3]);
			}
		}
	}
	// printf("Check-CountXSz\n");
	for (i = 0; i < size; i++) {
		for (k = 0; k < 3; k++) {
			G[0] += d[i][k] * (z[i][k]*Count[i][k] + z[i][k+3]*Count[i][k+3]);
			G[1] += k*d[i][k] * (z[i][k]*Count[i][k] + z[i][k+3]*Count[i][k+3]);
			G[2] += EnCode[i]*d[i][k] * (z[i][k]*Count[i][k] + z[i][k+3]*Count[i][k+3]);
		}
	}
	// printf("CountXSz done\n");
	// printf("BP XSz: G1 = %lf, G2 = %lf, G3 = %lf\n", G[0], G[1], G[2]);
}


//////////////////////////////////////////
// XSz vector using count, penalized (Firth's corretion)
void CountXSzPenalize(double G[3], double Beta[3], double Count[20][6], double EnCode[20], double InvH[100][100], int size) {
	double Eta[20][3], Yhat, d[20][3], h[20][3], z[20][6]; // because the WeightedSum function only accept thi dim
	int i,k,y;
	// int len = size+(size*(size-1))/2;
	// printf("BP XSz: G1 = %lf, G2 = %lf, G3 = %lf\n", G[0], G[1], G[2]);
	for (i = 0; i < size; i++) {
		for (k = 0; k < 3; k++)  {// genotype enumerator --case come first in table 
			for (y = 0; y < 2; y++) { //trait enumerator (0,1), 1-case; 0-control
				Eta[i][k] = Beta[0] + k*Beta[1] + EnCode[i]*Beta[2];
				Yhat = LogitPredY(Eta[i][k]);
				d[i][k] = Yhat*(1 - Yhat);

				h[i][k] = d[i][k]*( (InvH[0][0] + k*InvH[0][1] + EnCode[i]*InvH[0][2]) +  
					k*(InvH[1][0] + k*InvH[1][1] + EnCode[i]*InvH[1][2]) +
					EnCode[i]*(InvH[2][0] + k*InvH[2][1] + EnCode[i]*InvH[2][2]));

				z[i][k+y*3] = (d[i][k]!=0.0) ? (Eta[i][k] + ((1-y)-Yhat + h[i][k]*(1/2 - Yhat))/d[i][k]) : 0.0;

				// printf("%d-%dth:Eta = %lf, d = %lf, h = %lf, z = %lf\n", i,k,Eta[i][k], d[i][k], h[i][k], z[i][k+y*3]);
			}
		}
	}
	// printf("Check-CountXSz\n");
	for (i = 0; i < size; i++) {
		for (k = 0; k < 3; k++) {
			G[0] += d[i][k] * (z[i][k]*Count[i][k] + z[i][k+3]*Count[i][k+3]);
			G[1] += k*d[i][k] * (z[i][k]*Count[i][k] + z[i][k+3]*Count[i][k+3]);
			G[2] += EnCode[i]*d[i][k] * (z[i][k]*Count[i][k] + z[i][k+3]*Count[i][k+3]);
		}
	}
	// printf("CountXSz done\n");
	// printf("BP XSz: G1 = %lf, G2 = %lf, G3 = %lf\n", G[0], G[1], G[2]);
}
//////////////////////////////////////////



// Function to compute adjusted standard error from hessian matrix
double adjustSE(double Beta[3], double Count[20][3], double EnCode[20], int size) {
	double tmpM[100][100] = {0};
	double tmpInvM[100][100] = {0};
	double D; int i;
	CountHessian(tmpM, Beta, Count, EnCode, size);
	//int n = 1+size+(size*(size-1))/2;
	D = det(tmpM, 3);
	if (D == 0) {
		for (i = 0; i < 3; i++){
			tmpM[i][i] += 0.01*tmpM[i][i];
		}
		D = det(tmpM, 3);
	}
	inverse(tmpM, tmpInvM, 3, D);
	// printf("invH11: %lf\n",tmpInvM[1][1]);
	double SE = sqrt(tmpInvM[1][1]);
	return (SE);
}


// Fake logistic regression with counts (fake IRLS??)
Stat CountLogit(double Beta0, double GenoMatOut[20][6], double EnCode[20], int size, int flag) {
	int i,j,k,iter = 0;
	double tmpCoeff[3], coeff[3];
	memset(coeff, 0, 3 * sizeof(double));
	memset(tmpCoeff, 0, 3 * sizeof(double));
	coeff[0] = Beta0;
	double Dcoeff = 0;
	double H[100][100], InvH[100][100], D, G[3], HCount[20][3];
	Stat res; // struct defines in CountConstruct.h

	for (i = 0; i < size; i++) {
		for (k = 0; k < 3; k++)
			HCount[i][k] = GenoMatOut[i][k] + GenoMatOut[i][k+3];
	}
	// printf("Check1\n");
	while (iter < iterMax) {
		memset(G, 0, 3 * sizeof(double));
		Dcoeff = 0;
		CountHessian(H, coeff, HCount, EnCode, size);
		D = det(H, 3);
		if (D == 0) {
			for (i = 0; i < 3; i++){
				H[i][i] += 0.01*H[i][i];
			}
			D = det(H, 3);
		}
		inverse(H, InvH, 3, D);
		if (flag == 0)
			CountXSz(G, coeff, GenoMatOut, EnCode, size);
		else 
			CountXSzPenalize(G, coeff, GenoMatOut, EnCode, InvH, size);
		
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				tmpCoeff[i] += InvH[i][j]*G[j];
				InvH[i][j] = 0;
				H[i][j] = 0;
			}
			// printf("Check2\n");
			Dcoeff += fabs(tmpCoeff[i]-coeff[i]);
			coeff[i] = tmpCoeff[i];
			tmpCoeff[i] = 0;
		}
		if (Dcoeff < 1e-4) // same convergence condition as plink, different from R glm
			break;
		iter++;
	}
	res.OR = exp(coeff[1]);
	// printf("OR = %lf\n", res.OR);
	res.SE = adjustSE(coeff, HCount, EnCode, size);
	// printf("SE = %lf\n", res.SE);
	res.P = WaldP(coeff[1], res.SE);
	return(res);
}


Stat MonoCodeP(double GenoMatOut[20][6], int size) {
	Stat res;
	int i, k;
	double a = 0 ,b = 0 ,c = 0 ,d = 0; //a:case aff; b:case unaff; c:control aff; d:control unaff
	for (i = 0; i < size; i++) {
			a += (GenoMatOut[i][1] + 2*GenoMatOut[i][2]);
			b += (GenoMatOut[i][1] + 2*GenoMatOut[i][0]);
			c += (GenoMatOut[i][4] + 2*GenoMatOut[i][5]);
			d += (GenoMatOut[i][4] + 2*GenoMatOut[i][3]);
	}
	res.OR = (a*d)/(b*c);
	res.SE = sqrt(1/a + 1/b + 1/c + 1/d);
	res.P = WaldP(log(res.OR), res.SE);
	return(res);
}


// End of Count based Logistic regression 
