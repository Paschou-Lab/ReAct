#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SupportFunc.h"

double LogLikelihood(double x2, double y2, double xy, double rho, double N) {
	double l0 = - N/2 * log(1-pow(rho,2)) - (x2 - 2*xy*rho + y2)/(2 * (1 - pow(rho,2)));
	double l1 = N*log(bivnor(-Zthres,-Zthres,rho) + bivnor(Zthres,Zthres,rho) - bivnor(Zthres,-Zthres,rho) - bivnor(-Zthres,Zthres,rho));
	return(l0-l1);
}

//Check outside the function that OR[k] != 0.
void UpdateSumZ(double OR[2], double SE[2]) {
	double Z1, Z2;
	Z1 = (OR[0]*SE[0] != 0) ? log(OR[0])/SE[0] : 0.0f;
	Z1 = (fabs(Z1) <= Zthres) ? Z1 : 0.0f;
	Z2 = (OR[1]*SE[1] != 0) ? log(OR[1])/SE[1] : 0.0f;
	Z2 = (fabs(Z2) <= Zthres) ? Z2 : 0.0f;
	if (Z1*Z2 != 0) {
		N12 += 1.0;
		SumZ12+= Z1*Z2;
		SqrSumZ1 += pow(Z1,2);
		SqrSumZ2 += pow(Z2,2);
	}
}

double GetCorrR() {
	long int m,n;
	double rho, tmpLogL, maxLogL, bestRho;
	int i, k;
	bestRho = 0.0;
	maxLogL = LogLikelihood(SqrSumZ1, SqrSumZ2, SumZ12, bestRho, N12);
	for (rho = 0.001; rho < 1.0; rho = rho + 0.001) {
		tmpLogL = LogLikelihood(SqrSumZ1, SqrSumZ2, SumZ12, rho, N12);
		bestRho = (tmpLogL > maxLogL) ? rho : bestRho;
		maxLogL = (tmpLogL > maxLogL) ? tmpLogL : maxLogL;
	}
	return(bestRho);
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









