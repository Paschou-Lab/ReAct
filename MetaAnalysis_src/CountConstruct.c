#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CountConstruct.h"


// Function to reconstruct the allele count from each summary statistics
Frq GroupFreq(double se,long int nCase,long int nControl,double or,double freq) {
	Res res, TmpRes1, TmpRes2;
	Frq resFrq;
	double tmp2;
	double x = (double)2*nCase;
	double y = (double)2*nControl;
	double w = pow(se,2);
	double z = or;
	int i;

	for (i = 0; i < 50; i++) { //inflate w(se) by 1.001, for at max 49 times (usually can be done within 5 iterations. maximum inflation is 0.050195, ~5%), if tmp2 still < 0 then give up
		w = w*pow(1.001,i);
		tmp2 = pow((2*z*y*(1-z)-w*x*y*z),2) - 4*y*z*(y*z+x)*(pow((z-1),2) + w*x*z);
		if (tmp2 >= 0)
			break;
	} //this magical for loop deals with no solution due to discriminant < 0
	if (tmp2 < 0) {
		resFrq.pCa = 0.0f;
		resFrq.pCon = 0.0f;
		resFrq.pPop = 0.0f;
	} //If this happens, it means the missing rate of this SNP is more than 20%. need to check the input quality
	else {
		double tmp1 = w*x*y*z - 2*y*z*(1-z);
		tmp2 = pow(tmp2,0.5);
		double tmp3 = 2*(pow((z-1),2) + w*x*z);	
		double d1 = (tmp1-tmp2)/tmp3;
		double c1 = y-d1;
		double b1 = x*d1/(z*y-z*d1+d1);
		double a1 = x-b1;
		double frq1 = c1/(c1+d1);
		double d2 = (tmp1+tmp2)/tmp3;
		double c2 = y-d2;
		double b2 = x*d2/(z*y-z*d2+d2);
		double a2 = x-b2;
		double frq2 = c2/(c2+d2);
		int flag1 = 0;
		int flag2 = 0;

		if ((a1>0) && (b1>0) && (c1>0) && (d1>0)) {
			flag1 = 1;	
			TmpRes1.res11 = a1;
			TmpRes1.res12 = b1;
			TmpRes1.res21 = c1;
			TmpRes1.res22 = d1;
		}
		if ((a2>0) && (b2>0) && (c2>0) && (d2>0)) {
            flag2 = 1;
            TmpRes2.res11 = a2;
            TmpRes2.res12 = b2;
            TmpRes2.res21 = c2;
            TmpRes2.res22 = d2;
        }

		if (flag1 && (!flag2))
			res = TmpRes1;
		else if ((!flag1) && flag2)
			res = TmpRes2;
		else if (flag1 && flag2) {
			if (fabs(freq-frq1) < fabs(freq-frq2))
				res = TmpRes1;
			else
				res = TmpRes2;
		}
		else
			res.res11 = res.res12 = res.res21 = res.res22 = 0.0f;

		resFrq.pCa = res.res11/(res.res11+res.res12);
		resFrq.pCon = res.res21/(res.res21+res.res22);
		resFrq.pPop = (res.res11+res.res21)/(res.res11+res.res12+res.res21+res.res22);
		resFrq.pCa = isnan(resFrq.pCa) ? 0.0f : resFrq.pCa;
		resFrq.pCon = isnan(resFrq.pCon) ? 0.0f : resFrq.pCon;
		resFrq.pPop = isnan(resFrq.pPop) ? 0.0f : resFrq.pPop;
	}
	return resFrq;
}

//order: Case AA0/Aa1/aa2, Control AA0/Aa1/aa2
void GCount(double Gcount[6], double dfCase, double dfControl, Frq dfFrq) {
	if (fabs(dfFrq.pPop) > epsilon) {
		Gcount[0] = dfCase*pow((1-dfFrq.pCa),2);
    	Gcount[1] = dfCase*2*dfFrq.pCa*(1-dfFrq.pCa);
    	Gcount[2] = dfCase*pow(dfFrq.pCa,2);
		Gcount[3] = dfControl*pow((1-dfFrq.pCon),2);
		Gcount[4] = dfControl*2*dfFrq.pCon*(1-dfFrq.pCon);
		Gcount[5] = dfControl*pow(dfFrq.pCon,2);
	}
	else 
		Gcount[0] = Gcount[1] = Gcount[2] = Gcount[3] = Gcount[4] = Gcount[5] = 0.0;
}

int isEmpty(double Gcount[6]) {
	if ((Gcount[0] + Gcount[1] + Gcount[2] + Gcount[3] + Gcount[4] + Gcount[5]) == 0.0)
		return(1);
	else
		return(0);
}

void GCountSub(double Gcount[6], double Gshr[6], double scaler) {
	int i;
	for (i = 0; i < 6; i++)
		Gcount[i] -= scaler*Gshr[i];
}

void GCountPrint(){
	int i,m,k;
	for (i = 0; i < size; i++) {
		printf("Df %d: ", i);
		for (m = 0; m < 2; m++) {
			for (k = 0; k < 3; k++) {
				printf("%lf\t",GenoMatOut[i][3*m+k]);				
			}
		}
		printf("Encode = %lf\n", EnCode[i]);
	}
}

//Y: 0-control, 1-case; X: 0,1,2-genotype; Z: 0-df2, 1-df1, 2-overlap
void ObsCount(int size, Frq FreqMat[20], long int nCase[20], long int nControl[20]) {
	int i, m, k;
	int ExistFlag[20], AllFlag;
	double tmpSumCa, tmpSumCon, minFrq, maxFrq;	
	memset(ExistFlag, 0, sizeof(ExistFlag[0]) * 20);
	AllFlag = 0;
	minFrq = maxFrq = FreqMat[0].pPop;

	for (i = 0; i < size; i++) {
		minFrq = ((FreqMat[i].pPop < minFrq) ? FreqMat[i].pPop : minFrq);
		maxFrq = ((FreqMat[i].pPop > maxFrq) ? FreqMat[i].pPop : maxFrq);
		GCount(GenoMatOut[i], nCase[i], nControl[i], FreqMat[i]);
	}

	for (i = 0; i < size; i++) {
		EnCode[i] = (FreqMat[i].pPop-minFrq)/(maxFrq-minFrq) + 1.0;
	}
	
	for (i = 0; i < size; i++) {
		ExistFlag[i] = 1-isEmpty(GenoMatOut[i]);
		AllFlag += ExistFlag[i];
	}
	if (AllFlag == size) {
		for (i = 0; i < size; i++) {
			for (m = 0; m < 2; m++) {
				for (k = 0; k < 3; k++) {
					GenoMatOut[i][3*m+k] *= DefFacDF[m][i];
				}
			}
		}
	}
	else {
		for (i = 0; i < size; i++) {
			if (ExistFlag[i]) {
		        tmpSumCa = 0;
		        tmpSumCon = 0;
		        for (m = 0; m < size; m++) {
		            tmpSumCa += CaCa[i][m]*ExistFlag[m];
		            tmpSumCon += ConCon[i][m]*ExistFlag[m];
		        }
		        for (k = 0; k < 3; k++) {
					GenoMatOut[i][k] *= CaCa[i][i]/tmpSumCa;
					GenoMatOut[i][k+3] *= ConCon[i][i]/tmpSumCon;
				}
		    }
		}
	}
}
