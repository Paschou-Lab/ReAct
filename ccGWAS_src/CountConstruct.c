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
	} //If this happens, it means the missing rate of this SNP is more than 20%. need to check the input quality
	else {
		double tmp1 = w*x*y*z - 2*y*z*(1-z);
		tmp2 = pow(tmp2,0.5);
		double tmp3 = 2*(pow((z-1),2) + w*x*z);	
		// printf("tmp1 = %lf, tmp2 = %lf, tmp3 = %lf\n", tmp1, tmp2, tmp3);	
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

int FreqNotNan(Frq GroupFreq) {
	int flag = 1;
	if (GroupFreq.pPop == 0.0f) 
		flag = 0;
	return(flag);
}

Stat ORstat(double Frq1, double N1, double Frq2, double N2, double DefFac[2]) {
	Stat res;
	double a = 0 ,b = 0 ,c = 0 ,d = 0; //a:case aff; b:case unaff; c:control aff; d:control unaff
	a = Frq1 * 2*N1 * DefFac[0];
	b = (1-Frq1) * 2*N1 * DefFac[0];
	c = Frq2 * 2*N2 * DefFac[1];
	d = (1-Frq2) * 2*N2 * DefFac[1];
	res.OR = (a*d)/(b*c);
	res.SE = sqrt(1/a + 1/b + 1/c + 1/d) ;
	res.P = WaldP(log(res.OR), res.SE);
	return(res);
}

