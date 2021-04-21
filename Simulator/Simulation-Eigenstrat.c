#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gen_beta.h"

double dseed;
typedef struct Par {
	int nPop; // max 10
	double Fst;
	int nSample[10];
	long int nSNP;
	int nCausalShr; //For heter only 
	int nCausal[10]; // fist element indicate #causal shared by all populations, 
	// then followed by #causal held by each population individually, starting from the second population
	// so the causal SNPs owned by population 1 is actually shared by all populations
	double risk;
	int Exp;
} param;


double uniformRNG(double *dseed, double Rmin, double Rmax) {
	double d2p31m = 2147483647,
		d2p31  = 2147483711;
	*dseed = 16807*(*dseed) - floor(16807*(*dseed)/d2p31m) * d2p31m;
	double rnd = fabs((*dseed / d2p31));
	return(Rmin + rnd*(Rmax-Rmin)); 
}

// for control sample, risk - 1; for cases, causal allele, risk = r
int SimGeno(double freq, double risk) {
	double rnd = uniformRNG(&dseed, 0.0 , 1.0);
	if (rnd <= freq/(freq+risk*(1-freq)))
		return(1);
	else
		return(2);
}

//Arguements: 1-label(mode), 2-seed, 3-nPop, 4-Fst, 5-nSNP, 6-r, (7-nSample, 8-nCausal)*nPop, (7+2nPop)-Exp identifier, (8+2nPop)-output file name
int main(int argc,char* argv[]) {
	param par;
	if (strcmp(argv[1],"Homo") == 0) {
		par.nCausalShr = 0;
		if (argc == 11) {
			dseed = atof(argv[2]);
			par.nPop = atoi(argv[3]);
			par.Fst = atof(argv[4]);
			par.nSNP = atoi(argv[5]);
			par.risk = atof(argv[6]);
			par.nSample[0] = atoi(argv[7]);
			par.nCausal[0] = atoi(argv[8]);
			par.Exp = atoi(argv[9]);
		}
		else {
			printf("Wrong number of parametrs\n");
			exit(0);
		}
		char tmp1[200], tmp2[200], tmp3[200];

		FILE *PedFile;
		strcpy(tmp1, argv[10]);
		strcpy(tmp2, argv[10]);
		strcpy(tmp3, argv[10]);
		PedFile = fopen(strcat(tmp1,".tped"),"w");
		FILE *FamFile;
		FamFile = fopen(strcat(tmp2,".tfam"),"w");
		FILE *LogFile;
		LogFile = fopen(strcat(tmp3,".log"),"w");

		int geno, j = 0, k = 0;
		double r, freq;
		long int BP = 0, i = 0;
		gen_beta_param betaPar;
		for (i = 0; i < par.nSNP; i++) {
			BP += (int)uniformRNG(&dseed, 1, 1000);
			fprintf(PedFile, "%d rs%d.%ld 0 %ld", par.Exp, par.Exp, i+1, BP); //first four columns, CHR(Exp identified), rsID(Exp.SNPindex), molPos(0), BP;
			double p = uniformRNG(&dseed, 0.1, 0.9);
			double alpha = p*(1-par.Fst)/par.Fst;
			double beta = (1-p)*(1-par.Fst)/par.Fst;
			r = (i < par.nCausal[0] ? par.risk : 1.0); // if causal , else risk = 1.0
			for (j = 0; j < par.nPop; j++) {
				gen_beta_initialize(&betaPar, alpha, beta);
				freq = gen_beta(&betaPar);
				//Simulating genotype for cases first
				for (k = 0; k < par.nSample[0]; k++) {
					fprintf(PedFile, " %d", SimGeno(freq, r));
					fprintf(PedFile, " %d", SimGeno(freq, r));
					if (i == 0)
						fprintf(FamFile, "Pop-%d\tSample-ca%d\t0\t0\t0\t2\n", j+1, k+1);
				}
				//Simulating genotype for controls
				for (k = 0; k < par.nSample[0]; k++) {
					fprintf(PedFile, " %d", SimGeno(freq, 1.0));
					fprintf(PedFile, " %d", SimGeno(freq, 1.0));
					if (i == 0)
						fprintf(FamFile, "Pop-%d\tSample-con%d\t0\t0\t0\t1\n", j+1, k+1);
				}
			}
			fprintf(PedFile, "\n");
		}
		fprintf(LogFile, "%lf\n", dseed);
		fclose(PedFile);
		fclose(FamFile);
		fclose(LogFile);
	}



	if (strcmp(argv[1],"Heter") == 0) {
		dseed = atof(argv[2]);
		par.nPop = atoi(argv[3]);
		int k = 0;
		if (argc == 10+2*par.nPop) {
			par.Fst = atof(argv[4]);
			par.nSNP = atoi(argv[5]);
			par.risk = atof(argv[6]);
			par.nCausalShr = atof(argv[7]);
			for (k = 0; k < par.nPop; k++) {
				par.nSample[k] = atoi(argv[8+2*k]);
				par.nCausal[k] = atoi(argv[9+2*k]);
			}
			par.Exp = atoi(argv[8+2*par.nPop]);
		}
		else {
			printf("Wrong number of parametrs\n");
			exit(0);
		}
		char tmp1[200], tmp2[200], tmp3[200];

		FILE *PedFile;
		strcpy(tmp1, argv[9+2*par.nPop]);
		strcpy(tmp2, argv[9+2*par.nPop]);
		strcpy(tmp3, argv[9+2*par.nPop]);
		PedFile = fopen(strcat(tmp1,".tped"),"w");
		FILE *FamFile;
		FamFile = fopen(strcat(tmp2,".tfam"),"w");
		FILE *LogFile;
		LogFile = fopen(strcat(tmp3,".log"),"w");

		int geno, j = 0;
		double r, freq;
		long int BP = 0, i = 0;
		long int Lower[10];
		Lower[0] = par.nCausalShr;
		for (i = 1; i < par.nPop; i++)
			Lower[i] = Lower[i-1] + par.nCausal[i-1];
		gen_beta_param betaPar;
		for (i = 0; i < par.nSNP; i++) {
			BP += (int)uniformRNG(&dseed, 1, 1000);
			fprintf(PedFile, "%d rs%d.%ld 0 %ld", par.Exp, par.Exp, i+1, BP); //first four columns, CHR(Exp identified), rsID(Exp.SNPindex), molPos(0), BP;
			double p = uniformRNG(&dseed, 0.1, 0.9);
			double alpha = p*(1-par.Fst)/par.Fst;
			double beta = (1-p)*(1-par.Fst)/par.Fst;
			for (j = 0; j < par.nPop; j++) {
				gen_beta_initialize(&betaPar, alpha, beta);
				freq = gen_beta(&betaPar);
				if (i < par.nCausalShr)
					r = par.risk;
				else if (i >= Lower[j] && i < Lower[j]+par.nCausal[j])
					r = par.risk;
				else
					r = 1.0;
				//Simulating genotype for cases first
				for (k = 0; k < par.nSample[j]; k++) {
					fprintf(PedFile, " %d", SimGeno(freq, r));
					fprintf(PedFile, " %d", SimGeno(freq, r));
					if (i == 0)
						fprintf(FamFile, "Pop-%d\tSample-ca%d\t0\t0\t0\t2\n", j+1, k+1);
				}
				//Simulating genotype for controls
				for (k = 0; k < par.nSample[j]; k++) {
					fprintf(PedFile, " %d", SimGeno(freq, 1.0));
					fprintf(PedFile, " %d", SimGeno(freq, 1.0));
					if (i == 0)
						fprintf(FamFile, "Pop-%d\tSample-con%d\t0\t0\t0\t1\n", j+1, k+1);
				}
			}
			fprintf(PedFile, "\n");
		}
		fprintf(LogFile, "%lf\n", dseed);
		fclose(PedFile);
		fclose(FamFile);
		fclose(LogFile);
	}
	return(0);
}





