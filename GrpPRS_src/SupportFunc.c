#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdflib.h"
#include "SupportFunc.h"


// Hash function to carry out over the base pair position of a SNP
long int hashFunc(long int i) {
	return i%10000;
}

// Function to push a SNP into the hash table
void hashSNPpush(Data SNP) {
	long int ind = hashFunc(SNP.Pos);
	hashTable[SNP.CHR-1][ind][hashLen[SNP.CHR-1][ind]] = SNP;
	hashLen[SNP.CHR-1][ind]++;
}

// Function to search for a SNP by position and rsID from the hash table and return the index off its position or 0
Data *hashSNPsearch(long int CHR, long int Pos, char *SNPID) {
	Data *tmp;
	tmp = malloc(sizeof(Data));
	long int ind = hashFunc(Pos);
	long int i;
	for (i = 0; i < hashLen[CHR-1][ind]; i++) {
		tmp = &hashTable[CHR-1][ind][i];
		if (strcmp(SNPID, (*tmp).SNP) == 0)
			return tmp;
	}
	return NULL;
}

double GetPval(double tStat, double df) {
    double p, q;
    int st = 0; // error variable
    int w = 1; // function variable
    double bnd = 1;
    cdft( &w,&p,&q,&tStat,&df, &st, &bnd);
    return(2*q);
}











