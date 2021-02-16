#ifndef Preparation_H
#define Preparation_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct data {
	char SNP[50], Aff[50], Unaff[50];
	int CHR;
	long int Pos, nCase[2], nControl[2];
	double OR[2], SE[2], Freq;
};
typedef struct data Data;

struct FrqW {
	double w1Ca, w1Con, w2Ca,w2Con;
};
typedef struct FrqW ShrFrqW;

double Zthres;

//mandatory  parameters
char Input[2][5000];
char Output[5000];
long int nCase[2];
long int nControl[2];
long int nCaCa[4];
long int nCaCon[4];
long int nConCon[4];
int size;

double CaCa[2][2];
double CaCon[2][2];
double ConCon[2][2];

// input from the first input file will be stored in a 2d hash table, with a modular hash function over the chromosomal position
Data *hashTable[23][10000];
int hashLen[23][10000];
char buffer[10000];
char bufftmp[10000];


// Hash function to carry out over the base pair position of a SNP
long int hashFunc(long int i);

// Function to push a SNP into the hash table
void hashSNPpush(Data SNP);

// Function to search for a SNP by position and rsID from the hash table and return the pointer to the SNP
Data *hashSNPsearch(long int CHR, long int Pos, char *SNPID);

// Function to read the parameters from the input text file
void ReadParam(char *ParIn);

#endif