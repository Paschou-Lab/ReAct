#ifndef Prep_H
#define Prep_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SupportFunc.h"
#include "CountConstruct.h"

char TargetIn[maxDF][5000];
long int nCase[maxDF];
long int nControl[maxDF];
char Output[5000];
char Base[5000];
char buffer[10000];
char bufftmp[10000];
int nInfile;
double thresP;

void ReadParam(char *ParIn);

void ReadBase(char *InputBase, FILE *LogFile);

#endif