//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#ifndef _configuration_H_
#define _configuration_H_

class classOutputs;

class classPrintForces;

class classTensor4;

class classTensor6;

class classLargeMatrix;

class classMatrix;

class classCompactMatrix;

class classCRSMatrix;

class classSpatialPoint;

class classNodes;

class classElements;

class classGPs;

class classDirichlet;

class classNeumann;

class classNeumannBCs;

class classNodalForces;

class classTractionRamp;

class classPressureRamp;

class classStochasticPressureRamp;

class classHertzianRamp;

class classStochasticHertzianRamp;

class classPressureInst;

class classCurrentInst;

class classHeatFluxtInst;

class classGravity;

class constitutiveModels;

class solvers;

class classExtraDof;

class classFSI;

#include <math.h>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <list>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <unistd.h>
#include <iomanip>

using namespace std;


#ifdef PARALLEL
#include <mpi.h>
#endif
#define SET_SEED srand48
#define GET_RAND drand48
#define SIZE_NAME_LIB 1020
#define SCALE 0.001
#define RED1 "\033[22;31m"
#define GREEN1 "\033[22;32m"
#define YELLOW1 "\033[01;33m"
#define DEFAULT "\033[0m"

#define WARNING_COLOR YELLOW1
#define INFO_COLOR GREEN1
#define ERROR_COLOR RED1

#define INFO(fmt, args...)                        \
  printf(INFO_COLOR"INFO ::: " fmt DEFAULT"\n",            \
     ##args); fflush(stdout)
#define ERROR(fmt, args...)                        \
  fprintf(stderr,ERROR_COLOR"ERROR: (%s) [%s, Function= %s at Line %d] ::: " fmt DEFAULT"\n", \
      Timestamp(),__BASE_FILE__, __FUNCTION__,__LINE__,##args); fflush(stdout)
#define WARNING(fmt, args...)                        \
  fprintf(stderr,WARNING_COLOR"WARNING ::: " fmt DEFAULT"\n",                \
     ##args); fflush(stdout)

extern char *Timestamp();

#ifdef MM

extern "C" {
void
__maxent_MOD_drivermaxent(int *ndim, int *nsddim, char *scheme, char *priorwt, double xyz[], double p[], double dmax[],
                          double D[], int *maxit, double *eps, int *printflag, int *ierror, double *charlengthscale,
                          double *dettol, double phi[], double dphi[], double ddphi[]);

}

#endif //MM

#endif
