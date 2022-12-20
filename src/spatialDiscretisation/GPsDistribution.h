//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _GPsDistribution_H_
#define _GPsDistribution_H_

#include "configuration.h"
#include "classNodes.h"
#include "classElements.h"
#include "maths.h"
#include "init.h"
#include "classGPs.h"
#include "geometry.h"

extern void
GPsDistributionFEM(vector<classGPs *> &GPs, vector<classNodes *> &nodes, vector<classElements *> &elements, int ndim,
                   vector<int> &nodesFEM, vector<int> &elementsFEM);

extern void
GPsDistributionMM(vector<classGPs *> &GPs, vector<classNodes *> &nodes, vector<classElements *> &elements, int ndim,
                  vector<int> &nodesMM, vector<int> &elementsMM, vector<classGPs *> &mmGPs, int orderd);

//generalised BC GP definition
extern void
defineBCGPs(vector<classGPs *> &GPs, vector<classElements *> &elements, vector<classNodes *> &nodes, int ndim,
            bool flagSub);

extern void GPsLine1DValues(int n, int &nGPs, double xw[][2]);

extern void GPsTriangle2DValues(int n, int &nGPs, double xw[][3]);

extern void GPsTetrahedron3DValues(int n, int &nGPs, double xw[][4]);

extern void GPsSquare2DValues(int n, int &nGPs, double xw[][3]);

extern void GPsHexahedron3DValues(int n, int &nGPs, double xw[][4]);

#endif
