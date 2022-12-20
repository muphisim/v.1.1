//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _supportDomains_H_
#define _supportDomains_H_

#include "configuration.h"
#include "classNodes.h"
#include "classElements.h"
#include "classGPs.h"
#include "geometry.h"

extern void supportDomain(double *nodesXYZ, int ndim, int numNodes, vector<double> &rmax, double &charlengthscale, int S0);

extern void findNeighbours(vector<double> GP, int me, double *nodesXYZ, int ndim, int numNodes, vector<int> &neighbours,
                           const vector<double> &rmax);

#endif
