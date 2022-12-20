//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _init_H_
#define _init_H_

#include "configuration.h"
#include "GPsDistribution.h"
#include "shapeFunctions.h"
#include "classNodes.h"
#include "classElements.h"
#include "classGPs.h"

extern void initGPsDefinition(vector<classGPs *> &GPs, vector<classNodes *> &nodes, vector<classElements *> &elements, int ndim,
                  int sim, vector<int> &nodesFEM, vector<int> &nodesMM, vector<int> &nodesShare,
                  vector<int> &elementsFEM, vector<int> &elementsMM, vector<constitutiveModels *> &consMod);

extern void initInternalVariables(vector<classGPs *> &GPs, vector<classElements *> &elem);

#endif
