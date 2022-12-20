//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _boundaryConditions_H_
#define _boundaryConditions_H_

#include "configuration.h"
#include "classGPs.h"
#include <petscsnes.h>
#include <petscksp.h>

extern void NeumannBCManagement(const vector<double> &uK, vector<classGPs *> &GPs, vector<classNodes *> &nodes, vector<classElements *> &element,
                                vector<classNeumannBCs *> &NeumannBCs, int &ndim, vector<double> &force, double timeRun,
                                double dt, bool tangentEstimation, map<pair<int,int>,double>& stiffPart);

extern void deleteNeumannBCs(vector<classNeumannBCs *> &NeumannBCs);

extern void findActiveDof(vector<classDirichlet *> DirichletNodes, int ndim, int numDof, vector<int> &activeDof);

extern void activeSearching(vector<classDirichlet *> DirichletNodes, int ndim, int numDof, vector<int> &activeDof);

#ifdef PARALLEL
extern void findActiveDof(vector<classDirichlet*> DirichletNodes, int ndim, int numDof, vector<int> &activeDof, vector<int> nod_local, vector<int> &nDof_local, int nNod);
#endif
#endif
