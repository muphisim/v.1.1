//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _input_H_
#define _input_H_

#include "configuration.h"
#include "classNodes.h"
#include "classNeumannBCs.h"
#include "classElements.h"
#include "constitutiveModels.h"
#include "output.h"
#include "solvers.h"
#include "dataExtraction.h"

extern void heatInstBCManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                 vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &epart,
                                 vector<int> &elemLocal, int surface_order);

extern void dispBCManagement(vector<classNodes *> &nodes, vector<classDirichlet *> &DirichletNodes, int ndim,
                             ifstream &ifs, vector<int> &npart, vector<int> &nodLocal, vector<int> &activeDofMapLocal,
                             vector<int> &nActiveDof, vector<int> &ISactiveDof); //Kamel, vector<int> &nActiveDof instead of int &nActiveDof
extern void dispBCManagement(vector<classNodes *> &nodes, vector<classDirichlet *> &DirichletNodes, int ndim,
                             ifstream &ifs, vector<int> &npart, vector<int> &nodLocal, vector<int> &activeDofMapLocal,
                             vector<int> &nActiveDof, vector<int> &ISactiveDof, double duration,
                             vector<constitutiveModels *> constitutiveModel, bool flagInstantaneous); //Kamel, vector<int> &nActiveDof instead of int &nActiveDof || new variable constitutiveModel
extern void nodesManagement(vector<classNodes *> &nodes, int &ndim, ifstream &ifs, vector<int> npartShare);

extern void elementsManagement(vector<classElements *> &elements, int ndim, ifstream &ifs, vector<classNodes *> &nodes,
                               vector<int> &epart, vector<int> &nodLocal, int &surface_order);

extern void nodesShareManagement(vector<int> &nodesShare, ifstream &ifs, vector<int> &npartShare);

extern void nodesMMManagement(vector<int> &nodesMM, ifstream &ifs, vector<classNodes *> &nodes, vector<int> &npartShare);

extern void nodesMMManagement(vector<int> &nodesMM, ifstream &ifs, vector<classNodes *> &nodes, vector<int> &npartShare,
                              vector<int> nodLocal);

extern void nodesFEMManagement(vector<int> &nodesFEM, ifstream &ifs, vector<classNodes *> &nodes,
                               vector<int> &npartShare);

extern void nodesFEMManagement(vector<int> &nodesFEM, ifstream &ifs, vector<classNodes *> &nodes,
                               vector<int> &npartShare, vector<int> nodLocal);

extern void elementsMMManagement(vector<int> &elementsMM, ifstream &ifs, vector<classElements *> &elements,
                                 vector<int> &epart);

extern void elementsMMManagement(vector<int> &elementsMM, ifstream &ifs, vector<classElements *> &elements,
                                 vector<int> &epart, vector<int> elemLocal);

extern void elementsFEMManagement(vector<int> &elementsFEM, ifstream &ifs, vector<classElements *> &elements,
                                  vector<int> &npartShare);

extern void elementsFEMManagement(vector<int> &elementsFEM, ifstream &ifs, vector<classElements *> &elements,
                                  vector<int> &epart, vector<int> elemLocal);

extern void materialManagement(vector<classElements *> &elements, ifstream &ifs, int ndim,
                               vector<constitutiveModels *> &consMod, vector<int> &epart, vector<int> &elemLocal);

extern void initialConditionsManagement(ifstream &ifs, vector<double> &iniValueExtraDof);

extern void dataExtractionManagement(ifstream &ifs, vector<classElements *> &elements,
                             vector<classNodes *> &nodes, vector<int> &npart,
                             vector<int> &nodLocal, vector<int> &epart, vector<int> &elemLocal);
extern void solverManagement(vector<solvers *> &solver, ifstream &ifs, vector<classElements *> &elements,
                             vector<classNodes *> &nodes, int nDof, int ndim, vector<int> &npart,
                             vector<int> &nodLocal, vector<int> &epart, vector<int> &elemLocal);

extern void NeumannCommonManagement(vector<classElements *> &elements, vector<classNodes *> &nodes,
                                    vector<classNeumannBCs *> &NeumannBCs, int ndim, ifstream &ifs, vector<int> &npart,
                                    vector<int> &nodLocal, vector<int> &epart, vector<int> &elemLocal,
                                    int surface_order);

extern void outManagement(classOutputs *&outputs, ifstream &ifs);

extern void forceOutManagement(classPrintForces *&forceOutputs, ifstream &ifs);

extern void forceBCManagement(vector<classNodes *> &nodes, vector<classNeumannBCs *> &NeumannBCs, int ndim,
                              ifstream &ifs, vector<int> &npart, vector<int> &nodLocal);

extern void FSIManagement(vector<classNodes *> &nodes, vector<classNeumannBCs *> &NeumannBCs, ifstream &ifs, vector<int> npart, vector<int> &nodLocal);

extern void NodeDataExtractionManagement(vector<classNodes *> &nodes, ifstream &ifs, vector<int> npart, vector<int> &nodLocal);

extern void ElementDataExtractionManagement(vector<classElements *> &elements, ifstream &ifs, vector<int> &epart, vector<int> &elemLocal);

#ifdef PARALLEL
#include <mpi.h>
#include "petscvec.h"
#include <petscsnes.h>
#include <petscksp.h>
#include "petsc.h"
#ifdef METIS
#include <metis.h>
extern void metiselementsManagement(idx_t *elements, ifstream &ifs, idx_t *nodes);
extern void metiselementsManagement(idx_t *elements, ifstream &ifs, idx_t *nodes, vector<int> elemFEM, vector<int> &mapFEM2global,int nNodeFEM);
extern void metiselementsManagement(idx_t *elements, ifstream &ifs, idx_t *nodes, vector<int> elemFEM, vector<int> &mapFEM2global, int nNodeFEM, vector<int> &mapMM2global,idx_t *elementsMM, idx_t *nodesMM, int nNodeMM);
void resequence(idx_t *nodes, int nodeNum, int vecLength);
#endif
extern void infoManagement(int &nNod, int &nElem, int &nodesPerElement, ifstream &ifs);
extern void infoManagement(int &nNod, int &nElem, int &nodesPerElement, ifstream &ifs, vector<int> &Elems, vector<int> &nodes);
extern void nodeOwnership(vector<int> &epart, vector<int> &npartShare, vector<int> &nodLocal, ifstream &ifs);
extern void FEM_node_elements(int &nNod, int &nElem, ifstream &ifs, vector<int> &elementsFEM, vector<int> &nodesFEM, vector<int> &allElems, vector<int> &allNodes);
extern void FEM_MM_node_elements(ifstream &ifs, vector<int> &elementsFEM, vector<int> &nodesFEM, vector<int> &allElems, vector<int> &allNodes,vector<int> &nodesMM);
void readPredef_part(ifstream &ifs,vector<int> &npart, vector<int> &epart);
bool meshPartitionMethod(ifstream &ifs);
#endif

#endif
