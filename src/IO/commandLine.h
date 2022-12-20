//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _commandLine_H_
#define _commandLine_H_

#include <assert.h>
#include "configuration.h"
#include "inputFile.h"
#include "classNodes.h"
#include "classElements.h"

class GeneralOptions
{
    public:
        static std::string workingDirectoryName;
        static std::string inputDirectoryName;
        static std::string outputDirectoryName;
        static std::string inputFileName;
        static std::string petscOptionsFile;
        static int commRank;
        static int commSize;
        static double MMScalingFactor;
        static int MMIntegrationOrder;
        static bool MMAdaptiveSupportRadius;
        static bool noImplicitIterative;
        
    public:
    
        GeneralOptions(){}
        ~GeneralOptions(){}
        static void initialise(int argc, char **argv);
        static void finalise();
        static string getInputFilePath();
        static string getOutputDirectoryPath();
};

extern void inputManagement(vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, vector<int> &nodFEM, vector<int> &nodMM,
                            vector<int> &nodShare, vector<int> &elemFEM, vector<int> &elemMM,
                            vector<classElements *> &geoShape, vector<constitutiveModels *> &consMod,
                            vector<int> &epart, vector<int> &npartShare, vector<int> &nodLocal, vector<int> &elemLocal,
                            vector<int> &npart, vector<solvers *> &solver, vector<double> &initialTemp);

extern void inputManagement(vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, vector<int> &nodFEM, vector<int> &nodMM,
                            vector<int> &nodShare, vector<int> &elemFEM, vector<int> &elemMM,
                            vector<classElements *> &geoShape, vector<constitutiveModels *> &consMod,
                            vector<int> &epart, vector<int> &npartShare, vector<int> &nodLocal, vector<int> &elemLocal,
                            vector<int> &npart, vector<solvers *> &solver, bool flagProc4MM, vector<double> &initialTemp);


#include "petscvec.h"
#include <petscsnes.h>
#include <petscksp.h>
#include "petsc.h"

#ifdef PARALLEL
#include <mpi.h>

//extern void modelInfo(int argc, char **argv, int &nNod, int &nElem, int &nodesPerElement, vector<int> &elementsFEM, vector<int> &nodesFEM, vector<int> &nodesMM);
extern void modelInfo(int &nNod, int &nElem, int &nodesPerElement, vector<int> &elementsFEM,
                      vector<int> &nodesFEM, vector<int> &nodesMM, bool &flagMETIS);

extern void modelInfo(int &nNod, int &nElem, int &nodesPerElement, int &nNodeFEM, int &nElemFEM,
                      vector<int> &elementsFEM, vector<int> &nodesFEM);

#ifdef METIS
extern void metisInput(idx_t *nod, idx_t *elem, bool flagCombined, vector<int> ElemFEM,
                       vector<int> &mapFEM2global, int nNodeFEM);

extern void metisInput(idx_t *nod, idx_t *elem, bool flagCombined, vector<int> ElemFEM,
                       vector<int> &mapFEM2global, int nNodeFEM, vector<int> &mapMM2global, idx_t *elementsMM,
                       idx_t *nodesMM, int nNodeMM);
#endif 
extern void getownership(vector<int> &npart, vector<int> &epart, vector<int> &nodLocal,
                         vector<int> &elemLocal, vector<int> &npartShare);


extern void parallelManagement(int &nNodGlobal, int &nNod, int &nElem, vector<int> &epart,
                               vector<int> &npartShare, vector<int> &nodLocal, vector<int> &elemLocal,
                               vector<int> &npart, bool &flagProc4MM);

void partitionAssignment(vector<int> &npart, vector<int> &epart, vector<int> &arr_flagFM,
                         vector<int> elemFEM);

#endif

#endif
