///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// $  $   $  $$   $$$$  $  $   $ /////////////////////////////
/////////////////////////// $  $   $ $  $  $     $  $$ $$ /////////////////////////////
/////////////////////////// $$$$    $$$$   $$$$  $  $ $ $ /////////////////////////////
/////////////////////////// $        $        $  $  $   $ /////////////////////////////
/////////////////////////// $        $     $$$$  $  $   $ /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


/*
  Authors:
  see Authors.txt
*/

/*!
 * \mainpage MuPhiSim
 * \brief MuPhiSim: multiphysics simulation platform
 * \copyright (C) 2022  This software is licenced under the terms stipulated in the license.txt file located in its root directory
 */           

/*!
 * \file MuPhiSim.cpp
 * \brief Main file for MuPhiSim.
 */
#include "configuration.h"
#include "commandLine.h"
#include "init.h"
#include "boundaryConditions.h"
#include "solvers.h"
#include "output.h"

/*! \brief Main function of MuPhiSim
 */
int main(int argc, char **argv) {
    
    //
    GeneralOptions::initialise(argc,argv);
    //
    double startTime = getCPUTime();
    double startWall = getTimeOfDay();
    int nNodGlobal = 0;
    int rank = GeneralOptions::commRank;
    int sim = 0;   /* Integer to identify what kind of simulation you are executing. The values are the following: MM=1, FEM=2 and both=3 */
    // These arrays define the finite element discretisation of the domain. The elements will be "integration cells"" for the MM domain
    vector<classNodes *> nod; /*  Nodes in the domain*/
    vector<classElements *> elem; /* Elements in domain*/
    vector<classGPs *> GPs; /* GPs in the whole domain. All GPs have a neighbourhood (for FEM, the neihborghood is composed by the nodes of the element)*/
    int ndim = 0, nmulti = 0, nGPs = 0, nElem = 0, nNod = 0; /* Spatial dimension, number of extra dimension, number of Gauss points, number of elements and number of nodes respectively*/
    vector<classNeumannBCs *> NeumannBCs; /* Neumann (natural) boundary conditions*/
    vector<int> nodFEM; /* Array with nodal labels that build the FEM domain*/
    vector<int> nodMM; /* Array with the nodal labels that build the MM domain*/
    vector<int> nodShare; /* Array with the nodal labels that are the interface of both domains*/
    vector<int> elemFEM; /* Array with the element labels that build the FEM domain*/
    vector<int> elemMM; /* Array with the element labels that build the MM domain. In this case, they become the integration cells of the MM method*/
    int nNodMM; /* Number of nodes of the MM domain*/
    int nNodFEM;/* Number of nodes of the FEM domain*/
    int nDof;/* Number of mechanical degrees of freedom*/
    int nExtraDof;/* Number of extra degrees of freedom*/
    // These values are assigned after the reading of the input file
    vector<classElements *> geoShape; /* These elements are the ones that close the geometric shape. The closed polygon will be built by segments in 2D. In 3D, there are several posibilities. Initially, they will be triangles*/
    vector<constitutiveModels *> consMod;// Constitutive models in the domain
    vector<int> npartShare;//Processor ownership of the nodes in the whole model npartShare[node i]=1 if the nodes is mine or shared with other processor, 0 otherwise (empty for sequential simulations)
    vector<int> nodLocal; //Mapping of local nodes to global nodes numbering nodlocal[local node index]= global node index (empty for sequential simulations)
    vector<int> elemLocal;// Mapping of local elements to global elements numbering -> elemlocal[local element index] = global element index (empty for sequential simulations)
    vector<int> npart;//Processor ownership mapping of the nodes in the whole model npart[node i]=rank (empty for sequential simulations)
    vector<int> epart; //Processor ownership mapping of the elements in the whole model epart[element i]=rank (empty for sequential simulations)
    vector<solvers *> solver;  // Array for all solvers
    vector<double> iniValueExtraDof; //This vector defines the initial values of extra fields. If the user doesn't define any value, it is initialised as zero
    

#ifdef PARALLEL
    int nprcs = GeneralOptions::commSize;
    bool flagProc4MM = false;

    parallelManagement(nNodGlobal, nNod, nElem, epart, npartShare, nodLocal, elemLocal, npart, flagProc4MM);

    if (rank == 0) {
        INFO("The simulation is being executed in parallel with %i processors", nprcs);
    }

    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
    inputManagement(nod, elem, NeumannBCs, ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, geoShape,
                    consMod, epart, npartShare, nodLocal, elemLocal, npart, solver, flagProc4MM, iniValueExtraDof);
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#else
    INFO("The simulation is being executed sequentially");
    inputManagement(nod, elem, NeumannBCs, ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, geoShape,
                    consMod, epart, npartShare, nodLocal, elemLocal, npart, solver, iniValueExtraDof);
    nNodGlobal = nod.size();
#endif
    // Here we should define the integration cells in case it is not provided in the input file
    nNodMM = nodMM.size();
    nNodFEM = nodFEM.size();
    if (nNodFEM == 0) {
        sim = 1;// Pure MM simulation
    } else if (nNodMM == 0) {
        sim = 2; // Pure FEM simulation
    } else {
        sim = 3; // Both
    }

    if (rank == 0) {
        INFO("Locating the integration points and building the shape functions");
    }
    initGPsDefinition(GPs, nod, elem, ndim, sim, nodFEM, nodMM, nodShare, elemFEM, elemMM, consMod);
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
    for (int i = 1; i < consMod.size(); i++) {
        nmulti += consMod[i]->getNbrDofConsMod();
    }
    nNod = nod.size();
    nElem = elem.size();
    nGPs = GPs.size();
    nDof = nNod * consMod[0]->getNbrDofConsMod(); // Only considering the displacements so far
    nExtraDof = nNod * nmulti; // Considering all of extra variables

    if (rank == 0) {
        INFO("MM ele %lu nod %lu GPs", elemMM.size(), nodMM.size());
        INFO("FEM ele %lu nod %lu", elemFEM.size(), nodFEM.size());
        INFO("Building the shape functions= OK ");
    }

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Wait all the processors to reach this point
#endif
    if (sim == 1 || sim == 3) {
        // The shape functions to extract the real displacements in MM regions
        maxEntSFatNod(nod, nod, ndim, nodMM, false, false);
    }

    // Initialising internal variables. Just once. For the other solvers the previous state is stored within the class
    initInternalVariables(GPs, elem);

    // init u, v, a


    if (rank == 0)
    {
      INFO("done preprocessing before calling solver (Wall %gs, CPU %gs)", getTimeOfDay()-startWall,getCPUTime() - startTime);
    }
    
    vector<double> u(nDof, 0);
    vector<double> v(nDof, 0);
    vector<double> a(nDof, 0);
    vector<double> vv(nExtraDof, 0);
    int ini = 0;
    int finalV = 0;
    if (iniValueExtraDof.size() != 0) 
    { //vv is previously initialised to zero. If the initial value has been defined in the input file, this is changed in vv
        for (int i = 1; i < consMod.size(); i++) {
            nmulti = consMod[i]->getNbrDofConsMod();
            finalV += nNod * nmulti;
            std::fill(vv.begin() + ini, vv.begin() + finalV, iniValueExtraDof[i-1]);
            ini = finalV;
        }
    }
    vector<double> extF(nDof, 0);
    solver[0]->initValues(u, v, a, vv, extF);// This function initialises u, v, a. Even if they do not exist
    for (int i = 0; i < solver.size(); i++) {
        if (rank == 0) {
            INFO("This is solver from %f to %f and type=%s", solver[i]->getT0(), solver[i]->getTf(),
                 solver[i]->getType().c_str());
        }
        solver[i]->calculate(NeumannBCs, GPs, sim, nod, elem, nodFEM, nodMM, nodShare, elemFEM, elemMM, geoShape,
                             consMod, nodLocal, nNodGlobal, nNod);

        solver[i]->getFinalValues(u, v, a, vv, extF);// This function initialises u, v, a. Even if they do not exist
        if (i != solver.size() - 1) {
            if (rank == 0) {
                if ((solver[i + 1]->getT0() - solver[i]->getTf()) > 1E-12) {
                    ERROR("This is not continuity in time between solvers");
                    ERROR("Final time %f, while init time of the next one %f", solver[i]->getTf(),
                          solver[i + 1]->getT0());
                    exit(-1);
                }
            }
            solver[i + 1]->initValues(u, v, a, vv, extF);
        }
    }
    for (int i = 0; i < solver.size(); i++) {
        delete solver[i];
    }

    for (int i = 0; i < GPs.size(); i++) {
        delete GPs[i];
    }
    for (int i = 0; i < consMod.size(); i++) {
        delete consMod[i];
    }
    for (int i = 0; i < elem.size(); i++) {
        delete elem[i];
    }
    for (int i = 0; i < nod.size(); i++) {
        delete nod[i];
    }

    deleteNeumannBCs(NeumannBCs);

    for (int i = 0; i < geoShape.size(); i++) {
        delete geoShape[i];
    }
    if (rank == 0) {
        INFO("The simulation is finished. See the results in the output folder");
    }

    if (rank == 0)
    {
      INFO("Total time consuming: Wall %gs, CPU %gs", getTimeOfDay()-startWall,getCPUTime()-startTime);
    }
    
    classExtractData::finalise();
    GeneralOptions::finalise();
    
    return 0;
}
