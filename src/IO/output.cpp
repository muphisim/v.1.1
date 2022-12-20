
//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file output.cpp
  \brief All functions that are related to the outputs of the simulations should be here
*/
#include "output.h"
#include "commandLine.h"

/*! Constructor for classOutputs. It should not be called.

*/
classOutputs::classOutputs() {
    _numFiles = 100; /* Number of outputs, defaul value*/
};

classPrintForces::classPrintForces() {
    this->_iterationsCompleted = 0;
    this->_firstTimeForcesOutput = true;
};

static int cont = 0;

static void vtkPlotting2D3D(vector<classNodes *>& nod, int nNod, int nElem, vector<classElements *>& elem, int ndim,
                            const vector<double>& ficDisp, const vector<double>& U, const vector<double>& vK, const vector<double>& aK,
                            vector<classGPs *>& GPs, classOutputs *outputs, double dt, double timeRun, int rank,
                            const vector<int> &nodMM, int sim);

/*! \brief This function is called to plot the outputs of the simulation from the solvers
  @param[in] nod Array with all nodes in the domain
  @param[in] elem Array with all elements (or background cells) in the domain
  @param[in] ndim Dimension of the problem
  @param[in] nodFEM Array with the labels that correspond to the FEM nodes in the nodes general array
  @param[in] nodMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[in] nodShare Array with the labels that correspond to the shared nodes in the nodes general array
  @param[in] elemFEM Array with the labels that correspond to the FEM elements in the elements general array
  @param[in] elemMM Array with the labels that correspond to the MM background cells in the elements general array
  @param[in] uK Array with the displacements of all nodes (they can be real (FEM) or fictitious (MM))
  @param[in] sim Type of simulation that is being executed
  @param[in] nDof Number of degrees for freedom that are activated
  @param[in] GPs Array with all GPs in the domain
  @param[in] outputs Object with all information about the outputs
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] rank Label of the processor
*/
extern void outputsManagement(vector<classNodes *> &nod, vector<classElements *> &elem, int &ndim,const vector<int> &nodFEM,
                              const vector<int> &nodMM, const vector<int> &nodShare, const vector<int> &elemFEM, const vector<int> &elemMM,
                              const vector<double> &uK, const vector<double> &vK, const vector<double> &aK, int sim, int nDof,
                              vector<classGPs *>& GPs, classOutputs *outputs, double dt, double timeRun, int rank) {

    const vector<double>& ficDisp = uK;
    const vector<double>& V = vK;
    const vector<double>& A = aK;
    vector<double> U(nDof, 0);   // Non real disp will be real at FEM and unreallistic at MM
    int nNMM = nodMM.size();
    int nNFEM = nodFEM.size();
    int nNod = nod.size();
    int nElem = elem.size();

    realDisplacementManagement(nod, ndim, nodMM, nodFEM, nodShare, ficDisp, sim, nNod, nNMM, nNFEM, U);
    vtkPlotting2D3D(nod, nNod, nElem, elem, ndim, ficDisp, U, V, A, GPs, outputs, dt, timeRun, rank, nodMM, sim);

}

extern void outputsManagement(vector<classNodes *> &nod, vector<classElements *> &elem, int &ndim, const vector<int> &nodFEM,
                              const vector<int> &nodMM, const vector<int> &nodShare, const vector<int> &elemFEM, const vector<int> &elemMM,
                              const vector<double> &uK, const vector<double> &vK, const vector<double> &aK, const vector<double> &vv,
                              const vector<double> &Fv, int sim, int nDof, vector<classGPs *>& GPs, classOutputs *outputs,
                              double dt, double timeRun, int rank) {

    const vector<double>& ficDisp = uK;
    const vector<double>& V = vK;
    const vector<double>& A = aK;
    vector<double> U(nDof + vv.size() + Fv.size(), 0);   // Non real disp will be real at FEM and unreallistic at MM
    vector<double> fic_unknown(nDof + vv.size() + Fv.size(), 0);
    std::copy(ficDisp.begin(), ficDisp.end(), fic_unknown.begin());
    std::copy(vv.begin(), vv.end(), fic_unknown.begin() + nDof);
    std::copy(Fv.begin(), Fv.end(), fic_unknown.begin() + nDof + vv.size());

    int nNMM = nodMM.size();
    int nNFEM = nodFEM.size();
    int nNod = nod.size();
    int nElem = elem.size();
    realDisplacementManagement(nod, ndim, nodMM, nodFEM, nodShare, fic_unknown, sim, nNod, nNMM, nNFEM, U);
    vtkPlotting2D3D(nod, nNod, nElem, elem, ndim, ficDisp, U, V, A, GPs, outputs, dt, timeRun, rank, nodMM , sim);

}

/*! \brief This function is called to plot the equivalent forces requested by the user
@param[in] timeRun double describing current simulation time
@param[in] outputNodes vector containing all nodes of which the output should be written
@param[in] Forces vector containing the force values to be written
@param[in] ndim dimension of simulation
@param[in] direction vector containing the direction (x,y or z) of the nodes to be written
@param[in] firsttime boolean to check whether this is the first time function is called
@param[in] t0 double initial time
@param[in] tf double final time
*/
extern void
writeExternalForces(double timeRun, const vector<int> &outputNodes, const vector<double> &Forces, int ndim, const vector<int> &direction,
                    bool firsttime, double t0, double tf) {

    std::ofstream myfile;
    string name = GeneralOptions::getOutputDirectoryPath()+ "/externalForces" + std::to_string(t0) + "-" + std::to_string(tf) + ".csv";
    if (firsttime) {
        myfile.open(name);
        myfile << "Node,Direction,Force,Time,\n";
    } else {
        myfile.open(name, std::ofstream::out | std::ofstream::app);
    }


    for (int i = 0; i < outputNodes.size(); i++) {
        if (direction[i] >= ndim) {
            ERROR("ERROR: Direction of force requested exceeds the simulation dimension");
            exit(0);
        }
        double currentForce;
        currentForce = Forces[outputNodes[i] + direction[i] * Forces.size() /
                                               ndim]; //Forces.size()/ndim should equal to the number of nodes in the domain.
        string lineOut = std::to_string(outputNodes[i] + 1) + "," + std::to_string(direction[i]) + "," +
                         std::to_string(currentForce) + "," + std::to_string(timeRun) + ",\n";
        myfile << lineOut;
    }
    myfile.close();
};

/*! \brief This function is called to process the forces to be written requested by the user
@param[in] intForce vector containing the internal forces
@param[in] totalDofs int value of number of degrees of freedom
@param[in] dt double timestep
@param[in] inertiaForce vector containing intertial forces
@param[in] type int describing whether simulation is static or dynamic
@param[in] forceOutputs classPrintForces data structure that contains parameters of the output to be written
@param[in] _t0 double is initial time of simulation
@param[in] _tf double final time of simulation
@param[in] timeRun double current time of the simulation
@param[in] _ndim int dimension of the simulation
@param[in] rank int indicating which processor is calling the function
@param[in] _numDof int total number of DOF on the processor
@param[in]_nNodGlobal int number of nodes in total
@param[in] nNodLocal int number of nodes on the processor calling the function
@param[in] nodLocal vector containing the indecies of nodes on the processor
@param[in] nprcs int number of processors working on the simulation
*/
extern void printForcesOut(const vector<double> &intForce, int totalDofs, double dt, const vector<double> &inertiaForce, int type,
                           classPrintForces *&forceOutputs, double _t0, double _tf, double timeRun, int _ndim, int rank,
                           int _numDof, int nNodGlobal, int nNodLocal, const vector<int> &nodLocal, int nprcs) {
    if (forceOutputs != nullptr) {//not mandatory parameter so check if object created
        double period = forceOutputs->getPeriod();
        vector<int> outputNodes = forceOutputs->getNodesToWrite();
        vector<int> directionsToWrite = forceOutputs->getDirectionsToWrite();
        int periodIterationsCompleted = forceOutputs->getIterationsCompleted();
        vector<double> forceOut(totalDofs, 0);
#ifdef PARALLEL
        vector<int> tempV; // temporary vector with mapping to global
        int dofGlobal;
        for (int j = 0; j < _numDof; j++) {
            locDofToGlobDof(nNodGlobal, nNodLocal, j, dofGlobal, nodLocal, false);
            tempV.push_back(dofGlobal);
        }
        // Get total length
        int *lenPerProc;//global array storing lenth of intForce of every process
        if (rank == 0) {
            lenPerProc = new int[nprcs];
        }
        int localSize = intForce.size();
        MPI_Gather(&localSize, 1, MPI_INT, lenPerProc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int globalSize = 0;
        if (rank == 0) {
            for (int i = 0; i < nprcs; i++) {
                globalSize += lenPerProc[i];
            }
        }
        //combine all the local indicies and internal forces into one array on rank 0
        int *globalIndexing;
        double *intForcesArray;
        double *intForcesGlobalUnOrd;
        intForcesArray = new double[intForce.size()];
        int *locToGlobLoc;
        locToGlobLoc = new int[tempV.size()];

        for (int i = 0; i < _numDof; i++) {
            intForcesArray[i] = intForce[i];//copy vector onto array
            locToGlobLoc[i] = tempV[i];
        }
        if (rank == 0) {
            intForcesGlobalUnOrd = new double[globalSize];//(double*)malloc(length*sizeof(double));
            globalIndexing = new int[globalSize];
        }

        // Gathering forces and indexing from all processors
        int *displc;
        if (rank == 0) {
            displc = new int[nprcs];
            displc[0] = 0;
            for (int i = 1; i < nprcs; i++)
                displc[i] = displc[i - 1] + lenPerProc[i - 1];
        }
        MPI_Gatherv(intForcesArray, intForce.size(), MPI_DOUBLE, intForcesGlobalUnOrd, lenPerProc, displc, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);//in rank order
        MPI_Gatherv(locToGlobLoc, tempV.size(), MPI_INT, globalIndexing, lenPerProc, displc, MPI_INT, 0,
                    MPI_COMM_WORLD);

        // Contructing a global internal forces
        vector<double> internalForceGlobal;
        if (rank == 0) {
            internalForceGlobal.resize(totalDofs, 0);
            for (int j = 0; j < globalSize; j++) {
                internalForceGlobal[globalIndexing[j]] += intForcesGlobalUnOrd[j];//bug fix
            }
        }
        if (type == 1) {
            //combine all accelerations into one array on rank 0
            double *inertiaForcesArray;
            double *inertiaForcesGlobalUnOrd;
            inertiaForcesArray = new double[inertiaForce.size()];

            for (int i = 0; i < _numDof; i++) {
                inertiaForcesArray[i] = inertiaForce[i];//copy vector onto array
            }
            if (rank == 0) {
                inertiaForcesGlobalUnOrd = new double[globalSize];//(double*)malloc(length*sizeof(double));
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gatherv(inertiaForcesArray, inertiaForce.size(), MPI_DOUBLE, inertiaForcesGlobalUnOrd, lenPerProc,
                        displc, MPI_DOUBLE, 0, MPI_COMM_WORLD);//in rank order

            if (rank == 0) {
                vector<double> inertiaForceGlobal(totalDofs, 0);
                //INFO("inertiaForce SIZE=%i", inertiaForceGlobal.size());
                for (int j = 0; j < globalSize; j++) {
                    inertiaForceGlobal[globalIndexing[j]] += inertiaForcesGlobalUnOrd[j];//bug fix
                }
                vectorArithmetic(internalForceGlobal, inertiaForceGlobal, forceOut, totalDofs, true);
                delete[] inertiaForcesGlobalUnOrd;
            }
            delete[] inertiaForcesArray;
        } else if (type == 0) {
            forceOut = internalForceGlobal;
        }
        if (rank == 0)//delete common variables
        {
            delete[] lenPerProc;
            delete[] intForcesGlobalUnOrd;
            delete[] globalIndexing;
            delete[] displc;
        }
        delete[] intForcesArray;
        delete[] locToGlobLoc;
        MPI_Barrier(MPI_COMM_WORLD);
#else
        if (type == 1) {
            vectorArithmetic(intForce, inertiaForce, forceOut, totalDofs, true);
        } else if (type == 0) {
            forceOut = intForce;
        }

#endif
        if (rank == 0) {
            if (period < dt) {
                period = dt;
            }
            if (timeRun >= periodIterationsCompleted * period + _t0) {
                writeExternalForces(timeRun, outputNodes, forceOut, _ndim, directionsToWrite,
                                    forceOutputs->isFirstTime(), _t0, _tf);
                forceOutputs->firstTimeComplete();
                forceOutputs->setIterationCompleted();
            }
        }
    }
};


/*! \brief This print all outputs in the corresponding .vtk files labelled with the current time
  @param[in] nod Array with all nodes in the domain
  @param[in] nNod Number of nodes
  @param[in] nElem Number of elements
  @param[in] elem Array with all elements (or background cells) in the domain
  @param[in] ndim Dimension of the problem
  @param[in] ficDisp Array with the displacements of all nodes (they can be real (FEM) or fictitious (MM))
  @param[in] U Real displacements
  @param[in] GPs Array with all GPs in the domain
  @param[in] outputs Object with all information about the outputs
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] rank Label of the processor
*/
static void vtkPlotting2D3D(vector<classNodes *>& nod, int nNod, int nElem, vector<classElements *>& elem, int ndim,
                            const vector<double>& ficDisp, const vector<double>& U, const vector<double>& V, const vector<double>& A,
                            vector<classGPs *>& GPs, classOutputs *outputs, double dt, double timeRun, int rank,
                            const vector<int> &nodMM, int sim) {
                                
    // Print stresses and strains
    int nGPs = GPs.size();
    // int nExtraDof = (U.size()/nNod) - ndim; // case U=[uK;vvK]
    int nExtraDof = GPs[0]->getConstitutiveManager().getNumExtraDofPerNode();
    int nStochasticDof = (GPs[0]->getConstitutiveManager().getNumMechDofPerNode()- nod[0]->getDimension())/nod[0]->getDimension();
    // This is to plot the internal variables at element level. This is dangerous in cases with several constitutive models (with different internal variables). The internal variales that are being plotted are the ones indicated in the input file (at the positions specified). If a consitutive model does not have that internal variable, a zero will printed. The interpretation of these values is under the user responsibility
    const vector<int>& intVars = outputs->getOutVars();
    int nIntVars = intVars.size();
    vector<int> indices(nGPs, 0);
    vector<double> ULong(nNod * 3, 0);
    vector<double> VLong(nNod * 3, 0);
    vector<double> ALong(nNod * 3, 0);
    int iter2 = 0;
    ofstream myfile2;
    string str;
    str.append(GeneralOptions::getOutputDirectoryPath()+"/MuPhiSimOutputs_");
    str.append(to_string(rank));
    str.append("_.");
    str.append(to_string(cont));
    str.append(".vtk");
    cont++;
    
    ofstream outfileSeries;
    ifstream infileSeries;
    string seriesName;
    seriesName.append(GeneralOptions::getOutputDirectoryPath()+"/allOutputs_");
    seriesName.append(to_string(rank));
    seriesName.append(".vtk.series");
    
    infileSeries.open(seriesName.c_str(),std::ios_base::in);
    if (!infileSeries.is_open() || cont == 1)
    {
        outfileSeries.open(seriesName.c_str(), std::ios_base::out);
        outfileSeries << "{"<< endl;
        outfileSeries << "\"file-series-version\" : \"1.0\","<< endl;
        outfileSeries << "\"files\" : ["<< endl;
        outfileSeries << "{ \"name\" : "<<"\""<< "MuPhiSimOutputs_" << rank<< "_."<< cont-1<< ".vtk\", \"time\" :" << timeRun << "}" << endl;
        outfileSeries << "]"<< endl;
        outfileSeries << "}"<< endl;
        outfileSeries.close();
    }
    else
    {
        std::vector<string> allLines;
        string line;
        while (std::getline(infileSeries, line))
        {
            allLines.push_back(line);
        }
        outfileSeries.open(seriesName.c_str(), std::ios_base::out);
        if (allLines.size() <5)
        {
            ERROR("wrong output file");
        }
        else
        {
            for (int i=0; i < allLines.size()-2; i++)
            {
                outfileSeries  << allLines[i];
                if (i == allLines.size()-3)
                {
                    outfileSeries << "," << endl;
                }
                else
                {
                    outfileSeries << endl;
                }
            }
            outfileSeries << "{ \"name\" : "<<"\""<< "MuPhiSimOutputs_" << rank<< "_."<< cont-1<< ".vtk\", \"time\" :" << timeRun << "}" << endl;
            for (int i = allLines.size()-2; i < allLines.size(); i++)
            {
                outfileSeries  << allLines[i] << endl;
            }
        }
        outfileSeries.close();
    };
    
    
    myfile2.open(str.c_str());
    myfile2 << "# vtk DataFile Version 1.0" << endl;
    myfile2 << ndim << "D Unstructured Grid of " << elem[0]->getType() << " elements" << endl;
    myfile2 << "ASCII" << endl;
    myfile2 << endl;
    myfile2 << "DATASET UNSTRUCTURED_GRID" << endl;
    myfile2 << "POINTS " << nNod << " float" << endl;

    for (int i = 0; i < nNod; i++) {
        vector<double> xyz;
        vector<double> xyzLong(3, 0);
        xyz = nod[i]->getXYZ();
        for (int j = 0; j < ndim; j++) {
            xyzLong[j] = xyz[j];
        }
        for (int j = 0; j < 3; j++) {
            myfile2 << xyzLong[j] << " ";
        }
        myfile2 << endl;
    }
    myfile2 << endl;
    if (ndim == 2) {
        int mult;
        if (elem[0]->getType() == "CPS6") {
            mult = 7;
        } else if (elem[0]->getType() == "CPS3") {
            mult = 4;
        } else if (elem[0]->getType() == "CPE4") {
            mult = 5;
        }
        myfile2 << "CELLS " << nElem << " " << nElem * (mult) << endl;
    } else {
        int mult;
        if (elem[0]->getType() == "C3D4") {
            mult = 5;
        } else if (elem[0]->getType() == "C3D10") {
            mult = 11;
        } else if (elem[0]->getType() == "C3D8") {
            mult = 9;
        }

        myfile2 << "CELLS " << nElem << " " << nElem * (mult) << endl;
    }
    vector<int> myNodes;
    int nNod2;
    for (int i = 0; i < nElem; i++) {
        elem[i]->initIntVars(nIntVars);
        myNodes = elem[i]->getMyNodes();
        nNod2 = myNodes.size();
        myfile2 << nNod2 << " ";
        for (int j = 0; j < nNod2; j++) {
            myfile2 << myNodes[j] << " ";
        }
        myfile2 << endl;
    }
    myfile2 << endl;

    myfile2 << "CELL_TYPES " << nElem << endl;
    if (ndim == 2) {
        for (int i = 0; i < nElem; i++) {
            if (elem[0]->getType() == "CPS6") {
                myfile2 << "22" << endl;// for 2d quad triangle
            } else if (elem[0]->getType() == "CPS3") {
                myfile2 << "5" << endl;
            } else if (elem[0]->getType() == "CPE4") {
                myfile2 << "9" << endl;
            }
        }
    } else {
        for (int i = 0; i < nElem; i++) {
            if (elem[0]->getType() == "C3D10") {
                myfile2 << "24" << endl;// for 3d quad tetra
            } else if (elem[0]->getType() == "C3D4") {
                myfile2 << "10" << endl; // Tetra
            } else if (elem[0]->getType() == "C3D8") {
                myfile2 << "12" << endl; // hexa
            }

        }
    }
    myfile2 << endl;

    myfile2 << "POINT_DATA " << nNod << endl;
    myfile2 << "VECTORS Displacement float " << endl;
    for (int i = 0; i < nNod; i++) {
        for (int j = 0; j < ndim; j++) {
            ULong[i + j * nNod] = U[i + j * nNod];
        }
    }
    for (int i = 0; i < nNod; i++) {
        for (int j = 0; j < 3; j++) {
            myfile2 << ULong[i + j * nNod] << " ";
        }
        myfile2 << endl;
    }
    if (sim ==1 || sim == 3)
    {
        myfile2 << "VECTORS DisplacementDifference float " << endl;
        for (int i = 0; i < nNod; i++) {
            for (int j = 0; j < ndim; j++) {
                ULong[i + j * nNod] = U[i + j * nNod]-ficDisp[i + j * nNod];
            }
        }
        for (int i = 0; i < nNod; i++) {
            for (int j = 0; j < 3; j++) {
                myfile2 << ULong[i + j * nNod] << " ";
            }
            myfile2 << endl;
        }
    }
    
    
    if (nStochasticDof >= 1) {
        for (int k = 1; k <= (nStochasticDof); k++) {
            myfile2 << endl;
            myfile2 << "VECTORS StochasticComponent" << k << " float" << endl;
            for (int i = 0; i < nNod; i++) {
                for (int j = 0; j < ndim; j++) {
                    myfile2 << U[nNod * k * ndim + i + j * nNod] << " ";
                    if (ndim == 2) {
                        if (j == 1) {
                            myfile2 << 0 << " ";
                        }
                    }
                }
                myfile2 << endl;
            }
        }
    }
    myfile2 << endl;

    myfile2 << "VECTORS Velocity float " << endl;
    for (int i = 0; i < nNod; i++) {
        for (int j = 0; j < ndim; j++) {
            VLong[i + j * nNod] = V[i + j * nNod];
        }
    }
    for (int i = 0; i < nNod; i++) {
        for (int j = 0; j < 3; j++) {
            myfile2 << VLong[i + j * nNod] << " ";
        }
        myfile2 << endl;
    }
    myfile2 << endl;

    myfile2 << "VECTORS Acceleration float " << endl;
    for (int i = 0; i < nNod; i++) {
        for (int j = 0; j < ndim; j++) {
            ALong[i + j * nNod] = A[i + j * nNod];
        }
    }
    for (int i = 0; i < nNod; i++) {
        for (int j = 0; j < 3; j++) {
            myfile2 << ALong[i + j * nNod] << " ";
        }
        myfile2 << endl;
    }
    myfile2 << endl;
    if (nExtraDof >= 1) {
        for (int k = 1; k <= nExtraDof; k++) {
            myfile2 << "SCALARS ExtraDof" << k << " float" << endl;
            myfile2 << "LOOKUP_TABLE default" << endl;
            for (int i = 0; i < nNod; i++) {
                myfile2 << U[nNod * (ndim + k - 1) + i] << endl;
            }
            myfile2 << endl;
        }
    }
//MM outputs
    int nNMM = nodMM.size();
    if (nNMM > 0) {
        myfile2 << "SCALARS RMAX float" << endl;
        myfile2 << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < nNod; i++) {

            myfile2 << nod[i]->getRmax() << endl;
        }
        myfile2 << endl;

        myfile2 << "SCALARS Neighbour_Number float" << endl;
        myfile2 << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < nNod; i++) {
            vector<int> neigh = nod[i]->getNeighbours();
            myfile2 << neigh.size() << endl;
        }
        myfile2 << endl;
    }

    myfile2 << "CELL_DATA " << nElem << endl;
    vector<double> delta(ndim * ndim, 0);
    for (int i = 0; i < ndim; i++) {
        delta[i * ndim + i] = 1.0;
    }

    vector<int> cellCount(nElem, 0);
    vector<double> Fv(nElem * nExtraDof, 0); // ExtraDof active part

    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) {
        classGPs *GaussPoint = *ite;
        indices[iter2] = GaussPoint->getCell();
        vector<double> F(ndim * ndim, 0);
        GaussPoint->getF(ficDisp, nNod, ndim, F);
        vector<double> Fprev(ndim * ndim, 0);
        vector<double> intVars23;
        vector<double> C(ndim * ndim, 0);
        intVars23 = GaussPoint->getCurrentState().getInternalVariables();
        classStates& currState = GaussPoint->getCurrentState();
        vector<double> &_E = currState.getGL();
        vector<double> &_Cauchy = currState.getCauchy();
        vector<double> &_volumetricStrain = currState.getVolumetricStrain();
        vector<double> &_equivalentStrain = currState.getEquivalentStrain();
        vector<double> &_VMS = currState.getVMS();
        vector<double> values(nIntVars, 0);
        for (int w = 0; w < nIntVars; w++) {
            if (intVars[w] < intVars23.size()) {
                values[w] = intVars23[intVars[w]];
            } else {
                values[w] = 0.;
            }
        }

        elem[GaussPoint->getCell()]->sumIntVars(values);
        elem[GaussPoint->getCell()]->sumCauchy(_Cauchy);
        elem[GaussPoint->getCell()]->sumE(_E);
        elem[GaussPoint->getCell()]->sumVolumetricStrain(_volumetricStrain);
        elem[GaussPoint->getCell()]->sumVMS(_VMS);
        elem[GaussPoint->getCell()]->sumEquivalentStrain(_equivalentStrain);
        if (nExtraDof >= 1) {
            for (int k = 1; k <= nExtraDof; k++) {
                double w = GaussPoint->getW();
                double J0 = GaussPoint->getJ();
                const vector<double>& FCurr = GaussPoint->getCurrentState().getDeformationGradient();
                double J1 = determinantTensor(FCurr, ndim);
                Fv[(GaussPoint->getCell()) + nElem * (k-1)] += w * J0 * J1 * U[nNod * (ndim + nExtraDof) + iter2 + nElem * (k-1)];
            }
        }

        cellCount[GaussPoint->getCell()] += 1;
        iter2++;
    }

    vector<double> presElem(nElem, 0);
    if (ndim == 2) {
        myfile2 << "TENSORS Cauchy float " << endl;
        for (int k = 0; k < nElem; k++) {
            vector<double> Cauchy = elem[k]->getCauchy();
            for (int i = 0; i < Cauchy.size(); i++) {
                Cauchy[i] = Cauchy[i] / cellCount[k];
            }
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (j == 2 || i == 2) {
                        myfile2 << "0 ";
                    } else {
                        myfile2 << Cauchy[i * ndim + j] << " ";
                    }
                }
                myfile2 << endl;
            }
            myfile2 << endl;
            presElem[k] = (Cauchy[0] + Cauchy[3]) / ndim;
            if(nStochasticDof==0) {
                elem[k]->resetCauchy();
            }
        }
        myfile2 << "TENSORS E float " << endl;
        for (int k = 0; k < nElem; k++) {
            vector<double> E = elem[k]->getE();
            for (int i = 0; i < E.size(); i++) {
                E[i] = E[i] / cellCount[k];
            }
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (j == 2 || i == 2) {
                        myfile2 << "0 ";
                    } else {
                        myfile2 << E[i * ndim + j] << " ";
                    }
                }
                myfile2 << endl;
            }

            myfile2 << endl;
            if(nStochasticDof==0) {
                elem[k]->resetE();
            }
        }
    } else {
        myfile2 << "TENSORS Cauchy float " << endl;
        for (int k = 0; k < nElem; k++) {
            vector<double> Cauchy = elem[k]->getCauchy();
            for (int i = 0; i < Cauchy.size(); i++) {
                Cauchy[i] = Cauchy[i] / cellCount[k];
            }
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    myfile2 << Cauchy[i * ndim + j] << " ";
                }
                myfile2 << endl;
            }
            myfile2 << endl;
            presElem[k] = (Cauchy[0] + Cauchy[4] + Cauchy[8]) / ndim;
            if(nStochasticDof==0) {
                elem[k]->resetCauchy();
            }
        }

        myfile2 << "TENSORS E float " << endl;
        for (int k = 0; k < nElem; k++) {
            vector<double> E = elem[k]->getE();
            for (int i = 0; i < E.size(); i++) {
                E[i] = E[i] / cellCount[k];
            }
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    myfile2 << E[i * ndim + j] << " ";
                }
                myfile2 << endl;
            }
            myfile2 << endl;
            if(nStochasticDof==0) {
                elem[k]->resetE();
            }
        }
    }
    if (nStochasticDof >= 1) {
        for (int s = 1; s <= (nStochasticDof); s++) {
            if (ndim == 2) {
                myfile2 << endl;
                myfile2 << "TENSORS CauchyStochasticComponent" << s << " float" << endl;
                for (int k = 0; k < nElem; k++) {
                    vector<double> Cauchy = elem[k]->getCauchy();
                    for (int i = 0; i < Cauchy.size(); i++) {
                        Cauchy[i] = Cauchy[i] / cellCount[k];
                    }
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            if (j == 2 || i == 2) {
                                myfile2 << "0 ";
                            } else {
                                myfile2 << Cauchy[s * ndim * ndim + i * ndim + j] << " ";
                            }
                        }
                        myfile2 << endl;
                    }
                    myfile2 << endl;
                    if (s == nStochasticDof) {
                        elem[k]->resetCauchy();
                    }
                }
                myfile2 << "TENSORS E_stochastic_Component" << s << " float" << endl;
                for (int k = 0; k < nElem; k++) {
                    vector<double> E = elem[k]->getE();
                    for (int i = 0; i < E.size(); i++) {
                        E[i] = E[i] / cellCount[k];
                    }
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            if (j == 2 || i == 2) {
                                myfile2 << "0 ";
                            } else {
                                myfile2 << E[s * ndim * ndim + i * ndim + j] << " ";
                            }
                        }
                        myfile2 << endl;
                    }

                    myfile2 << endl;
                    if (s == nStochasticDof) {
                        elem[k]->resetE();
                    }
                }
            }else {
                myfile2 << "TENSORS CauchyStochasticComponent" << s << " float" << endl;
                for (int k = 0; k < nElem; k++) {
                    vector<double> Cauchy = elem[k]->getCauchy();
                    for (int i = 0; i < Cauchy.size(); i++) {
                        Cauchy[i] = Cauchy[i] / cellCount[k];
                    }
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            myfile2 << Cauchy[s * ndim * ndim + i * ndim + j] << " ";
                        }
                        myfile2 << endl;
                    }
                    myfile2 << endl;

                    if(s==nStochasticDof) {
                        elem[k]->resetCauchy();
                    }
                }

                myfile2 << "TENSORS E_stochastic_Component" << s << " float" << endl;
                for (int k = 0; k < nElem; k++) {
                    vector<double> E = elem[k]->getE();
                    for (int i = 0; i < E.size(); i++) {
                        E[i] = E[i] / cellCount[k];
                    }
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            myfile2 << E[s*ndim*ndim + i * ndim + j] << " ";
                        }
                        myfile2 << endl;
                    }
                    myfile2 << endl;
                    if(s==nStochasticDof) {
                        elem[k]->resetE();
                    }
                }
            }
        }
    }
    myfile2 << "SCALARS VMS float " << endl;
    myfile2 << "LOOKUP_TABLE default" << endl;
    for (int k = 0; k < nElem; k++) {
        vector<double> VMSdet = elem[k]->getVMS();
        double VMSFinal = VMSdet[0]/ cellCount[k];
        myfile2 << VMSFinal << endl;
        if(nStochasticDof==0) {
            elem[k]->resetVMS();
        }
    }
    if (nStochasticDof >= 1) {
        for (int s = 1; s <= (nStochasticDof); s++) {
            myfile2 << "SCALARS VMSStochasticComponent" << s << " float" << endl;
            myfile2 << "LOOKUP_TABLE default" << endl;
            for (int k = 0; k < nElem; k++) {
                vector<double> VMSdet = elem[k]->getVMS();
                double VMSFinal = VMSdet[s] / cellCount[k];
                myfile2 << VMSFinal << endl;
                if (s == nStochasticDof) {
                    elem[k]->resetVMS();
                }
            }
        }
    }
    myfile2 << endl;

    myfile2 << "SCALARS VolumetricStrain float " << endl;
    myfile2 << "LOOKUP_TABLE default" << endl;
    for (int k = 0; k < nElem; k++) {
        vector<double> volumetricStraindet = elem[k]->getVolumetricStrain();
        double volumetricStrainFinal = volumetricStraindet[0]/ cellCount[k];
        myfile2 << volumetricStrainFinal << endl;
        if(0==nStochasticDof) {
            elem[k]->resetVolumetricStrain();
        }
    }
    if (nStochasticDof >= 1) {
        for (int s = 1; s <= (nStochasticDof); s++) {
            myfile2 << "SCALARS VolumetricStrainStochasticComponent" << s << " float" << endl;
            myfile2 << "LOOKUP_TABLE default" << endl;
            for (int k = 0; k < nElem; k++) {
                vector<double> volumetricStraindet = elem[k]->getVolumetricStrain();
                double volumetricStrainFinal = volumetricStraindet[s] / cellCount[k];
                myfile2 << volumetricStrainFinal << endl;
                if (s == nStochasticDof) {
                    elem[k]->resetVolumetricStrain();
                }
            }
        }
    }
    myfile2 << endl;

    myfile2 << "SCALARS EquivalentStrain float " << endl;
    myfile2 << "LOOKUP_TABLE default" << endl;
    for (int k = 0; k < nElem; k++) {
        vector<double> equivalentStraindet = elem[k]->getEquivalentStrain();
        double equivalentStrainFinal = equivalentStraindet[0]/ cellCount[k];
        myfile2 << equivalentStrainFinal << endl;
        if(0==nStochasticDof) {
            elem[k]->resetEquivalentStrain();
        }
    }
    if (nStochasticDof >= 1) {
        for (int s = 1; s <= (nStochasticDof); s++) {
            myfile2 << "SCALARS EquivalentStrainStochasticComponent" << s << " float" << endl;
            myfile2 << "LOOKUP_TABLE default" << endl;
            for (int k = 0; k < nElem; k++) {
                vector<double> equivalentStraindet = elem[k]->getEquivalentStrain();
                double equivalentStrainFinal = equivalentStraindet[s] / cellCount[k];
                myfile2 << equivalentStrainFinal << endl;
                if (s == nStochasticDof) {
                    elem[k]->resetEquivalentStrain();
                }
            }
        }
    }
    myfile2 << endl;
    myfile2 << "SCALARS Pressure float " << endl;
    myfile2 << "LOOKUP_TABLE default" << endl;
    for (int k = 0; k < nElem; k++) {
        myfile2 << presElem[k] << endl;
    }
    myfile2 << endl;

    if (nExtraDof >= 1) {
        for (int k = 1; k <= nExtraDof; k++) {
            myfile2 << "SCALARS ExtraDofFlux" << k << " float" << endl;
            myfile2 << "LOOKUP_TABLE default" << endl;
            for (int i = 0; i < nElem; i++) {
                //double FvFinal=Fv[i]/cellCount[i]; // Current needs to be summed up and not averaged
                double FvFinal = Fv[i + nElem * (k-1)];
                myfile2 << FvFinal << endl;
            }
            myfile2 << endl;
        }
    }

    for (int j = 0; j < nIntVars; j++) {
        myfile2 << "SCALARS IntVars" << intVars[j] + 1 << " float " << endl;
        myfile2 << "LOOKUP_TABLE default" << endl;
        for (int k = 0; k < nElem; k++) {
            vector<double> intVarss = elem[k]->getIntVars();
	    for (int i = 0; i < intVarss.size(); i++){
                intVarss[i] = intVarss[i] / cellCount[k];
            }
            myfile2 << intVarss[j] << endl;
        }
        myfile2 << endl;
    }

    for (int k = 0; k < nElem; k++) {
        elem[k]->resetIntVars();
    }
    myfile2.close();

}

/*! \brief This function organises the calculation of the real displacements of all nodes
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the problem
  @param[in] nodMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[in] nodFEM Array with the labels that correspond to the FEM nodes in the nodes general array
  @param[in] nodShare Array with the labels that correspond to the shared nodes in the nodes general array
  @param[in] ficDisp Array with the displacements of all nodes (they can be real (FEM) or fictitious (MM))
  @param[in] sim Type of simulation that is being executed
  @param[in] nNod Number of nodes
  @param[in] nNMM Number of MM nodes
  @param[in] nNFEM Number of FEM nodes
  @param[out] U Real displacements
*/
extern void realDisplacementManagement(vector<classNodes *> &nod, int &ndim, const vector<int> &nodMM, const vector<int> &nodFEM,
                                       const vector<int> &nodShare, const vector<double> &ficDisp, int sim, int nNod, int nNMM,
                                       int nNFEM, vector<double> &U) {
    int sMM = nodMM.size() * ndim;
    vector<double> UMM(sMM, 0);
    int ind2 = 0;
    U = ficDisp;
    if (sim == 1) {
        realDisplacements(ficDisp, UMM, nNod, nod, nodMM, ndim);
        for (int i = 0; i < nNMM; i++) {
            for (int j = 0; j < ndim; j++) {
                ind2 = nodMM[i];
                U[ind2 + j * nNod] = UMM[i + j * nodMM.size()];
            }
        }
    } else if (sim == 3) {
        for (int i = 0; i < nNFEM; i++) {
            for (int j = 0; j < ndim; j++) {
                ind2 = nodFEM[i];
                U[ind2 + j * nNod] = ficDisp[ind2 + j * nNod];
            }
        }
        realDisplacements(ficDisp, UMM, nNod, nod, nodMM, ndim);
        for (int i = 0; i < nNMM; i++) {
            for (int j = 0; j < ndim; j++) {
                ind2 = nodMM[i];
                U[ind2 + j * nNod] = UMM[i + j * nodMM.size()];
            }
        }
        
    }
}

/*! \brief This function calculates the real displacements of the MM domain
  @param[in] uKK fictitious displacements of the MM nodes
  @param[out] realUMM Real displacements of the MM nodes
  @param[in] numNodes Number of total nodes
  @param[in] nodes Array with all nodes in the domain
  @param[in] nodesMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[in] ndim Dimension of the problem
*/
extern void realDisplacements(const vector<double> &uKK, vector<double> &realUMM, int numNodes, vector<classNodes *> &nodes,
                              const vector<int> &nodesMM, int ndim) {

    int numNodes2;
    numNodes2 = nodesMM.size();
    int ind = 0, ind2 = 0;
    for (int j = 0; j < numNodes2; j++) {
        ind2 = nodesMM[j];
        classNodes *node = nodes[ind2];
        vector<double> phi = node->getPhi();
        vector<int> neighbours = node->getNeighbours();
        int numNeighbours = neighbours.size();
        for (int i = 0; i < numNeighbours; i++) {
            for (int j = 0; j < ndim; j++) {
                realUMM[ind + j * numNodes2] += phi[i] * uKK[neighbours[i] + j * numNodes];
            }
        }
        ind++;
    }
}

