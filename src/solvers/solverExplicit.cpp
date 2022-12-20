//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file solverExplicit.cpp
  \brief This file contains all functions related to the explicit solver
*/
#include "solverExplicit.h"
#include "output.h"
#include "boundaryConditions.h"
#include "FSIAdapter.h"
#include "dataExtraction.h"
#include "nonlinearSystem.h"


/*! \brief Constructor of the class
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of degrees of fredom
  @param[in] DNod Dirichlet nodes
  @param[in] scaleFactor Value to be multiplied by the critical time step (<1 for Explicit, >= 1 Implicit)
  @param[in] outputs Object with all data needed for the outputs
  @param[in] t0 Initial time for this solver
  @param[in] tf Final for this solver
  @param[in] solverType String to identify the type of solver
  @param[in] activeDofMapLocal Array for parallel simulations with the local active Dof of a processor (empty in sequential simulations)
  @param[in] nActiveDof Global number of active Dof in parallel simulations (zero in sequential simulations)
  @param[in] ISactiveDof Array for parallel simulations with the mapping for the global displacements vector (empty in sequential simulations)
*/
classSolverExplicit::classSolverExplicit(vector<classNodes *> &nod, int ndim, int numDof,
                                         vector<classDirichlet *> &DNod, double scaleFactor, classOutputs *&outputs,
                                         double t0, double tf, string solverType,  const vector<int> &activeDofMapLocal,
                                         const vector<int>& nActiveDof, const vector<int> &ISactiveDof, int flag_AV,
                                         const vector<double>& para_AV, classPrintForces *&forceOutputs, bool flagActivation,
                                         string &petsc_solver, string &petsc_precon, bool flagFSI)
        : solvers(nod, ndim, numDof, DNod, scaleFactor, outputs, t0, tf, solverType, activeDofMapLocal, nActiveDof,
                  ISactiveDof, forceOutputs, flagActivation, petsc_solver, petsc_precon, flagFSI) {
    _flagAV = flag_AV;
    _c1 = para_AV[0];
    _cl = para_AV[1];
}

/*! \brief This function executes the solver. It MUST be implemented in all new solvers added
  @param[in] NeumannBCs Array with all Neumann BCs. They are global for all solvers
  @param[in] GPs Array with all GPs in the domain
  @param[in] sim Type of simulation executed (see the main for further information)
  @param[in] nod Array with all nodes in the domain
  @param[in] elem Array with all elements (or background cells) in the domain
  @param[in] nodFEM Array with the labels that correspond to the FEM nodes in the nodes general array
  @param[in] nodMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[in] nodShare Array with the labels that correspond to the shared nodes in the nodes general array
  @param[in] elemFEM Array with the labels that correspond to the FEM elements in the elements general array
  @param[in] elemMM Array with the labels that correspond to the MM background cells in the elements general array
  @param[in] geoShape Array with all elements that close the boundary of the domain. These elements will be defined in a dimension ndim-1
  @param[in] consMod Array with all constitutive models in the domain
  @param[in] nodLocal Index mapping from local to global node in parallel simulations (empty for sequential)
  @param[in] nNodGlobal Number of nodes of the global-whole model in parallel simulations (zero for sequential)
  @param[in] nNodLocal Local number of nodes in parallel simulations (zero for sequential)
*/
void classSolverExplicit::calculate(vector<classNeumannBCs *> &NeumannBCs, vector<classGPs *> &GPs, int sim,
                                    vector<classNodes *> &nod, vector<classElements *> &elem, const vector<int> &nodFEM,
                                    const vector<int> &nodMM, const vector<int> &nodShare, const vector<int> &elemFEM,
                                    const vector<int> &elemMM, vector<classElements *> &geoShape,
                                    vector<constitutiveModels *> &consMod, const vector<int> &nodLocal, int nNodGlobal,
                                    int nNodLocal) {

    double startTime = getCPUTime();
    double startWall = getTimeOfDay();
    
    int numMechDofPerNode = consMod[0]->getNbrDofConsMod();
    int numEtraDofPerNode = 0;
    for (int i=1; i< consMod.size(); i++)
    {
        numEtraDofPerNode += consMod[i]->getNbrDofConsMod();
    }
    int numExtraDof = _numNodes*numEtraDofPerNode;
    int totalNumberDofsPerNode = numMechDofPerNode+numEtraDofPerNode;
    //
    bool extraDofFlag(false), StochasticFlag(false);
    int activeNumberDispDofs(0), activeNumberStochasticDofs(0), activeNumberExtraDofs(0);
    if (numMechDofPerNode == _ndim)
    {
        StochasticFlag = false;
        activeNumberDispDofs = _nActiveDof[0];
        activeNumberStochasticDofs = 0;
        activeNumberExtraDofs = _nActiveDof[1];
        if (numEtraDofPerNode>0)
            extraDofFlag = true;
    }
    else
    {
        StochasticFlag = true;
        extraDofFlag = false;
        activeNumberDispDofs = _nActiveDof[0];
        activeNumberStochasticDofs = _nActiveDof[1];
        activeNumberExtraDofs = 0;
    }
    int numGPs = GPs.size();
    vector<double> extraDofAdditionalOutput(numGPs *numEtraDofPerNode, 0);
    vector<double> allField(_numNodes*totalNumberDofsPerNode, 0);
    vector<double> allAcc(_numNodes*totalNumberDofsPerNode, 0); // all acceleration
    vector<double> extForce(_numNodes*totalNumberDofsPerNode, 0);
    vector<double> intForce(_numNodes*totalNumberDofsPerNode, 0);
    vector<double> inerForce(_numNodes*totalNumberDofsPerNode, 0);
    //temporary data for view
    vector<double> uK(_numNodes*numMechDofPerNode, 0);
    vector<double> vK(_numNodes*numMechDofPerNode, 0);
    vector<double> aK(_numNodes*numMechDofPerNode, 0);
    vector<double> vvK(_numNodes*numEtraDofPerNode, 0);
    //
    vector<int> scatter_all(_numNodes*totalNumberDofsPerNode,0);
    vector<int> scatter_uK(_numNodes*numMechDofPerNode,0);
    vector<int> scatter_vvK(_numNodes*numEtraDofPerNode,0);
    for (int j = 0; j < _numNodes*numMechDofPerNode; j++) 
    {
        int dofGlobal=j;
#ifdef PARALLEL  
        locDofToGlobDof(nNodGlobal, nNodLocal, j, dofGlobal, nodLocal, false);
#endif
        scatter_uK[j] = (dofGlobal);
        scatter_all[j] = dofGlobal;
    }
    for (int j = _numNodes*numMechDofPerNode; j < _numNodes*numMechDofPerNode+numExtraDof; j++) 
    {
        int dofGlobal = j;
#ifdef PARALLEL  
        locDofToGlobDof(nNodGlobal, nNodLocal, j, dofGlobal, nodLocal, false);
#endif     
        scatter_vvK[j-_numNodes*numMechDofPerNode] = (dofGlobal);
        scatter_all[j] = dofGlobal;
    };
    
    // create system
    nonlinearSystem* nlsys = new ExplicitNonlinearSystem();
    // allocate solution
    
    int totalDofs = nNodGlobal*totalNumberDofsPerNode;
    nlsys->allocateSolution(totalDofs,scatter_all);
    // add solution previous solver
    setLastSolutionsToSystem(nlsys,scatter_uK,scatter_vvK,_u,_v,_a,_vv);
    // save last step as previous state in the system
    
    // get unknown
    nlsys->getFromSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allField);
    computeKinematicVariables(GPs,allField);
    nextStep(GPs,nlsys,true);
    //alllocate system
    nlsys->allocateSystem(_activeDofMapLocal,_ISactiveDof, 
                    activeNumberDispDofs, activeNumberStochasticDofs, activeNumberExtraDofs);
                    
    //compute mass matrix if necessary
    computeLumpedMassVector(GPs,*nlsys);

    // Time specifications
    double initTime = _t0;
    double totalTime = _tf;
    int nStep = 0;
    long long int numStepsT = 1000; // Default value
    double dt = 0;
    double timeRun = 0;
    timeRun = _t0; // Local
    dt = this->calculateCriticalTimeStep(nod, elem);
    dt = dt * _scaleFactor;
#ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

    numStepsT = long((totalTime - initTime) / dt);
    int times = _outputs->getNumFiles(); // Number of outputs in time
    long int sampler = long(numStepsT / times);
    int sampler2 = 0;

    if (timeRun < 1E-12) 
    {
        getSolutionsFormSystem(nlsys, _numNodes*numMechDofPerNode, numExtraDof,uK, vK, aK, vvK);
        if (extraDofFlag) {
            outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, vvK,
                              extraDofAdditionalOutput,
                              sim, _numDof, GPs, _outputs, dt, timeRun, _comRank);
        } else {
            outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, sim, _numDof, GPs,
                              _outputs, dt, timeRun, _comRank);
        }
    }


    if (_comRank == 0) {
        INFO("Solver used is EXPLICIT");
        INFO("The initial time step is %g", dt);
    }


#ifdef PARALLEL
    if (_comRank == 0)
#endif
    {
        INFO("done preprocessing in solver (Wall %gs, CPU %gs)", getTimeOfDay()-startWall, getCPUTime() - startTime);
    }
    
#ifdef FSI
    bool FSIstepDone =true;
    double timeRunFSI = timeRun; 
    double dtFSISubstepping = 0.;
    double dtPrev = dt;
    double dtFSI = 0.;
    if (_flagFSI)
    {
      _adapter = new FSIAdapter();
      dtFSI = _adapter->initialize(NeumannBCs);
    }
#endif   

    ////////////////////////////// TIME LOOP //////////////////////////////
    while (timeRun < totalTime - 1E-12) {
#ifdef FSI
        if (_flagFSI && FSIstepDone)
        {
            _adapter->readBlockVectorData();
            timeRunFSI += dtFSI;
            dtFSISubstepping = 0.;
            dt = dtPrev;
        }
#endif 
        if (timeRun + dt > totalTime) {
            dt = totalTime - timeRun;
            timeRun = totalTime;
        } else {
            timeRun += dt;
        }
        
#ifdef FSI
        if (_flagFSI)
        {
          if (dtFSISubstepping+dt > dtFSI - 1E-12)
          {
              FSIstepDone = true;
              timeRun -= dt;
              dtPrev = dt;
              dt = dtFSI - dtFSISubstepping;
              dtFSISubstepping = dtFSI;
              if (timeRun + dt > totalTime) {
                dt = totalTime - timeRun;
                timeRun = totalTime;
              } else {
                  timeRun += dt;
              }
          }
          else
          {
            dtFSISubstepping += dt;
            FSIstepDone = false;
          }
        }
#endif     
        //if flagActivation is true, the function to update the status of the GPs (active/inactive) is called
        if (_flagActivation) {
            updateActivation(GPs, timeRun, _activeDofMapLocal, _nActiveDof, _ISactiveDof, nNodGlobal*numMechDofPerNode, numExtraDof,
                             nodLocal, nNodGlobal, nNodLocal, _comSize, _comRank);
            //This function update the status of the GPs (active/inactive) and the number of active DOFs

        }
        sampler2++;
        //if(sampler2>=sampler){
        dt = this->calculateCriticalTimeStep(nod, elem, uK);//dt is calculated every timestep
        dt = dt * _scaleFactor;
        //}
        if (dt < 1E-20) {
            ERROR("Negative jacobian");
            exit(-1);
        }
#ifdef PARALLEL
        MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

        nStep++;
        
        nlsys->dynamicPrediction(timeRun,dt,scatter_uK);
        // apply displacement BC
        applyDisplacementBC(nlsys, _numNodes*numMechDofPerNode, timeRun, dt,nodLocal, nNodGlobal,nNodLocal);
        // get unknown
        nlsys->getFromSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allField);
        nlsys->getFromAccelerationSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allAcc);
        // compute constitutive
        computeState(GPs,allField,dt,timeRun,false);
        // compute resudial
        double stepRatio = 1;
#ifdef FSI
        stepRatio = dtFSISubstepping/dtFSI;
#endif 
        computeRightHandSide(nlsys, inerForce, intForce,extForce,totalNumberDofsPerNode,allField, allAcc, 
                            GPs,  nod, elem, NeumannBCs, _ndim, timeRun, dt,stepRatio);
        nlsys->systemSolve(timeRun,dt, _petsc_solver, _petsc_precon);
        getSolutionsFormSystem(nlsys, _numNodes*numMechDofPerNode, numExtraDof,uK, vK, aK, vvK);
#ifdef FSI
        if (_flagFSI && FSIstepDone)
        {
            _adapter->updateDisplacementInterface(dtFSI,uK,_numNodes);
            dtFSI = _adapter->advance(dtFSI);
        }
#endif  
        printForcesOut(intForce, totalDofs, dt, inerForce, 1, _forceOutputs, _t0, _tf, timeRun, _ndim, _comRank,
                       _numDof, nNodGlobal, nNodLocal, nodLocal, _comSize);
        if (!std::isfinite(uK[1]) || (extraDofFlag && !std::isfinite(vvK[1]))) { // A criteria more accurate is welcome
            ERROR("There are NANs in the displacements. Try a smaller scale factor");
            exit(-1);
        } else {
            if ((sampler2 >= sampler) || (timeRun == totalTime)) {
                if (extraDofFlag) {
                    outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, vvK,
                                      extraDofAdditionalOutput, sim, _numDof, GPs, _outputs, dt, timeRun, _comRank);
                } else {
                    outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, sim,
                                      _numDof, GPs, _outputs, dt, timeRun, _comRank);
                }
                sampler2 = 0;
                if (_comRank == 0) {
                    INFO("Time %g out of %g", timeRun, totalTime);
                }
            }
            nextStep(GPs,nlsys,false);
        }
        classExtractData::extractField(GPs,uK, vK,vvK, intForce, extForce, _numNodes,timeRun);

#ifdef FSI        
        if (_flagFSI && FSIstepDone)
        {
            if (!_adapter->isCouplingOngoing())
            {
              INFO("FSI coupling is done!");
              break;
            }
        }
#endif

    }  // end of the time loop

    _u = uK;
    _v = vK;
    _a = aK;
    _vv = vvK;
    _extF = extForce;
    
#ifdef FSI
    if (_flagFSI)
    {
        _adapter->finalize();
        delete _adapter;
    }
#endif

#ifdef PARALLEL
    if (_comRank == 0)
#endif
    {
        INFO("done solving (Wall %gs, CPU %gs)", getTimeOfDay() - startWall, getCPUTime() - startTime);
    }
    delete nlsys;
}

void classSolverExplicit::computeInternalForce(const vector<classGPs *>& GPs, vector<double>& fint, double timeRun, double dt) const
{    
    vector<int> vDofs;
    mVector vForce;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        //
        GaussPoint->getConstitutiveManager().computeInternalForceExplicit(_numNodes,GaussPoint,vDofs,vForce,timeRun,dt);
        int numDofs = vDofs.size();
        for (int j=0; j< numDofs; j++)
        {
            fint[vDofs[j]] += vForce(j);
        }
    }
    
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        // Definition of vectors
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        const vector<double>& phi = GaussPoint->getPhi();
        const vector<double>& dphi = GaussPoint->getDphi();
        const vector<int>& neighbours = GaussPoint->getNeighbours(); // GaussPointoral list of nodes in which a given GPs lies
        vector < constitutiveModels * >& GaussPointCons = GaussPoint->getConstitutiveManager().getConsModel();
        
        const vector < double >& FCurr = GaussPoint->getCurrentState().getDeformationGradient();
        const vector < double >& FPrev = GaussPoint->getPreviousState().getDeformationGradient();
        double rho = GaussPointCons[0]->getRho();
        vector<double > Piola_vis(_ndim * _ndim, 0.);
        bool withVisco = false;
        if (_flagAV == 1) { //bulk viscosity
            bulkViscosity(_elementSize, FCurr, FPrev, Piola_vis, rho, dt);
            withVisco = true;
        } else if (_flagAV == 2) { //dev viscosity
            artificialViscosity(_elementSize, FCurr, FPrev, Piola_vis, rho, dt);
            withVisco = true;
        }
        if (withVisco)
        {
            double wJ = GaussPoint->getWeightJ();
            // Assembly
            int numNeighbours = neighbours.size();
            for (int a = 0; a < numNeighbours; a++) 
            {
                for (int i = 0; i < _ndim; i++) 
                {
                    int indi11 = neighbours[a] + _numNodes * i;
                    for (int J = 0; J < _ndim; J++) 
                    {
                        fint[indi11] += wJ * Piola_vis[i * _ndim + J] * dphi[a * _ndim + J];
                    }
                }
            }
        }
    }
   
};

/*! \brief This function compute the lumped mass vector
  @param[in] GPs Array with all GPs in the domain
  @param[in] KK stiffness matrix
*/
void classSolverExplicit::computeLumpedMassVector(const vector<classGPs *>& GPs, nonlinearSystem& nlsys) const
{
    if (!nlsys.isAllocated()) 
    {
        ERROR("the nonlinear system must be allocated");
        exit(-1);
    }
    vector<int> vDofs;
    mVector mass;
    nlsys.zeroLumpedMassVector();
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
        {
            exit(-1);
            continue;
            
        }
        GaussPoint->getConstitutiveManager().computeLumpedMassVector(_numNodes,GaussPoint,vDofs,mass);
        int numDofs = vDofs.size();
        for (int j=0; j< numDofs; j++)
        {
            int row = _activeDofMapLocal[vDofs[j]];
            if (row > -1)
            {
                nlsys.addToLumpedMassVector(row,mass(j));
            };
        };
    };
}

//-----------------------------------Deviatoric Artificial viscosity with the input of FDot-------------------------------------//
/*! \brief This function computes an additional Piola-Kirchhoff tensor containing bulk viscosity 
  @param[in] elementSize Minimum distance inside the element
  @param[in] F1 Current transformation gradient
  @param[in] F0 Previous transformation gradient
  @param[out] Piola_vis Additional Piola-Kirchhoff tensor containing deviatoric viscosity
  @param[in] dt Time step
*/
void classSolverExplicit::artificialViscosity(double elementSize, const vector<double> &F1, const vector<double> &F0,
                                              vector<double> &Piola_vis, double rho, double dt) const{
    double J1 = determinantTensor(F1, _ndim);
    double J0 = determinantTensor(F0, _ndim);
    double elementSize1 = pow(J1, 1. / .3) * elementSize;

    double du = elementSize1 * (log10(J1) - log10(J0)) / dt; // dun+1 = hn+1*(logJn+1 - logJn)/dt
    double deta = 0;
    double Rho_1 = rho / J1;

    if (du < 0) {
        deta = std::max(0., -3. / 4. * Rho_1 * elementSize1 * (_c1 * du - _cl * _c0));
    } else {
        deta = 0;
    }

    vector < double > Cauchy_vis(_ndim * _ndim,
                                 0.);
    vector < double > Fdot(_ndim * _ndim,
                           0.);
    sumTensorTensor(F1, F0, Fdot, _ndim, -1);
    multiSTensorScalar(Fdot, _ndim, _ndim, 1 / dt);
    vector < double > inv_F(_ndim * _ndim,
                            0.);
    inverse(F1, _ndim, inv_F);

    vector < double > multi(_ndim * _ndim,
                            0.);
    multTensorTensor3(Fdot, _ndim, _ndim, inv_F, _ndim, multi);
    vector<double> sym_multi = symSTensor3(multi, _ndim);
    vector<double> dev_multi = devSTensor3(sym_multi, _ndim);
    multiSTensorScalar(dev_multi, _ndim, _ndim, 2. * deta);
    multSTensor3SecondTranspose(Cauchy_vis, _ndim, _ndim, inv_F, _ndim, Piola_vis);
    multiSTensorScalar(Piola_vis, _ndim, _ndim, J1);    //Cauchy->PK1
}

//--------------------------------Bulk viscosity ----------------------------------------------------//
/*! \brief This function computes an additional First Piola-Kirchhoff tensor containing bulk viscosity 
  @param[in] elementSize Minimum distance inside the element
  @param[in] F1 Current transformation gradient
  @param[in] F0 Previous transformation gradient
  @param[out] Piola_vis Additional Piola-Kirchhoff tensor containing bulk viscosity
  @param[in] dt Time step
*/
void
classSolverExplicit::bulkViscosity(double elementSize, const vector<double>& F1,const vector<double>& F0, vector<double> &Piola_vis,
                                   double rho, double dt) const {
    vector < double > Fdot(_ndim * _ndim,
                           0);
    if (dt == 0) {
        fill(Fdot.begin(), Fdot.end(), 0.);
    } else {
        sumTensorTensor(F1, F0, Fdot, _ndim, -1);
        multiSTensorScalar(Fdot, _ndim, _ndim, 1 / dt);
    }
    double trD = getTrace(Fdot, _ndim);

    double J1 = determinantTensor(F1, _ndim);
    double Rho_1 = rho / J1;

    double qBV = Rho_1 * elementSize * (_c1 * _c1 * elementSize * abs(trD) * std::min(0., trD) +
                                        _cl * _c0 * trD); //c1 qrd coef, cl linear coef
    vector < double > Cauchy_vis(_ndim * _ndim,
                                 0.);
    for (int i = 0; i < _ndim; i++) {
        Cauchy_vis[i + i * _ndim] = -1 * qBV;
    }

    vector < double > inv_F(_ndim * _ndim,
                            0.);
    inverse(F1, _ndim, inv_F);
    multSTensor3SecondTranspose(Cauchy_vis, _ndim, _ndim, inv_F, _ndim, Piola_vis);
    multiSTensorScalar(Piola_vis, _ndim, _ndim, -1 * J1);    //Cauchy->PK1
}


void classSolverExplicit::applyDisplacementBC(nonlinearSystem* nlsys, int numMechDofs, double timeRun, double dt, 
                        const vector<int> &nodLocal, int nNodGlobal, int nNodLocal)
{
    for (int i = 0; i < _DNod.size(); i++) 
    {
        double U_diri=0;
        double V_diri=0;
        int ind = _DNod[i]->getMe();
        V_diri = _DNod[i]->getUU();
        int flagDirichlet = _DNod[i]->getflagDiri();
        if (flagDirichlet == 1) { //if flagDirichlet=1, Dirichlet BC are applied instantaneously and are applied just one time
            U_diri = V_diri;
            _DNod[i]->resetflagDiri();
            V_diri = V_diri / dt;
        } else if (flagDirichlet == -1) { //if flagDirichlet=-1, Dirichlet BC are applied as ramp from t0 to tf
            U_diri = V_diri * dt;
        } else {
            V_diri = 0;
        }
        int dofGlobal = ind;
#ifdef PARALLEL
        locDofToGlobDof(nNodGlobal, nNodLocal, ind, dofGlobal, nodLocal, false);
#endif
        nlsys->addToSolution(dofGlobal,U_diri);
        if (ind < numMechDofs)
        {
            nlsys->setToVelocitySolution(dofGlobal,U_diri);
            nlsys->setToAccelerationSolution(dofGlobal,0.);
        };
    };
};

void classSolverExplicit::setLastSolutionsToSystem(nonlinearSystem* nlsys, 
                                const vector<int>& scatter_uK,
                                const vector<int>& scatter_vvK,
                                const vector<double>& uK,
                                const vector<double>& vK,
                                const vector<double>& aK,
                                const vector<double>& vvK) const
{
    if (scatter_uK.size() > 0)
    {
        nlsys->setToSolution(scatter_uK,uK);
        nlsys->setToVelocitySolution(scatter_uK,vK);
        nlsys->setToAccelerationSolution(scatter_uK,aK);
    }
    
    if (scatter_vvK.size() >0)
    {
        nlsys->setToSolution(scatter_vvK,vvK);
    }
}

void classSolverExplicit::computeRightHandSide(nonlinearSystem* nlsys, vector<double>&inerForce, vector<double>&intForce, vector<double>&extForce,
                            int totalNumberDofsPerNode, const vector<double> &allField, const vector<double>& allAcc, 
                            vector<classGPs *> &GPs, 
                            vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, double timeRun,
                            double dt, double FSIstepRatio)
{
    
    setAll(intForce,0);
    setAll(extForce,0);
    setAll(inerForce,0);
    computeInternalForce(GPs,intForce,timeRun,dt);
    map<pair<int,int>,double> mapTmp;
    NeumannBCManagement(allField, GPs, nod, elem, NeumannBCs, ndim, extForce, timeRun, dt, 
                        false, mapTmp);
    //computeInertialForce(GPs, allAcc,inerForce);
        // 
#ifdef FSI
    if (_flagFSI)
    {
        _adapter->updateForceInterface(extForce,_numNodes, FSIstepRatio);
    }
#endif      
    nlsys->zeroRHS();
    MPI_Barrier(MPI_COMM_WORLD);
    for (int row=0; row< _numNodes*totalNumberDofsPerNode; row++)
    {
        double val = -intForce[row] + extForce[row];
        nlsys->addToRHS(row,val);
    }
};

void classSolverExplicit::getSolutionsFormSystem(nonlinearSystem* nlsys, 
                                int numMechDofs,  int numExtraDof,
                                vector<double>& uK,
                                vector<double>& vK,
                                vector<double>& aK,
                                vector<double>& vvK) const
                                
{
   
    if (numMechDofs > 0)
    {
        nlsys->getFromSolutionByRange(0,numMechDofs,uK);
        nlsys->getFromVelocitySolutionByRange(0,numMechDofs,vK);
        nlsys->getFromAccelerationSolutionByRange(0,numMechDofs,aK);
    }
    if (numExtraDof >0)
    {
        nlsys->getFromSolutionByRange(numMechDofs,numMechDofs+numExtraDof,vvK);
    }
    
}