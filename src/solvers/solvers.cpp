
//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file solvers.cpp
  \brief This file contains all common functions of the solvers. See MuPhiSim documentation related to the temporal discretization
*/
#include "solvers.h"
#include <petscsnes.h>
#include <petscksp.h>
#include "classNodes.h"
#include "classElements.h"
#include "classGPs.h"
#include "dataExtraction.h"
#include <sys/resource.h>
#include <chrono>

double getCPUTime() {
    static struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return (double) r.ru_utime.tv_sec + 1.e-6 * (double) r.ru_utime.tv_usec;
};

double getTimeOfDay()
{
    auto t = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t.time_since_epoch()).count();
}



/*! \brief Constructor of the class
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of degrees of fredom
  @param[in] DNod Dirichlet nodes
  @param[in] scaleFactor Value to be multiplied by the critical time step (<1 for Explicit, >= 1 Implicit)
  @param[in] outputs Object with all data needed for the outputs
  @param[in] forceOutputs Object with all data needed for the force outputs
  @param[in] t0 Initial time for this solver
  @param[in] tf Final for this solver
  @param[in] solverType String to identify the type of solver
  @param[in] activeDofMapLocal Array for parallel simulations with the local active Dof of a processor
  @param[in] nActiveDof Global number of active Dof in parallel simulations (zero in sequential simulations)
  @param[in] ISactiveDof Array for parallel simulations with the mapping for the global displacements vector (empty in sequential simulations)
*/

solvers::solvers(vector<classNodes *> &nod, int ndim, int numDof, vector<classDirichlet *> &DNod, double scaleFactor,
                 classOutputs *&outputs, double t0, double tf, string solverType, const vector<int> &activeDofMapLocal,
                 const vector<int>& nActiveDof, const vector<int> &ISactiveDof, classPrintForces *&forceOutputs,
                 bool flagActivation, string &petsc_solver, string &petsc_precon, bool flagFSI):
                 _timeSteppingPlan(1e-12), _comRank(0), _comSize(1), _maxFailedSteps(50),
                 _tangentByPerturbation(false), _perturbationTol(0), _numIterationsIncreaseTimeStep(5),
                 _maxTimeStep(1e10),
                 _adaptiveTimeStep(NULL), _trueDispField(NULL)
                 {
    _numNodes = nod.size();
    _numDof = numDof;
    _ndim = ndim;
    _scaleFactor = scaleFactor;
    _numSteps = 0;
    _maxIterations = 50;
    _absTol = 1e-12;
    _relTol = 1e-6;
    _outputs = outputs; // Default outputs. U, VMS, Cauchy
    _forceOutputs = forceOutputs;
    _t0 = t0;
    _tf = tf;
    _DNod = DNod;
    _solverType = solverType;
    // Parallel arrays. Empty
    _activeDofMapLocal = activeDofMapLocal;
    _activeDofMapLocalInitial = activeDofMapLocal;
    _nActiveDof = nActiveDof;
    _ISactiveDof = ISactiveDof;
    _flagActivation = flagActivation;
    _petsc_solver = petsc_solver;
    _petsc_precon = petsc_precon;
    _flagFSI = flagFSI;
    _stiffnessMatrixType = solvers::Current;
    _stiffnessMatrixFixedIteration = -1; // applicable only _stiffnessMatrixType = Step
#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &_comSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &_comRank);
#endif    

#ifdef FSI
    if (flagFSI)
    {
      _adapter = new FSIAdapter();
    };
#endif 

};

solvers::~solvers()
{
    if (_adaptiveTimeStep != NULL)
    {
        delete _adaptiveTimeStep;
        _adaptiveTimeStep=NULL;
    }
    if (_trueDispField !=NULL)
    {
        delete _trueDispField;
        _trueDispField = NULL;
    }
}

/*! \brief This function is used to set the tangent by perturbation in implicit scheme
*/
void solvers::setTangentByPerturbation(bool fl, double tol)
{
    _tangentByPerturbation = fl;
    _perturbationTol = tol;
    if (_tangentByPerturbation and _comRank == 0)
    {
        INFO("Tangent is computed by perturbation with perturbed tolerance = %g",tol);
    }
}

/*! \brief This function is used to set number of allowable failed steps
*/
void solvers::setAllowableFailedSteps(int nstep)
{
  if (_comRank == 0)
  {
    INFO("number of allowable failed steps: %d",nstep);
  }
  _maxFailedSteps = nstep;
};

/*! \brief This function is used to set number of steps
*/
void solvers::setNumberOfSteps(int nstep)
{
  if (_comRank == 0)
  {
    INFO("number of steps: %d",nstep);
  }
  _numSteps = nstep;
};

/*! \brief This function is used to set absolute tolerance
*/
void solvers::setAbsoluteTolerance(double absTol)
{
    if (_comRank == 0)
    {
        INFO("absolute tolerance: %.16g",absTol);
    }
    _absTol = absTol;
}

/*! \brief This function is used to set relative tolerance
*/
void solvers::setRelativeTolerance(double relTol)
{
    if (_comRank == 0)
    {
        INFO("relative tolerance: %.16g",relTol);
    }
    _relTol = relTol;
}

/*! \brief This function is used to set maximal number of iteration in each time step
*/
void solvers::setMaximalNumberOfIterations(int maxIt)
{
    if (_comRank == 0)
    {
        INFO("max number of iterations: %d",maxIt);
    }
    _maxIterations = maxIt;
}

/*! \brief This function is used to set maximal time step
*/
void solvers::setMaximalTimeStep(double dtMax)
{
    if (_comRank == 0)
    {
        INFO("max time step: %g",dtMax);
    }
    _maxTimeStep = dtMax;
}

/*! \brief This function is used to set time stepping plan
*/
void solvers::setTimeSteppingPlan(const TimeStepping& plan)
{
    _timeSteppingPlan  = plan;
    if (_comRank == 0)
    {
        INFO("Time stepping plan is used");
    }
};

/*! \brief This function is used to set the number of iterations allowng to increase time step after reduction
*/
void solvers::setNumIterationsIncreaseTimeStep(int nIter)
{
    _numIterationsIncreaseTimeStep = nIter;
};

/*! \brief allocate true displacement field for error analysis
*/
void solvers::allocateTrueDisplacementFieldForErrorAnalysis(const char what[], const std::vector<double>& data)
{
    if (_trueDispField != NULL)
    {
        delete _trueDispField;
    } 
    _trueDispField = trueDisplacementField::allocate(what, data);
}

/*! \brief This function is used to set the adaptive timestep object
*/
void solvers::setAdaptiveTimeStep(const char what[], const vector<double>& data)
{
    if (_adaptiveTimeStep == NULL)
    {
        delete _adaptiveTimeStep;
    }
    _adaptiveTimeStep = AdaptiveTimeStep::createAdaptiveTimeStep(what,data);
};


/*! \brief This function is used to set stiffness type in iterative solver
*/
void solvers::setStiffnessMatrixType(solvers::StiffnessMatrixType type, int iteration)
{
    _stiffnessMatrixType = type;
    _stiffnessMatrixFixedIteration = iteration;
    
    if (_comRank == 0)
    {
        if (type == solvers::Current)
        {
            INFO("Stiffness matrix is updated every iteration in each time step");
        }
        else if (type == solvers::Initial)
        {
            INFO("Stiffness matrix is only updated at the beginning of each time step");
        }
        else if (type == solvers::Iteration)
        {
            INFO("Stiffness matrix is not updated after %dth iteration of each time step",_stiffnessMatrixFixedIteration);
        }
        else if (type == solvers::Interval)
        {
            INFO("Stiffness matrix is updated after %d iterations of each time step",_stiffnessMatrixFixedIteration);
        }
        
        else if (type == solvers::Beginning)
        {
            INFO("Stiffness matrix is re-used from the beginning of simulation");
        }
        else
        {
            ERROR("This option for Stiffness matrix is not implemented");
        }
    }
}


/*! \brief This function calculates the critical time step based on the sound speed in the material and the mesh
  @param[in] nod Array with all nodes in the domain
  @param[in] elem Array with all elements (or background cells) in the domain
  return Critical dt
*/
double solvers::calculateCriticalTimeStep(vector<classNodes *> &nod, vector<classElements *> &elem) {
    vector<double> allTimeSteps(elem.size(),0.);
    for (int i=0; i< elem.size(); i++)
    {
        double minDist = elem[i]->getCharLength();
        vector<constitutiveModels *>&tempCons = elem[i]->getConsModel();
        double soundSpeed  = tempCons[0]->soundSpeed();
        allTimeSteps[i] = minDist / soundSpeed;
    }
     
    return *(std::min_element(allTimeSteps.begin(),allTimeSteps.end()));
};

/*! \brief This function calculates the critical time step based on the sound speed in the material and the mesh
  @param[in] nod Array with all nodes in the domain
  @param[in] elem Array with all elements (or background cells) in the domain
  @param[in] uK Displacement of all nodes

  return Critical dt
*/
double solvers::calculateCriticalTimeStep(vector<classNodes *> &nod, vector<classElements *> &elem,
                                        const vector<double> &uK) {
    vector<double> allTimeSteps(elem.size(),0.);
    for (int i=0; i< elem.size(); i++)
    {
        double minDist = elem[i]->getCharLength();
        const vector<int>& eleNodes = elem[i]->getMyNodes();
        int numNodes = eleNodes.size();
        vector<vector<double> > xyzDeformed(numNodes, vector<double>(_ndim,0));
        for (int j = 0; j < numNodes; j++) 
        {
            const vector<double>& initialCoods = nod[eleNodes[j]]->getXYZ();//get undeformed position
            xyzDeformed[j] = initialCoods;
            for (int k= 0; k < _ndim; k++) 
            {
                xyzDeformed[j][k] += uK[_numNodes * k + eleNodes[j]];
            }
        }
        double dist = elem[i]->characteristicLength(_ndim, xyzDeformed);//gets the inscribed n-sphere for the deformed configuration
        vector<constitutiveModels *>&tempCons = elem[i]->getConsModel();
        double soundSpeed  = tempCons[0]->soundSpeed();
        allTimeSteps[i] = std::min(minDist,dist)/soundSpeed;
    }
    return *(std::min_element(allTimeSteps.begin(),allTimeSteps.end()));
};


/*! \brief This function updates the GPs and the system for the next time step
  @param[inout] GPs Array with all GPs in the domain
  @param[inout] lsys nonlinear system
  @param[in] resetInitialGP true if the initial state is reset with current state
*/
void solvers::nextStep(vector<classGPs *> &GPs, nonlinearSystem* lsys,  bool resetInitialGP) const
{
    int numGPs = GPs.size();
    for (int i=0; i < numGPs; i++)
    {
        GPs[i]->nextStep();
    }
    lsys->nextStep();
    if (resetInitialGP)
    {
        for (int i=0; i < numGPs; i++)
        {
            GPs[i]->getInitialState() = (GPs[i]->getCurrentState());
        }
    }
};

/*! \brief This function updates the GPs and the system for the next time step
  @param[inout] GPs Array with all GPs in the domain
  @param[inout] lsys nonlinear system
*/
void solvers::resetToPreviousStep(vector<classGPs *> &GPs, nonlinearSystem* lsys) const
{
    int numGPs = GPs.size();
    for (int i=0; i < numGPs; i++)
    {
        GPs[i]->resetToPreviousState();
    }
    lsys->resetToPreviousSolution();
};

/*! \brief This function makes a check point of the GPs and the system
  @param[inout] GPs Array with all GPs in the domain
  @param[inout] lsys nonlinear system
*/
void solvers::makeACheckPoint(vector<classGPs *> &GPs, nonlinearSystem* lsys) const
{
    int numGPs = GPs.size();
    for (int i=0; i < numGPs; i++)
    {
        GPs[i]->makeACheckPoint();
    }
    lsys->makeACheckPoint();
}

/*! \brief This function restores a check point of the GPs and the system
  @param[inout] GPs Array with all GPs in the domain
  @param[inout] lsys nonlinear system
*/
void solvers::restoreACheckPoint(vector<classGPs *> &GPs, nonlinearSystem* lsys) const
{
    int numGPs = GPs.size();
    for (int i=0; i < numGPs; i++)
    {
        GPs[i]->restoreACheckPoint();
    }
    lsys->restoreACheckPoint();
}
    


/*! \brief This function updates the status (active/inactive) of the GPs and active DOFs vectors and number of active DOFs
  @param[inout] GPs Array with all GPs in the domain
  @param[in] timeRun Global time of the simulation
  @param[inout] activeDofMapLocal Array for parallel simulations with the local active Dof of a processor (empty in sequential simulations)
  @param[inout] nActiveDof Global number of active Dof
  @param[inout] ISactiveDof Array for parallel simulations with the mapping for the global displacements vector (empty in sequential simulations)
  @param[in] totalDofs Global number of mechanical DOF
  @param[in] numExtraDof Global number of extra Dof
  @param[in] nodLocal Index mapping from local to global node in parallel simulations (empty for sequential)
  @param[in] nNodGlobal Number of nodes of the global-whole model in parallel simulations (zero for sequential)
  @param[in] nNodLocal Local number of nodes in parallel simulations (zero for sequential)
  @param[in] nprcs int number of processors working on the simulation
  @param[in] rank int indicating which processor is calling the function
*/
void solvers::updateActivation(vector<classGPs *> &GPs, double timeRun, vector<int> &_activeDofMapLocal,
                               vector<int> &_nActiveDof, vector<int> &_ISactiveDof, int totalDofs, int numExtraDof,
                               const vector<int>& nodLocal, int nNodGlobal, int nNodLocal, int nprcs, int rank) {

    fill(_activeDofMapLocal.begin(), _activeDofMapLocal.end(), -2);
    //This loop updates the GPS status and _activeDofMapLocal vector
    int TotalGlobaDof =0;
    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) {
        classGPs *GaussPoint = *ite;
        const vector<int>& neighbours = GaussPoint->getNeighbours();
        int TotalConsDof = GaussPoint->getConstitutiveManager().getTotalNumDofPerNode();
#ifdef PARALLEL
        TotalGlobaDof=TotalConsDof*nNodGlobal; //Global numb of DOF (mechanical+extra) in parallel simulations
#endif

        //GPs status update (active/inactive)
        GaussPoint->getConstitutiveManager().checkActivation(GaussPoint,timeRun);
        if (GaussPoint->getActivate() == 1)
        {
            //Vector _activeDofMapLocal update. Globally for sequential simulations and locally for parallel simulations
            int sizeNeighbours = neighbours.size();
            for (int a = 0; a < sizeNeighbours; a++) {
                for (int ii = 0; ii < TotalConsDof; ii++) {
                    int indi1 = neighbours[a] + _numNodes * ii; 
                    _activeDofMapLocal[indi1] = _activeDofMapLocalInitial[indi1];//if the GP is active, we reactive the node again
                }
            }
        }
    }


#ifdef PARALLEL
    // The following code shares with the other processors the information related to _activeDofMapLocal of each processor (parallel simulations)
    int dofGlobal;
    int localDof;
    localDof = _numDof + numExtraDof;
    vector<int> tempV(localDof); //temporary vector with mapping to global
    for (int j = 0; j < localDof; j++) {
        locDofToGlobDof(nNodGlobal, nNodLocal, j, dofGlobal, nodLocal, false);
        tempV[j] = dofGlobal;
    }
    // Get total length
    vector<int> lenPerProc;//global array storing length of activeDofMapLocal of every process
    if (rank == 0) {
        lenPerProc.resize(nprcs);
    }
    int localSize = _activeDofMapLocal.size();
    MPI_Gather(&localSize, 1, MPI_INT, &lenPerProc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

    int globalSize = 0;
    if (rank == 0) {
        for (int i = 0; i < nprcs; i++) {
            globalSize += lenPerProc[i];
        }
    }
    //combine all the local indices and activeDofMapLocal into one array
    vector<int> globalIndexing;
    vector<int> activeMapArray(_activeDofMapLocal.size());
    vector<int> activeMapGlobalUnOrd;
    vector<int> locToGlobLoc(tempV.size());

    for (int i = 0; i < localDof; i++) {
        activeMapArray[i] = _activeDofMapLocal[i];//copy vector onto array
        locToGlobLoc[i] = tempV[i];
    }

    if (rank == 0) {
        activeMapGlobalUnOrd.resize(globalSize); 
        globalIndexing.resize(globalSize);
    }

    vector<int> displc;
    if (rank == 0) {
        displc.resize(nprcs);
        displc[0] = 0;
        for (int i = 1; i < nprcs; i++)
            displc[i] = displc[i - 1] + lenPerProc[i - 1];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(&activeMapArray[0], _activeDofMapLocal.size(), MPI_INT, &activeMapGlobalUnOrd[0], &lenPerProc[0], &displc[0],
               MPI_INT, 0, MPI_COMM_WORLD);//in rank order
    MPI_Gatherv(&locToGlobLoc[0], tempV.size(), MPI_INT, &globalIndexing[0], &lenPerProc[0], &displc[0], MPI_INT, 0, MPI_COMM_WORLD);
    vector<int> activeGlobal;
    if (rank == 0) {
        activeGlobal.resize(TotalGlobaDof, -3); // up to this point, a node can be active (>-1) in one processor and inactive (-2) in another, we initialize the vector to -3
        for (int j = 0; j < globalSize; j++) {
            activeGlobal[globalIndexing[j]] = max(activeMapGlobalUnOrd[j], activeGlobal[globalIndexing[j]]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

#endif
    //_nActiveDof and _ISactiveDof (just for parallel simulations) update
#ifdef PARALLEL
    if (rank == 0) {
        int count = 0;
        for (int b = 0; b < activeGlobal.size(); b++) { //can we use TotalGlobalDof?
            if (activeGlobal[b] > -1) {
                activeGlobal[b] = count;
                count++;
            }
        }
    }

   if (rank == 0) {
        MPI_Send(&activeGlobal[0], TotalGlobaDof, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else {
        activeGlobal.resize(TotalGlobaDof, 0);
        MPI_Recv(&activeGlobal[0], TotalGlobaDof, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank < (nprcs - 1)) {
            MPI_Send(&activeGlobal[0], TotalGlobaDof, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

      for (int i = 0; i < tempV.size(); i++) {
          _activeDofMapLocal[i] = activeGlobal[tempV[i]];
      }

    _nActiveDof = {0, 0}; //how to change that
    for (int b = 0; b < activeGlobal.size(); b++) {
        if (activeGlobal[b] > -1) {
            if (b < totalDofs) { //here we need the total num of DOF mechanical
                _nActiveDof[0] += 1;
            } else {
                _nActiveDof[1] += 1;
            }
        }
    }
    
    int totalActive = _nActiveDof[0] + _nActiveDof[1];

    _ISactiveDof.resize(totalActive, 0); //check
    int k = 0;
    for (int j = 0; j < activeGlobal.size(); j++) {
        if (activeGlobal[j] > -1) {
            _ISactiveDof[k] = j;
            k++;
        }
    }

    activeGlobal.clear();
#else
    _nActiveDof = {0, 0}; //how to change that
    int count = 0;
    for (int b = 0; b < _activeDofMapLocal.size(); b++) {
        if (_activeDofMapLocal[b] > -1) {
            _activeDofMapLocal[b] = count;
            count++;
            if (b < _numDof) { //here we need the total num of DOF mechanical
                _nActiveDof[0] += 1;
            } else {
                _nActiveDof[1] += 1;
            }
        }
    }
#endif

}

/*! \brief This function preallocate stiffness matrix
  @param[inout] GPs Array with all GPs in the domain
  @param[in] KK_total stiffness matrix
*/
void solvers::preAllocateMatrix(const vector<classGPs *>& GPs, classMatrix &KK_total) const
{
    if (KK_total.isPreallocated()) return;
    vector<int> vDofs;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        
        GaussPoint->getLocalDofs(_numNodes,vDofs);
        int vDofsSize = vDofs.size();
        for (int ii = 0; ii < vDofsSize; ii++) 
        {
            int row = _activeDofMapLocal[vDofs[ii]];
            for (int jj = 0; jj < vDofsSize; jj++) 
            {
                int col = _activeDofMapLocal[vDofs[jj]];
                if (row > -1 and col > -1) 
                {
                    KK_total.insertSparsityPattern(row,col);
                }
            }
        }
    }
    KK_total.preallocate();
};

/*! \brief This function updates kinematic at GPs
  @param[inout] GPs Array with all GPs in the domain
  @param[in] uKK Array with the displacements of all nodes
*/
void solvers::computeKinematicVariables(vector<classGPs *>& GPs, const vector<double> &uKK) const
{
    // compute kinematic vars
    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        GaussPoint->getConstitutiveManager().computeKinematicVariables(GaussPoint,uKK,_numNodes);
    };
}

/*! \brief This function updates constitutive at GPs
  @param[inout] GPs Array with all GPs in the domain
  @param[in] uKK Array with the displacements of all nodes
*/
void solvers::computeState(vector<classGPs *>& GPs, const vector<double> &uKK, double dt, double timeRun, bool withStiffness) const
{
    // compute kinematic vars
    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        GaussPoint->getConstitutiveManager().computeKinematicVariables(GaussPoint,uKK,_numNodes);
    };
        
    // update constitutive laws
    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        if (_tangentByPerturbation)
        {
            GaussPoint->getConstitutiveManager().constitutive(GaussPoint,dt,timeRun,false);
        }
        else
        {
            GaussPoint->getConstitutiveManager().constitutive(GaussPoint,dt,timeRun,withStiffness);
        }
    };
};

/*! \brief This function compute the internal force vector
  @param[int] GPs Array with all GPs in the domain
  @param[int] aK acceleration vector
  @param[out] fint Array of inertial force vector
*/
void solvers::computeInertialForce(const vector<classGPs *>& GPs, const vector<double> &aK, vector<double>& finer) const
{
    vector<int> vDofs;
    mVector vForce;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        //
        GaussPoint->getConstitutiveManager().computeInertialForce(_numNodes,aK,GaussPoint,vDofs,vForce);
        int numDofs = vDofs.size();
        for (int j=0; j< numDofs; j++)
        {
            finer[vDofs[j]] += vForce(j);
        }
    }
}

