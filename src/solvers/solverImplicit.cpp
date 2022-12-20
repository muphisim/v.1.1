//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file solverImplicit.cpp
  \brief This file contains all functions related to the implicit non-linear static solver
*/
#include "solverImplicit.h"
#include "output.h"
#include "boundaryConditions.h"
#include "maths.h"
#include "dataExtraction.h"
#include "nonlinearSystem.h"
#include "commandLine.h"

/*! \brief Constructor of the class
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of mechanical degrees of fredom
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
classSolverImplicitBase::classSolverImplicitBase(vector<classNodes *> &nod, int ndim, int numDof,
                                                     vector<classDirichlet *> &DNod, double scaleFactor,
                                                     classOutputs *&outputs, double t0, double tf, string solverType,
                                                     const vector<int> &activeDofMapLocal, const vector<int>& nActiveDof,
                                                     const vector<int> &ISactiveDof, classPrintForces *&forceOutputs,
                                                     bool flagActivation, string &petsc_solver, string &petsc_precon, bool flagFSI)
        : solvers(nod, ndim, numDof, DNod,
                  scaleFactor, outputs, t0, tf,
                  solverType, activeDofMapLocal,
                  nActiveDof, ISactiveDof,
                  forceOutputs, flagActivation,
                  petsc_solver, petsc_precon, flagFSI){



}

/*! \brief This function compute the stiffness matrix
  @param[in] GPs Array with all GPs in the domain
  @param[in] KK stiffness matrix
*/
void classSolverImplicitBase::computeStiffness(double timeRun, double dt, vector<double>& uK,
                                               const  vector<classGPs *>& GPs, nonlinearSystem* nlsys) const
{
    double startTimeStiffness = getCPUTime();
    double startWallStiffness = getTimeOfDay();
    
    if (nlsys->getMatrix() == NULL)
    {
        ERROR("The stiffness matrix has not been allocated yet!");
        exit(-1);
    }
    classMatrix& KK = *(nlsys->getMatrix());
    if (!KK.isPreallocated())
    {
        preAllocateMatrix(GPs,KK);
    }
    nlsys->zeroMatrix();
    //
    vector<int> vDofs;
    mMatrix stiffness;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        
        GaussPoint->getConstitutiveManager().computeStiffnessMatrix(timeRun,dt,uK,_numNodes,GaussPoint,vDofs,stiffness,
                            _tangentByPerturbation,_perturbationTol);
        int numDofs = vDofs.size();
        for (int j=0; j< numDofs; j++)
        {
            int row = _activeDofMapLocal[vDofs[j]];
            if (row > -1)
            {
                for (int k=0; k< numDofs; k++)
                {
                    int col = _activeDofMapLocal[vDofs[k]];
                    if (col > -1) 
                    {
                        KK.addIJ(row, col, stiffness(j,k));
                    };
                };
            };
        };
    };
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Let petsc copy all the values to export before I do anything else
#endif
    KK.assembly();
    if (_comRank==0)
        INFO("Time for computeStiffness (Wall %gs, CPU %gs)", getTimeOfDay() - startWallStiffness, getCPUTime() - startTimeStiffness);
};

/*! \brief This function compute the mass matrix
  @param[in] GPs Array with all GPs in the domain
  @param[in] KK stiffness matrix
*/
void classSolverImplicitBase::computeMassMatrix(const vector<classGPs *>& GPs, classMatrix& M) const
{
    if (!M.isPreallocated())
    {
        preAllocateMatrix(GPs,M);
    }
    M.zeroEntries();
    //
    vector<int> vDofs;
    mMatrix mass;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        
        GaussPoint->getConstitutiveManager().computeMassMatrix(_numNodes,GaussPoint,vDofs,mass);
        int numDofs = vDofs.size();
        for (int j=0; j< numDofs; j++)
        {
            int row = _activeDofMapLocal[vDofs[j]];
            if (row > -1)
            {
                for (int k=0; k< numDofs; k++)
                {
                    int col = _activeDofMapLocal[vDofs[k]];
                    if (col > -1) 
                    {
                        M.addIJ(row, col, mass(j,k));
                    };
                };
            };
        };
    };
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD); //Let petsc copy all the values to export before I do anything else
#endif
    M.assembly();
};
/*! \brief This function compute the internal force vector
  @param[int] GPs Array with all GPs in the domain
  @param[out] fint Array of internal force vector
*/
void classSolverImplicitBase::computeInternalForce(const vector<classGPs *>& GPs, vector<double>& fint, double timeRun, double dt) const
{
    vector<int> vDofs;
    mVector vForce;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        //
        GaussPoint->getConstitutiveManager().computeInternalForce(_numNodes,GaussPoint,vDofs,vForce);
        int numDofs = vDofs.size();
        for (int j=0; j< numDofs; j++)
        {
            fint[vDofs[j]] += vForce(j);
        }
    }
};

        

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

void classSolverImplicitBase::calculate(vector<classNeumannBCs *> &NeumannBCs, vector<classGPs *> &GPs, int sim,
                                          vector<classNodes *> &nod, vector<classElements *> &elem, const vector<int> &nodFEM,
                                          const vector<int> &nodMM, const vector<int> &nodShare, const vector<int> &elemFEM,
                                          const vector<int> &elemMM, vector<classElements *> &geoShape,
                                          vector<constitutiveModels *> &consMod, const vector<int> &nodLocal, int nNodGlobal,
                                          int nNodLocal) {

    //
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
    nonlinearSystem* nlsys = createSystem();
    // allocate solution
    
    int totalDofs = nNodGlobal*totalNumberDofsPerNode;
    nlsys->allocateSolution(totalDofs,scatter_all);
    // add solution previous solver
    setLastSolutionsToSystem(nlsys,scatter_uK,scatter_vvK,_u,_v,_a,_vv);
    // daa at GP
    nlsys->getFromSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allField);
    computeKinematicVariables(GPs,allField); // some fields need to be initialised their kinematic variablem at GP
    nlsys->getFromAccelerationSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allAcc);
    
    // save last step as previous state in the system
    nextStep(GPs,nlsys,true);
    
    //alllocate system
    nlsys->allocateSystem(_activeDofMapLocal,_ISactiveDof, 
                    activeNumberDispDofs, activeNumberStochasticDofs, activeNumberExtraDofs);
    //compute mass matrix if necessary
    if (nlsys->withMassMatrix())
    {
        computeMassMatrix(GPs,*(nlsys->getMassMatrix()));
    }
    
    _timeSteppingPlan.setMaximumTimeStep(_maxTimeStep);
    
    // Time specifications
    double totalTime = _timeSteppingPlan.getEndTime(); // For this particular solver
    double initTime = _timeSteppingPlan.getStartTime();
    if (_timeSteppingPlan.getNumSteps() <= 0)
    {
        double dt = this->calculateCriticalTimeStep(nod, elem);
#ifdef PARALLEL
        double dtMin;
        MPI_Allreduce(&dt,&dtMin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        dt = dtMin;
#endif
        dt = dt * _scaleFactor;
        int numSteps = std::max(1,int((totalTime - initTime) / dt));
        _timeSteppingPlan.setNumberOfSteps(numSteps);
    }
    _timeSteppingPlan.initFromOptions();

    // Outputs
    int times = _outputs->getNumFiles();
    int sampler = _timeSteppingPlan.getNumSteps() / times;
    double savingTimeInterval=totalTime/times;
    double savingTimeIntervalRegister=0;
    int sampler2 = 0;
    int ite2 = 0;
    bool flagConver = true;

    if (_comRank == 0) {
        INFO("Solver used is Implicit Static");
    }
    if (_comRank == 0) {
        INFO("Compressed row storage format");
    }

    if (initTime < 1E-12)
    {
        // save initial view
        getSolutionsFormSystem(nlsys, _numNodes*numMechDofPerNode, numExtraDof,uK, vK, aK, vvK);
        if (extraDofFlag) {
            outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, vvK,
                              extraDofAdditionalOutput, sim, _numDof, GPs, _outputs, 0., initTime,
                              _comRank);
        } else {
            outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, sim,
                              _numDof, GPs, _outputs, 0., initTime, _comRank);
        }
    }

#ifdef PARALLEL
    if (_comRank == 0)
#endif
    {
        INFO("done preprocessing in solver (Wall %gs, CPU %gs)", getTimeOfDay()-startWall, getCPUTime() - startTime);
    }

    /////////////////////////////////////  TIME LOOP
    double dt = _timeSteppingPlan.getTimeStep();
    double timeRun = _timeSteppingPlan.getStartTime();
    bool tangentEstimation = true;
    int nStep = 0; // index value
    int numFailedSteps = 0;
#ifdef FSI
    bool FSIstepDone =true;
    double timeRunFSI = timeRun; 
    double dtFSISubstepping = 0.;
    double dtFSI = 0.;
    if (_flagFSI)
    {
      dtFSI = _adapter->initialize(NeumannBCs);
      _timeSteppingPlan.forceSteppingScheme(true);
      dt = dtFSI; // initialize time step by preCICE
      _timeSteppingPlan.setTimeStep(dt);
    };
#endif 
    while (_timeSteppingPlan.endCheck()) 
    {
#ifdef FSI
        if (_flagFSI && FSIstepDone)
        {
            if(_adapter->isActionRequired(_adapter->actionWriteIterationCheckpoint()))
            {
                makeACheckPoint(GPs,nlsys);
                _adapter->makeACheckPoint(timeRunFSI,dtFSI);
                _adapter->markActionFulfilled(_adapter->actionWriteIterationCheckpoint());
            };
            _adapter->readBlockVectorData();
            timeRunFSI += dtFSI;
            dtFSISubstepping = 0.;
        }
        
#endif
        
        if (flagConver) {
            timeRun = _timeSteppingPlan.getTimeRun();
            dt = _timeSteppingPlan.getTimeStep();
            if (_comRank == 0) {
                INFO("Time %g %g; time step %g", timeRun, totalTime,dt);
            }
            nStep++;
            sampler2++;
            savingTimeIntervalRegister += dt;
        } else {// Decreasing the time step in case of not convergence
            savingTimeIntervalRegister-=dt;
            _timeSteppingPlan.modifyTimeStepByFactor(0.5);
            timeRun = _timeSteppingPlan.getTimeRun();;
            dt =  _timeSteppingPlan.getTimeStep();
            savingTimeIntervalRegister+=dt;
            flagConver = true;
            numFailedSteps ++;
            if (numFailedSteps > _maxFailedSteps )
            {
                ERROR("maximum number of failed steps %d is reached, the simulation cannot continue",_maxFailedSteps);
                timeRun = _timeSteppingPlan.getLastTime();
                break;
            }
            if (_comRank == 0) {
                INFO("Time after resetting %g %g; time step %g", timeRun, totalTime,dt);
            }
            tangentEstimation = true; // tangent is re-estimated if time step is reduced
        }
#ifdef FSI
        dtFSISubstepping += dt;
#endif  
        //if flagActivation is true, the function to update the status of the GPs (active/inactive) is called and vectors and matrices using to solve the problem are resized
        if (_flagActivation) 
        {
            vector<int> nActiveDofOld = _nActiveDof; 
            updateActivation(GPs, timeRun, _activeDofMapLocal, _nActiveDof, _ISactiveDof, nNodGlobal*numMechDofPerNode, numExtraDof,
                             nodLocal, nNodGlobal, nNodLocal, _comSize, _comRank);
            
            int diff = 0;
            for (int ii=0; ii < _nActiveDof.size(); ii++)
            {
                diff += abs(nActiveDofOld[ii] - _nActiveDof[ii]);
            }
            if (diff != 0)
            {
                if (numMechDofPerNode == _ndim)
                {
                    activeNumberDispDofs = _nActiveDof[0];
                    activeNumberStochasticDofs = 0;
                    activeNumberExtraDofs = _nActiveDof[1];
                }
                else
                {
                    activeNumberDispDofs = _nActiveDof[0];
                    activeNumberStochasticDofs = _nActiveDof[1];
                    activeNumberExtraDofs = 0;
                }
                nlsys->clearSystem();
                nlsys->allocateSystem(_activeDofMapLocal,_ISactiveDof, 
                            activeNumberDispDofs, activeNumberStochasticDofs, activeNumberExtraDofs);
                tangentEstimation = true;
                
            }
        };
        nlsys->dynamicPrediction(timeRun,dt, scatter_uK);
        // apply displacement BC
        applyDisplacementBC(nlsys, _numNodes*numMechDofPerNode, timeRun, dt,nodLocal, nNodGlobal,nNodLocal);
        // get unknown from system
        nlsys->getFromSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allField);
        nlsys->getFromAccelerationSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allAcc);
        
        ///////////////////////////////////// NEWTON RAPHSON LOOP
        ite2 = 0; // starting iteration counting
        // initial residual
        double normResInit(0.), normResExtraInit(0.), normResStochInit(0.); // residual at iteratioon 0
        double normRes(0.), normResExtra(0.), normResStoch(0.); // residual at iteratioon k
        
        double timeCPU = getCPUTime();
        double timeWall = getTimeOfDay();
        
        while (true) 
        {
            ite2++;
            // check if tangent need to be estimated
            if (not(nStep ==1 and ite2 ==1))
            {
                if (_stiffnessMatrixType == solvers::Current)
                {
                    tangentEstimation = true;
                }
                else if (_stiffnessMatrixType == solvers::Beginning)
                {
                    tangentEstimation = false;
                    if ( _flagActivation and ite2 ==1)
                    {
                        tangentEstimation = true;
                    }
                }
                else if (_stiffnessMatrixType == solvers::Initial)
                {
                    if (ite2 ==1)
                    {
                        tangentEstimation = true;
                    }
                    else
                    {
                        tangentEstimation = false;
                    }
                }
                else if (_stiffnessMatrixType == solvers::Iteration)
                {
                    if (ite2 <=_stiffnessMatrixFixedIteration)
                    {
                        tangentEstimation = true;
                    }
                    else
                    {
                        tangentEstimation = false;
                    }
                }
                else if (_stiffnessMatrixType == solvers::Interval)
                {
                    if ((ite2-1)%_stiffnessMatrixFixedIteration==0)
                    {
                        tangentEstimation = true;
                    }
                    else
                    {
                        tangentEstimation = false;
                    }
                }
                else
                {
                    ERROR("This option has not been implemented yet!");
                }
            }
            
            // 
            // compute constitutive
            computeState(GPs,allField,dt,timeRun,tangentEstimation);
            // compute resudial
            double stepRatio = 1;
#ifdef FSI
            stepRatio = dtFSISubstepping/dtFSI;
#endif
            map<pair<int,int>,double> stiffBC;
            computeRightHandSide(nlsys, inerForce, intForce,extForce,totalNumberDofsPerNode,allField, allAcc, 
                                GPs,  nod, elem, NeumannBCs, _ndim, timeRun, dt, 
                                stiffBC,tangentEstimation,stepRatio);
            // add to RHS
            vector<double> normResall;
            nlsys->getNormInfRHS(normResall);
            
            normRes = normResall[0];
            normResStoch = normResall[1];
            normResExtra = normResall[2];
            // CONVERGENCE CRITERIA : Check of residual norm
            if ((!std::isfinite(normRes)) || (!std::isfinite(normResExtra)) || (!std::isfinite(normResStoch))) 
            {
                if (_comRank == 0)
                {
                    if (extraDofFlag and StochasticFlag) 
                    {
                        INFO("Residuals are not finite numbers : Mech-atol %.16g  Extra-atol %.16g Stoch-atol %.16g", normRes, normResExtra,normResStoch);
                    }
                    else if (extraDofFlag)
                    {
                        INFO("Residuals are not finite numbers : Mech-atol %.16g  Extra-atol %.16g", normRes, normResExtra);
                    }
                    else if (StochasticFlag)
                    {
                        INFO("Residuals are not finite numbers : Mech-atol %.16g Stoch-atol %.16g", normRes,normResStoch);
                    }
                    else
                    {
                        INFO("Residuals are not finite numbers : Mech-atol %.16g", normRes);
                    }
                }
                flagConver = false;
                // back to previous data
                resetToPreviousStep(GPs,nlsys);
                break;
            }
            else
            {
                double tolRes=0.;
                double tolExtra =0.;
                double tolStoch =0.;
                
                if (normResInit <= _absTol)
                {
                    normResInit = normRes;
                }
                if (normResInit >_absTol)
                {
                    tolRes = normRes/normResInit;
                }
                if (normResExtraInit <= _absTol)
                {
                    normResExtraInit = normResExtra;
                }
                if (normResExtraInit > _absTol)
                {
                    tolExtra = normResExtra/normResExtraInit;
                }
                if (normResStochInit <= _absTol)
                {
                    normResStochInit = normResStoch;
                }
                if (normResStochInit >_absTol)
                {
                    tolStoch= normResStoch/normResStochInit;
                }
                
                if (_comRank == 0) {
                    if (extraDofFlag and StochasticFlag) 
                    {
                        INFO("Iteration %i: Mech-atol %e Mech-rtol %e  Extra-atol %e Extra-rtol %e Stoch-atol %e Stoch-rtol %e", 
                                    ite2, normRes, tolRes, normResExtra,tolExtra ,normResStoch,tolStoch);
                    }
                    else if (extraDofFlag)
                    {
                        INFO("Iteration %i: Mech-atol %e Mech-rtol %e  Extra-atol %e Extra-rtol %e", 
                                    ite2, normRes, tolRes, normResExtra,tolExtra);
                    }
                    else if (StochasticFlag)
                    {
                        INFO("Iteration %i: Mech-atol %e Mech-rtol %e Stoch-atol %e Stoch-rtol %e", 
                                    ite2, normRes, tolRes, normResStoch,tolStoch);
                    }
                    else
                    {
                        INFO("Iteration %i: Mech-atol %e Mech-rtol %e", 
                                    ite2, normRes, tolRes);
                    }
                }
                
                if ((normRes < _absTol) && (normResExtra < _absTol) && (normResStoch < _absTol))
                {
                    if (_comRank == 0)
                    {
                        if (extraDofFlag and StochasticFlag) 
                        {
                            INFO("Converged Iteration %i on absolute residual: Mech-atol %.16g  Extra-atol %.16g Stoch-atol %.16g", ite2, normRes, normResExtra,normResStoch);
                        }
                        else if (extraDofFlag)
                        {
                            INFO("Converged Iteration %i on absolute residual: Mech-atol %.16g  Extra-atol %.16g", ite2, normRes, normResExtra);
                        }
                        else if (StochasticFlag)
                        {
                            INFO("Converged Iteration %i on absolute residual: Mech-atol %.16g Stoch-atol %.16g", ite2, normRes,normResStoch);
                        }
                        else
                        {
                            INFO("Converged Iteration %i on absolute residual: Mech-atol %.16g", ite2, normRes);
                        }
                        INFO("Time for step %i: Wall %gs, CPU %gs",nStep,getTimeOfDay()-timeWall,getCPUTime()-timeCPU);
                    }
                    flagConver = true;
                    break;
                }
                else
                {
                    if (normRes < _absTol)
                    {
                        tolRes = 0;
                    }
                    if (normResExtra < _absTol)
                    {
                        tolExtra = 0.;
                    }
                    if (normResStoch < _absTol)
                    {
                        tolStoch = 0.;
                    }
                    
                    double tol = std::max(tolRes,std::max(tolExtra,tolStoch));
                    if (tol < _relTol)
                    {
                        if (_comRank == 0) {
                            INFO("Converged Solution at Iteration %i, total computation time in step %i: Wall %gs, CPU %gs",ite2,nStep,getTimeOfDay()-timeWall,getCPUTime()-timeCPU);
                        }
                        flagConver = true;
                        break;
                    }
                }
            }
            
            
            if (tangentEstimation)
            {
                computeStiffness(timeRun, dt, allField,GPs,nlsys);
                if (stiffBC.size() > 0)
                {
                    for (map<pair<int,int>,double>::const_iterator itMa = stiffBC.begin(); itMa !=stiffBC.end(); itMa++)
                    {
                        const pair<int,int>& pos = itMa->first;
                        double val = itMa->second;
                        int rowR = _activeDofMapLocal[pos.first];
                        int colR = _activeDofMapLocal[pos.second];
                        if (rowR > -1 and colR > -1) 
                        {
                            nlsys->getMatrix()->addIJ(rowR,colR,-1.*val); 
                        }
                    }
                }
            }
            nlsys->getMatrix()->assembly();
            
            // solve system
            flagConver = nlsys->systemSolve(timeRun,dt, _petsc_solver, _petsc_precon);
            
            if ((!flagConver) or (ite2 == _maxIterations))
            {
                resetToPreviousStep(GPs,nlsys);
                if (ite2 == _maxIterations) 
                {
                    ERROR("ERROR: Too many iterations in Newmark Time=%g", timeRun);
                }
                flagConver = false;
                break; // for new smaller time step
            }
            else
            {
                // update field
                nlsys->getFromSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allField);
                nlsys->getFromAccelerationSolutionByRange(0,_numNodes*totalNumberDofsPerNode,allAcc);
            }
            
            if (GeneralOptions::noImplicitIterative)
            {
                if (_comRank==0)
                {
                    INFO("iteration is stopped from options");
                }
                
                break;
            }
            
        } // End of Newton Loop

#ifdef FSI
        if (!flagConver)
        {
          dtFSISubstepping -= dt;
        }
        
        if (flagConver && fabs(dtFSI- dtFSISubstepping) <1e-12)
        {
          FSIstepDone = true;
        }
        else
        {
          FSIstepDone = false;
        }
  
        if (_flagFSI && FSIstepDone)
        {
            _adapter->updateDisplacementInterface(dtFSI,uK,_numNodes);
            dtFSI = _adapter->advance(dtFSI);
            if(_adapter->isActionRequired(_adapter->actionReadIterationCheckpoint()))
            {
                restoreACheckPoint(GPs,nlsys);
                _adapter->restoreACheckPoint(timeRunFSI,dtFSI);
                _adapter->markActionFulfilled(_adapter->actionReadIterationCheckpoint());
                _timeSteppingPlan.setTime(timeRunFSI,dtFSI);
                continue;
            }; 
        }
#endif  
        if (!flagConver)
        {
            continue;
        }

        // get solution
        getSolutionsFormSystem(nlsys, _numNodes*numMechDofPerNode, numExtraDof,uK, vK, aK, vvK);
        // save to file
        if (savingTimeIntervalRegister >= savingTimeInterval || sampler2 == 0 || (fabs(timeRun -totalTime) < 1E-12*totalTime)) {
            double startTimeF = getCPUTime();
            double startWallF = getTimeOfDay();

            if (extraDofFlag) {
                outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, vvK,
                                  extraDofAdditionalOutput, sim, _numDof, GPs, _outputs, dt, timeRun,
                                  _comRank);
            } else {
                outputsManagement(nod, elem, _ndim, nodFEM, nodMM, nodShare, elemFEM, elemMM, uK, vK, aK, sim,
                                  _numDof, GPs, _outputs, dt, timeRun, _comRank);
            }
            sampler2 = 0;
            savingTimeIntervalRegister=0;
            if (_comRank==0)
                INFO("Time to call outputsManagement wall %g cpu %g",getTimeOfDay()-startWallF,getCPUTime()-startTimeF);
        }
        
        printForcesOut(intForce, totalDofs, dt, inerForce, 0, _forceOutputs, _t0, _tf, timeRun, _ndim, _comRank, _numDof,
                       nNodGlobal, nNodLocal, nodLocal, _comSize);
        
        classExtractData::extractField(GPs,uK,vK,vvK, intForce, extForce, _numNodes,timeRun);
        if (_trueDispField != NULL)
        {
            double L2Error = _trueDispField->getL2Error(allField,GPs,nod);
            double H1Error = _trueDispField->getH1Error(allField,GPs,nod);
            INFO("error analysis was carried out: L2-error = %g H1-error %g",L2Error,H1Error);
            
            string fname= "";
            if (GeneralOptions::commSize==1)
            {
                fname = GeneralOptions::getOutputDirectoryPath()+"/errorData_step"+ std::to_string(nStep) +".csv";
            }
            else
            {
                fname = GeneralOptions::getOutputDirectoryPath()+"/errorData_step"+ std::to_string(nStep)+"_rank"+std::to_string(GeneralOptions::commRank)+".csv";
            }
            FILE* pFile = fopen(fname.c_str(),"w");
            fprintf(pFile,"Number Of Nodes, Number of Elements, L2 error, H1 error\n");
            fprintf(pFile,"%ld,%ld,%g,%g", nod.size(),elem.size(),L2Error,H1Error);
            fclose(pFile);
        }
        
        nextStep(GPs,nlsys,false);
        _timeSteppingPlan.nextStep();
        if (_adaptiveTimeStep != NULL)
        {
            double dtNew = _adaptiveTimeStep->getTimeStep(ite2, dt);
            _timeSteppingPlan.setTimeStep(dtNew);
        }
        
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
    }    // end of the time loop
    
    
    // For multisolvers
    _u = uK;
    _v = vK;
    _a = aK;
    _vv = vvK;
    _extF = extForce;
    
    
#ifdef PARALLEL
    if (_comRank == 0)
#endif
    {
        INFO("done solving (Wall %gs, CPU %gs)", getTimeOfDay() - startWall, getCPUTime() - startTime);
    }
#ifdef FSI
    if (_flagFSI)
    {
        _adapter->finalize();
        delete _adapter;
    }
#endif
    delete nlsys;
};


/*! \brief Constructor of the class
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of mechanical degrees of fredom
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
classSolverImplicitStatic::classSolverImplicitStatic(vector<classNodes *> &nod, int ndim, int numDof,
                                                     vector<classDirichlet *> &DNod, double scaleFactor,
                                                     classOutputs *&outputs, double t0, double tf, string solverType,
                                                     const vector<int> &activeDofMapLocal, const vector<int>& nActiveDof,
                                                     const vector<int> &ISactiveDof, classPrintForces *&forceOutputs,
                                                     bool flagActivation, string &petsc_solver, string &petsc_precon, bool flagFSI)
        : classSolverImplicitBase(nod, ndim, numDof, DNod,
                  scaleFactor, outputs, t0, tf,
                  solverType, activeDofMapLocal,
                  nActiveDof, ISactiveDof,
                  forceOutputs, flagActivation,
                  petsc_solver, petsc_precon, flagFSI) {

}

void classSolverImplicitStatic::setLastSolutionsToSystem(nonlinearSystem* nlsys, 
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
    }
    
    if (scatter_vvK.size() >0)
    {
        nlsys->setToSolution(scatter_vvK,vvK);
    }
}
void classSolverImplicitStatic::getSolutionsFormSystem(nonlinearSystem* nlsys, 
                                int numMechDofs,  int numExtraDof,
                                vector<double>& uK,
                                vector<double>& vK,
                                vector<double>& aK,
                                vector<double>& vvK) const
                                
{
    double startTime = getCPUTime();
    double startWall = getTimeOfDay();
    
    if (numMechDofs > 0)
    {
        nlsys->getFromSolutionByRange(0,numMechDofs,uK);
    }
    vK.resize(uK.size());
    setAll(vK,0);
    aK.resize(uK.size());
    setAll(aK,0);
    if (numExtraDof >0)
    {
        nlsys->getFromSolutionByRange(numMechDofs,numMechDofs+numExtraDof,vvK);
    }
    
    if (_comRank==0)
    {
        INFO("Time to call getSolutionsFormSystem wall %g cpu %g",getTimeOfDay()-startWall,getCPUTime()-startTime);
    }
    
}
void classSolverImplicitStatic::computeRightHandSide(nonlinearSystem* nlsys, vector<double>&inerForce, vector<double>&intForce, vector<double>&extForce,
                            int totalNumberDofsPerNode, const vector<double> &allField, const vector<double>& allAcc, 
                            vector<classGPs *> &GPs, 
                            vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, double timeRun,
                            double dt, map<pair<int,int>,double>& stiffBC, bool tangentEstimation, double FSIstepRatio)
{
    double startTime = getCPUTime();
    double startWall = getTimeOfDay();
    
    setAll(intForce,0);
    setAll(extForce,0);
    setAll(inerForce,0);
    computeInternalForce(GPs,intForce,timeRun,dt);
    NeumannBCManagement(allField, GPs, nod, elem, NeumannBCs, ndim, extForce, timeRun, dt, tangentEstimation, stiffBC);
        // 
#ifdef FSI
    if (_flagFSI)
    {
        _adapter->updateForceInterface(extForce,_numNodes, FSIstepRatio);
    }
#endif      
    nlsys->zeroRHS();
    for (int row=0; row< _numNodes*totalNumberDofsPerNode; row++)
    {
        double val = intForce[row] - extForce[row];
        nlsys->addToRHS(row,val);
    }
     if (_comRank==0)
    {
        INFO("Time to call computeRightHandSide wall %g cpu %g",getTimeOfDay()-startWall,getCPUTime()-startTime);
    }
};

void classSolverImplicitStatic::applyDisplacementBC(nonlinearSystem* nlsys, int numMechDofs, double timeRun, double dt, 
                        const vector<int> &nodLocal, int nNodGlobal, int nNodLocal)
{
    // apply dirichlet BC
    double U_diri=0;
    //Restore the dirichlet imposed displacements
    for (int i = 0; i < _DNod.size(); i++) {
        int ind = _DNod[i]->getMe();
        U_diri = _DNod[i]->getUU();
        int flagDirichlet = _DNod[i]->getflagDiri();
        if (flagDirichlet == 1) { //if flagDirichlet=1, Dirichlet BC are applied instantaneously and are applied just one time
            U_diri = U_diri;
            _DNod[i]->resetflagDiri();
        } else if (flagDirichlet == -1) { //if flagDirichlet=-1, Dirichlet BC are applied as ramp from t0 to tf
            U_diri = U_diri * dt;
        } else {
            U_diri = 0;
        }
        int dofGlobal = ind;
#ifdef PARALLEL
        locDofToGlobDof(nNodGlobal, nNodLocal, ind, dofGlobal, nodLocal, false);
#endif
        nlsys->addToSolution(dofGlobal,U_diri);
    };
}

/*! \brief Constructor of the class
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of mechanical degrees of fredom
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
classSolverImplicit::classSolverImplicit(vector<classNodes *> &nod, int ndim, int numDof,
                                         vector<classDirichlet *> &DNod, double scaleFactor, classOutputs *&outputs,
                                         double t0, double tf, string solverType, const vector<int> &activeDofMapLocal,
                                         const vector<int>& nActiveDof, const vector<int> &ISactiveDof,
                                         classPrintForces *&forceOutputs, bool flagActivation,
                                         string &petsc_solver, string &petsc_precon, bool flagFSI) : 
                                         classSolverImplicitBase(nod, ndim,
                                               numDof, DNod,
                                               scaleFactor,
                                               outputs, t0,
                                               tf, solverType,
                                               activeDofMapLocal,
                                               nActiveDof,
                                               ISactiveDof,
                                               forceOutputs,
                                               flagActivation,
                                               petsc_solver,
                                               petsc_precon, flagFSI) {

}

void classSolverImplicit::setLastSolutionsToSystem(nonlinearSystem* nlsys, 
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
void classSolverImplicit::getSolutionsFormSystem(nonlinearSystem* nlsys, 
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
void classSolverImplicit::computeRightHandSide(nonlinearSystem* nlsys, vector<double>&inerForce, vector<double>&intForce, vector<double>&extForce,
                            int totalNumberDofsPerNode, const vector<double> &allField, const vector<double>& allAcc, 
                            vector<classGPs *> &GPs, 
                            vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, double timeRun,
                            double dt, map<pair<int,int>,double>& stiffBC, bool tangentEstimation, double FSIstepRatio)
{
    
    setAll(intForce,0);
    setAll(extForce,0);
    setAll(inerForce,0);
    computeInternalForce(GPs,intForce,timeRun,dt);
    NeumannBCManagement(allField, GPs, nod, elem, NeumannBCs, ndim, extForce, timeRun, dt, 
                        tangentEstimation, stiffBC);
    computeInertialForce(GPs, allAcc,inerForce);
        // 
#ifdef FSI
    if (_flagFSI)
    {
        _adapter->updateForceInterface(extForce,_numNodes, FSIstepRatio);
    }
#endif      
    nlsys->zeroRHS();
    for (int row=0; row< _numNodes*totalNumberDofsPerNode; row++)
    {
        double val = inerForce[row]+ intForce[row] - extForce[row];
        nlsys->addToRHS(row,val);
    }
};

void classSolverImplicit::applyDisplacementBC(nonlinearSystem* nlsys, int numMechDofs, double timeRun, double dt, 
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