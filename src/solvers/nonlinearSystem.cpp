//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file nonlinearSystem.cpp
  \brief This file contains all functions related to the system
*/

#include "nonlinearSystem.h"
#include "maths.h"
#include "solvers.h"
#include "commandLine.h"

ImplicitStaticNonlinearSystem::ImplicitStaticNonlinearSystem() : _solutionAllocated(false), 
                _systemAllocated(false), _stiffnessMatrix(NULL), 
                _matrixModified(false), _kspIsAllocated(false), _solutionModified(false),
                _tmpStateAllocated(false){}
ImplicitStaticNonlinearSystem::~ImplicitStaticNonlinearSystem()
{ 
    clearSolution();
    clearSystem();
    if (_kspIsAllocated) KSPDestroy(&_ksp);
    _kspIsAllocated = false;
}

void ImplicitStaticNonlinearSystem::allocateSolution(int totalDofs, const vector<int>& scatterIS)
{
    if (_solutionAllocated)
    {
        clearSolution();
    }
    _solutionAllocated = true;
    _solutionModified = false;
    
    _call(VecCreate(PETSC_COMM_WORLD, &uKPetsci));
    _call(VecSetSizes(uKPetsci, PETSC_DECIDE, totalDofs));
    _call(VecSetFromOptions(uKPetsci));
    _call(VecDuplicate(uKPetsci, &uKPetsciPrev));
    _call(VecZeroEntries(uKPetsci));
    _call(VecZeroEntries(uKPetsciPrev));
#ifdef PARALLEL
    _scatterIS = scatterIS;
    _call(VecCreateSeq(PETSC_COMM_SELF, scatterIS.size(), &uKPetsci_Seq));
    IS From;
    _call(ISCreateGeneral(PETSC_COMM_SELF, scatterIS.size(), &_scatterIS[0], PETSC_COPY_VALUES, &From));
    _call(VecScatterCreate(uKPetsci, From, uKPetsci_Seq, NULL, &_uKScatter));
    _call(VecZeroEntries(uKPetsci_Seq));
#endif            
    
}
void ImplicitStaticNonlinearSystem::clearSolution()
{
    _solutionAllocated = false;
    _solutionModified = true;
    _call(VecDestroy(&uKPetsci));
    _call(VecDestroy(&uKPetsciPrev));
#ifdef PARALLEL
    _call(VecScatterDestroy(&_uKScatter));
    _call(VecDestroy(&uKPetsci_Seq));
    _scatterIS.clear();
#endif            
    if (_tmpStateAllocated)
    {
        _call(VecDestroy(&uKPetsciTmp));
#ifdef PARALLEL
        _call(VecDestroy(&uKPetsciTmp_Seq));
#endif            
    _tmpStateAllocated = false;
    }

}

void ImplicitStaticNonlinearSystem::allocateSystem(const vector<int>& activeDofMapLocal,
                             const vector<int>& ISactiveDof,
                             int activeDispNumDofs,
                             int activeStochasticNumDofs,
                             int activeExtraNumDofs)
{
    if (_systemAllocated)
        clearSystem();
    _systemAllocated = true;
    _activeDofMapLocal = activeDofMapLocal;
    _activeDispNumDofs = activeDispNumDofs; // 
    _activeStochasticNumDofs = activeStochasticNumDofs; //
    _activeExtraNumDofs = activeExtraNumDofs;  //
    int totalActiveDofs = activeDispNumDofs+activeStochasticNumDofs+activeExtraNumDofs;
    
    int rank=0;
#ifdef PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (rank == 0)
    {
        INFO("allocate a system of %d DOFs, mech %d, stoch %d extradof %d",totalActiveDofs,activeDispNumDofs,activeStochasticNumDofs,activeExtraNumDofs);
    }
    // allocate system 
    // create residual
    _call(VecCreate(PETSC_COMM_WORLD, &rhsTotalPetsci));
    _call(VecSetSizes(rhsTotalPetsci, PETSC_DECIDE, totalActiveDofs));
    _call(VecSetFromOptions(rhsTotalPetsci));
    // create solution
    _call(VecDuplicate(rhsTotalPetsci, &deltaTotalDofPetsci));
    //  matrix
    _matrixModified = true;
    if (withStiffnessMatrix())
    {
        _stiffnessMatrix = new classCRSMatrix(totalActiveDofs, totalActiveDofs, "CRS");
    }
    //
    if (activeDispNumDofs >0)
    {
        _call(VecCreate(PETSC_COMM_WORLD, &rhsDispPetsci));
        _call(VecSetSizes(rhsDispPetsci, PETSC_DECIDE, activeDispNumDofs));
        _call(VecSetFromOptions(rhsDispPetsci));
    }
    
    if (activeStochasticNumDofs >0)
    {
        _call(VecCreate(PETSC_COMM_WORLD, &rhsStochasticPetsci));
        _call(VecSetSizes(rhsStochasticPetsci, PETSC_DECIDE, activeStochasticNumDofs));
        _call(VecSetFromOptions(rhsStochasticPetsci));
        
    }
    
    if (activeExtraNumDofs >0)
    {
        _call(VecCreate(PETSC_COMM_WORLD, &rhsExtraPetsci));
        _call(VecSetSizes(rhsExtraPetsci, PETSC_DECIDE, activeExtraNumDofs));
        _call(VecSetFromOptions(rhsExtraPetsci));
    }
    
    
#ifdef PARALLEL
    int istart(0), iend(totalActiveDofs);
    _call(VecGetOwnershipRange(rhsTotalPetsci, &istart, &iend));
    ISMapLocal.resize(iend - istart,0);
    for (int j = 0; j < (iend - istart); j++) 
    {
        ISMapLocal[j] = ISactiveDof[istart + j];
    };
#else
    ISMapLocal.resize(totalActiveDofs);
    for (int i=0; i< activeDofMapLocal.size(); i++)
    {
        if (activeDofMapLocal[i] > -1)
        {
            ISMapLocal[activeDofMapLocal[i]] = i;
        }
    }
#endif

};

void ImplicitStaticNonlinearSystem::clearSystem()
{
    _systemAllocated = false;
    _call(VecDestroy(&rhsTotalPetsci));
    _call(VecDestroy(&deltaTotalDofPetsci));
    if (withStiffnessMatrix())
    {
        delete _stiffnessMatrix; _stiffnessMatrix = NULL;
    }
    _matrixModified = true;
    _activeDofMapLocal.clear();
    if (_activeDispNumDofs >0)
        _call(VecDestroy(&rhsDispPetsci));
    if (_activeStochasticNumDofs >0)
        _call(VecDestroy(&rhsStochasticPetsci));
    if (_activeExtraNumDofs >0)
        _call(VecDestroy(&rhsExtraPetsci));

    _activeDispNumDofs = 0; // 
    _activeStochasticNumDofs = 0; //
    _activeExtraNumDofs = 0;  //
    ISMapLocal.clear(); 
}

void ImplicitStaticNonlinearSystem::zeroRHS()
{
    _call(VecAssemblyBegin(rhsTotalPetsci));
    _call(VecAssemblyEnd(rhsTotalPetsci));
    _call(VecZeroEntries(rhsTotalPetsci));
    if (_activeDispNumDofs >0)
    {
        _call(VecAssemblyBegin(rhsDispPetsci));
        _call(VecAssemblyEnd(rhsDispPetsci));
        _call(VecZeroEntries(rhsDispPetsci));
        
    }
    if (_activeStochasticNumDofs>0)
    {
        _call(VecAssemblyBegin(rhsStochasticPetsci));
        _call(VecAssemblyEnd(rhsStochasticPetsci));
        _call(VecZeroEntries(rhsStochasticPetsci));
        
    }
    if (_activeExtraNumDofs>0)
    {
        _call(VecAssemblyBegin(rhsExtraPetsci));
        _call(VecAssemblyEnd(rhsExtraPetsci));
        _call(VecZeroEntries(rhsExtraPetsci));
        
    }
}

void ImplicitStaticNonlinearSystem::addToRHS(int row, double val)
{
    int rowR = _activeDofMapLocal[row];
    if (rowR > -1)
    {
        _call(VecSetValue(rhsTotalPetsci,rowR,val,ADD_VALUES));
        if (rowR < _activeDispNumDofs)
        {
            _call(VecSetValue(rhsDispPetsci,rowR,val,ADD_VALUES));
        }
        else if (rowR < _activeDispNumDofs+ _activeStochasticNumDofs)
        {
            int rowSto = rowR-_activeDispNumDofs;
            _call(VecSetValue(rhsStochasticPetsci,rowSto,val,ADD_VALUES));
        }
        else
        {
            int rowExtra = rowR-_activeDispNumDofs - _activeStochasticNumDofs;
            _call(VecSetValue(rhsExtraPetsci,rowExtra,val,ADD_VALUES));
        };
    };
};

void ImplicitStaticNonlinearSystem::getNormInfRHS(vector<double>& norm)
{
    norm.resize(3);
    norm[0] = 0;
    if (_activeDispNumDofs > 0)
    {
        _call(VecAssemblyBegin(rhsDispPetsci));
        _call(VecAssemblyEnd(rhsDispPetsci));
        _call(VecNorm(rhsDispPetsci, NORM_INFINITY, &norm[0]));
    }
    
    norm[1] = 0;
    if (_activeStochasticNumDofs > 0)
    {
        _call(VecAssemblyBegin(rhsStochasticPetsci));
        _call(VecAssemblyEnd(rhsStochasticPetsci));
        _call(VecNorm(rhsStochasticPetsci, NORM_INFINITY, &norm[1]));
    }
    
    norm[2] = 0;
    
    if (_activeExtraNumDofs >0)
    {
        _call(VecAssemblyBegin(rhsExtraPetsci));
        _call(VecAssemblyEnd(rhsExtraPetsci));
        _call(VecNorm(rhsExtraPetsci, NORM_INFINITY, &norm[2]));
    };
}

void ImplicitStaticNonlinearSystem::zeroMatrix()
{
    _matrixModified = true;
    _stiffnessMatrix->zeroEntries();
}

void ImplicitStaticNonlinearSystem::addToMatrix(int row, int col, double val)
{
    int rowR = _activeDofMapLocal[row];
    int colR = _activeDofMapLocal[col];
    if (rowR > -1 and colR > -1) 
    {
        _stiffnessMatrix->addIJ(rowR,colR,val);
         _matrixModified = true;
    }
}

bool ImplicitStaticNonlinearSystem::systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon)
{
    double timeCPU = getCPUTime();
    double timeWall = getTimeOfDay();
    // call linear system to find deltaTotalDofPetsci
    int value = _matrixModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _matrixModified = (value==1);
    
    if (not(_kspIsAllocated) or  _matrixModified)
    {
        if (_kspIsAllocated)
        {
            _call(KSPDestroy(&_ksp));
        }
        _kspIsAllocated = true;
        PC pc;           /* preconditioner context */
        Mat A;
        _stiffnessMatrix->getPetscMat(A);
        _call(KSPCreate(PETSC_COMM_WORLD, &_ksp));
        _call(KSPSetOperators(_ksp, A, A));
        _call(KSPGetPC(_ksp, &pc));
        /// Assigning the solver
        if ((petsc_solver == "lu") or (petsc_solver == "cholesky")) // Direct solvers
        {
            _call(KSPSetType(_ksp,KSPPREONLY));
            _call(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
            _call(PCSetType(pc, petsc_solver.c_str()));
        } 
        else
        {
            _call(KSPSetType(_ksp, petsc_solver.c_str())); // Iterative solvers
            _call(PCSetType(pc, petsc_precon.c_str()));
        }
        _call(PCSetFromOptions(pc));
        _call(KSPSetTolerances(_ksp, 1E-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
        _call(KSPSetReusePreconditioner(_ksp,PETSC_TRUE));
        _call(PCSetReusePreconditioner(pc,PETSC_TRUE));
	_call(KSPSetFromOptions(_ksp));                
    }
    
    _call(VecAssemblyBegin(rhsTotalPetsci));
    _call(VecAssemblyEnd(rhsTotalPetsci));
    _call(KSPSolve(_ksp, rhsTotalPetsci, deltaTotalDofPetsci));
    
    KSPConvergedReason reason;
    _call(KSPGetConvergedReason(_ksp, &reason));
    if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "PETSC has failed\n");
        PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
        ERROR("The linear system of equations did not converge");
        return false;
    }
    // update
    _call(VecScale(deltaTotalDofPetsci,-1.));
    double* localData;
    _call(VecGetArray(deltaTotalDofPetsci,&localData));
    int sizeIS = ISMapLocal.size();
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    _call(VecSetValues(uKPetsci, sizeIS, &ISMapLocal[0], localData, ADD_VALUES));
    _call(VecRestoreArray(deltaTotalDofPetsci, &localData));
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    //_call(VecView(uKPetsci,PETSC_VIEWER_STDOUT_WORLD));
    _solutionModified = true;
    _matrixModified = false;
    
    if (GeneralOptions::commRank==0)
        INFO("Time for solving linear system (Wall %gs, CPU %gs)",getTimeOfDay()-timeWall,getCPUTime()-timeCPU);
    return true;
}

void ImplicitStaticNonlinearSystem::setToSolution(int row, double val)
{
    _call(VecSetValues(uKPetsci, 1, &row, &val, INSERT_VALUES));
    _solutionModified = true;
}
void ImplicitStaticNonlinearSystem::setToSolution(const vector<int>& row, const vector<double>& val)
{
    int nsize = row.size();
    _call(VecSetValues(uKPetsci, nsize, &row[0], &val[0], INSERT_VALUES));
    _solutionModified = true;
}

void ImplicitStaticNonlinearSystem::addToSolution(int row, double val) 
{
    _call(VecSetValues(uKPetsci, 1, &row, &val, ADD_VALUES));
    _solutionModified = true;
}

void ImplicitStaticNonlinearSystem::getFromSolution(const vector<int>& rows, vector<double>& uK)
{
    int rowSize = rows.size();
    uK.resize(rowSize);
#ifdef PARALLEL           
    int value = _solutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _solutionModified = (value==1);
    if (_solutionModified)
    {
        _call(VecAssemblyBegin(uKPetsci));
        _call(VecAssemblyEnd(uKPetsci));
        _call(VecScatterBegin(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _solutionModified = false;
    }
    _call(VecGetValues(uKPetsci_Seq, rowSize, &rows[0], &uK[0]));
#else
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    _call(VecGetValues(uKPetsci, rowSize, &rows[0], &uK[0]));
#endif
}

void ImplicitStaticNonlinearSystem::getFromSolutionByRange(int startIndex, int endIndex, vector<double>& uK)
{
    int rowSize = endIndex - startIndex;
    uK.resize(rowSize);
#ifdef PARALLEL  
    int value = _solutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _solutionModified = (value==1);
    if (_solutionModified)
    {
        _call(VecAssemblyBegin(uKPetsci));
        _call(VecAssemblyEnd(uKPetsci));
        _call(VecScatterBegin(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _solutionModified = false;
    }
    double* data;
    MPI_Barrier(MPI_COMM_WORLD);
    _call(VecGetArray(uKPetsci_Seq, &data));
    uK.assign(data+startIndex,data+endIndex);
    _call(VecRestoreArray(uKPetsci_Seq, &data));   
#else
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    double* data;
    _call(VecGetArray(uKPetsci, &data));
    uK.assign(data+startIndex,data+endIndex);
    _call(VecRestoreArray(uKPetsci, &data));
#endif
}


void ImplicitStaticNonlinearSystem::nextStep()
{
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    _call(VecCopy(uKPetsci,uKPetsciPrev));
}
void ImplicitStaticNonlinearSystem::resetToPreviousSolution()
{
    _call(VecCopy(uKPetsciPrev,uKPetsci));
    _matrixModified = true;
};


void ImplicitStaticNonlinearSystem::makeACheckPoint()
{
    if (!_tmpStateAllocated)
    {
        _tmpStateAllocated = true;
        _call(VecDuplicate(uKPetsci,&uKPetsciTmp));
#ifdef PARALLEL  
        _call(VecDuplicate(uKPetsci_Seq,&uKPetsciTmp_Seq));
#endif
    };
    
    _call(VecCopy(uKPetsci,uKPetsciTmp));
#ifdef PARALLEL  
    _call(VecCopy(uKPetsci_Seq,uKPetsciTmp_Seq));
#endif
}

void ImplicitStaticNonlinearSystem::restoreACheckPoint()
{
    if (!_tmpStateAllocated)
    {
        ERROR("No checkpoint exists!");
        exit(-1);
    }
     _call(VecCopy(uKPetsciTmp,uKPetsci));
#ifdef PARALLEL  
    _call(VecCopy(uKPetsciTmp_Seq,uKPetsci_Seq));
#endif
}

DynamicNonlinearSystem::DynamicNonlinearSystem() : ImplicitStaticNonlinearSystem(), _massMatrix(NULL),
        _velocitySolutionModified(false), _accelerationSolutionModified(false)
{
}
DynamicNonlinearSystem::~DynamicNonlinearSystem(){}

void DynamicNonlinearSystem::allocateSolution(int totalDofs, const vector<int>& scatterIS)
{
    ImplicitStaticNonlinearSystem::allocateSolution(totalDofs,scatterIS);
    _call(VecDuplicate(uKPetsci, &vKPetsci));
    _call(VecDuplicate(uKPetsci, &vKPetsciPrev));
    _call(VecDuplicate(uKPetsci, &aKPetsci));
    _call(VecDuplicate(uKPetsci, &aKPetsciPrev));
    _call(VecZeroEntries(vKPetsci));
    _call(VecZeroEntries(vKPetsciPrev));
    _call(VecZeroEntries(aKPetsci));
    _call(VecZeroEntries(aKPetsciPrev));
    
#ifdef PARALLEL
    _call(VecCreateSeq(PETSC_COMM_SELF, scatterIS.size(), &vKPetsci_Seq));
    _call(VecCreateSeq(PETSC_COMM_SELF, scatterIS.size(), &aKPetsci_Seq));
    IS From;
    _call(ISCreateGeneral(PETSC_COMM_SELF, scatterIS.size(), &_scatterIS[0], PETSC_COPY_VALUES, &From));
    _call(VecScatterCreate(vKPetsci, From, vKPetsci_Seq, NULL, &_vKScatter));
    _call(VecScatterCreate(aKPetsci, From, aKPetsci_Seq, NULL, &_aKScatter));
    _call(VecZeroEntries(vKPetsci_Seq));
    _call(VecZeroEntries(aKPetsci_Seq));
    
    _call(VecDuplicate(uKPetsci_Seq,&uKPetsciPrev_Seq)); // previous local field
    _call(VecDuplicate(vKPetsci_Seq,&vKPetsciPrev_Seq)); // previous local velocity
    _call(VecDuplicate(aKPetsci_Seq,&aKPetsciPrev_Seq)); // previous local acceleration
    _call(VecZeroEntries(uKPetsciPrev_Seq));
    _call(VecZeroEntries(vKPetsciPrev_Seq));
    _call(VecZeroEntries(aKPetsciPrev_Seq));
#endif            
    _velocitySolutionModified = false;
    _accelerationSolutionModified = false;
};

void DynamicNonlinearSystem::clearSolution()
{
    if (_tmpStateAllocated)
    {
        _call(VecDestroy(&vKPetsciTmp));
        _call(VecDestroy(&aKPetsciTmp));
#ifdef PARALLEL
        _call(VecDestroy(&vKPetsciTmp_Seq));
        _call(VecDestroy(&aKPetsciTmp_Seq));
#endif       
    }
    
    if (_solutionAllocated)
    {
        ImplicitStaticNonlinearSystem::clearSolution();
        _call(VecDestroy(&vKPetsci));
        _call(VecDestroy(&vKPetsciPrev));
        _call(VecDestroy(&aKPetsci));
        _call(VecDestroy(&aKPetsciPrev));
        
#ifdef PARALLEL
        _call(VecScatterDestroy(&_vKScatter));
        _call(VecDestroy(&vKPetsci_Seq));
        _call(VecScatterDestroy(&_aKScatter));
        _call(VecDestroy(&aKPetsci_Seq));
        _call(VecDestroy(&uKPetsciPrev_Seq));
        _call(VecDestroy(&vKPetsciPrev_Seq));
        _call(VecDestroy(&aKPetsciPrev_Seq));
#endif       
    };
    
};

void DynamicNonlinearSystem::allocateSystem(const vector<int>& activeDofMapLocal,
                             const vector<int>& ISactiveDof,
                             int activeDispNumDofs,
                             int activeStochasticNumDofs,
                             int activeExtraNumDofs)
{
    ImplicitStaticNonlinearSystem::allocateSystem(activeDofMapLocal,ISactiveDof,
                    activeDispNumDofs,activeStochasticNumDofs,activeExtraNumDofs);
    int totalActiveDofs = activeDispNumDofs+activeStochasticNumDofs+activeExtraNumDofs;
    if (withMassMatrix())
    {
        _massMatrix = new classCRSMatrix(totalActiveDofs, totalActiveDofs, "CRS");
    }
};

void DynamicNonlinearSystem::clearSystem()
{
    if (_systemAllocated)
    {
        ImplicitStaticNonlinearSystem::clearSystem();
        if (withMassMatrix())
        {
            delete _massMatrix; _massMatrix = NULL;
        }
    }
};

void DynamicNonlinearSystem::setToVelocitySolution(int row, double val)
{
    _call(VecSetValues(vKPetsci, 1, &row, &val, INSERT_VALUES));
    _velocitySolutionModified = true;
};
void DynamicNonlinearSystem::setToVelocitySolution(const vector<int>& row, const vector<double>& val)
{
    int nsize = row.size();
    _call(VecSetValues(vKPetsci, nsize, &row[0], &val[0], INSERT_VALUES));
    _velocitySolutionModified = true;
};
void DynamicNonlinearSystem::addToVelocitySolution(int row, double val)
{
    _call(VecSetValues(vKPetsci, 1, &row, &val, ADD_VALUES));
    _velocitySolutionModified = true;
};

void DynamicNonlinearSystem::setToAccelerationSolution(int row, double val)
{
    _call(VecSetValues(aKPetsci, 1, &row, &val, INSERT_VALUES));
    _accelerationSolutionModified = true;
};
void DynamicNonlinearSystem::setToAccelerationSolution(const vector<int>& row, const vector<double>& val)
{
    int nsize = row.size();
    _call(VecSetValues(aKPetsci, nsize, &row[0], &val[0], INSERT_VALUES));
    _accelerationSolutionModified = true;
};
void DynamicNonlinearSystem::addToAccelerationSolution(int row, double val)
{
    _call(VecSetValues(aKPetsci, 1, &row, &val, ADD_VALUES));
    _accelerationSolutionModified = true;
};


void DynamicNonlinearSystem::getFromVelocitySolution(const vector<int>& rows, vector<double>& uK)
{
    int rowSize = rows.size();
    uK.resize(rowSize);
#ifdef PARALLEL           
    int value = _velocitySolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _velocitySolutionModified = (value==1);
    if (_velocitySolutionModified)
    {
        _call(VecAssemblyBegin(vKPetsci));
        _call(VecAssemblyEnd(vKPetsci));
        _call(VecScatterBegin(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _velocitySolutionModified = false;
    }
    _call(VecGetValues(vKPetsci_Seq, rowSize, &rows[0], &uK[0]));
#else
    _call(VecAssemblyBegin(vKPetsci));
    _call(VecAssemblyEnd(vKPetsci));
    _call(VecGetValues(vKPetsci, rowSize, &rows[0], &uK[0]));
#endif
}
void DynamicNonlinearSystem::getFromVelocitySolutionByRange(int startIndex, int endIndex, vector<double>& uK)
{
    int rowSize = endIndex - startIndex;
    uK.resize(rowSize);
#ifdef PARALLEL  
    int value = _velocitySolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _velocitySolutionModified = (value==1);
    if (_velocitySolutionModified)
    {
        _call(VecAssemblyBegin(vKPetsci));
        _call(VecAssemblyEnd(vKPetsci));
        _call(VecScatterBegin(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _velocitySolutionModified = false;
    }
    double* data;
    MPI_Barrier(MPI_COMM_WORLD);
    _call(VecGetArray(vKPetsci_Seq, &data));
    uK.assign(data+startIndex,data+endIndex);
    _call(VecRestoreArray(vKPetsci_Seq, &data));   
#else
    _call(VecAssemblyBegin(vKPetsci));
    _call(VecAssemblyEnd(vKPetsci));
    double* data;
    _call(VecGetArray(vKPetsci, &data));
    uK.assign(data+startIndex,data+endIndex);
    _call(VecRestoreArray(vKPetsci, &data));
#endif
}

void DynamicNonlinearSystem::getFromAccelerationSolution(const vector<int>& rows, vector<double>& uK)
{
    int rowSize = rows.size();
    uK.resize(rowSize);
#ifdef PARALLEL           
    int value = _accelerationSolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _accelerationSolutionModified = (value==1);
    if (_accelerationSolutionModified)
    {
        _call(VecAssemblyBegin(aKPetsci));
        _call(VecAssemblyEnd(aKPetsci));
        _call(VecScatterBegin(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _accelerationSolutionModified = false;
    }
    _call(VecGetValues(aKPetsci_Seq, rowSize, &rows[0], &uK[0]));
#else
    _call(VecAssemblyBegin(aKPetsci));
    _call(VecAssemblyEnd(aKPetsci));
    _call(VecGetValues(aKPetsci, rowSize, &rows[0], &uK[0]));
#endif
}
void DynamicNonlinearSystem::getFromAccelerationSolutionByRange(int startIndex, int endIndex, vector<double>& uK)
{
    int rowSize = endIndex - startIndex;
    uK.resize(rowSize);
#ifdef PARALLEL  
    int value = _accelerationSolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _accelerationSolutionModified = (value==1);
    if (_accelerationSolutionModified)
    {
        _call(VecAssemblyBegin(aKPetsci));
        _call(VecAssemblyEnd(aKPetsci));
        _call(VecScatterBegin(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _accelerationSolutionModified = false;
    }
    double* data;
    MPI_Barrier(MPI_COMM_WORLD);
    _call(VecGetArray(aKPetsci_Seq, &data));
    uK.assign(data+startIndex,data+endIndex);
    _call(VecRestoreArray(aKPetsci_Seq, &data));   
#else
    _call(VecAssemblyBegin(aKPetsci));
    _call(VecAssemblyEnd(aKPetsci));
    double* data;
    _call(VecGetArray(aKPetsci, &data));
    uK.assign(data+startIndex,data+endIndex);
    _call(VecRestoreArray(aKPetsci, &data));
#endif
}

void DynamicNonlinearSystem::nextStep()
{
    ImplicitStaticNonlinearSystem::nextStep();
    _call(VecAssemblyBegin(vKPetsci));
    _call(VecAssemblyEnd(vKPetsci));
    _call(VecCopy(vKPetsci,vKPetsciPrev));
    _call(VecAssemblyBegin(aKPetsci));
    _call(VecAssemblyEnd(aKPetsci));
    _call(VecCopy(aKPetsci,aKPetsciPrev));
#ifdef PARALLEL
    _call(VecCopy(uKPetsci_Seq,uKPetsciPrev_Seq));
    _call(VecCopy(vKPetsci_Seq,vKPetsciPrev_Seq));
    _call(VecCopy(aKPetsci_Seq,aKPetsciPrev_Seq));
#endif
}
void DynamicNonlinearSystem::resetToPreviousSolution()
{
    ImplicitStaticNonlinearSystem::resetToPreviousSolution();
    _call(VecCopy(vKPetsciPrev,vKPetsci));
    _call(VecCopy(aKPetsciPrev,aKPetsci));
    _velocitySolutionModified = true;
    _accelerationSolutionModified = true;
#ifdef PARALLEL
    _call(VecCopy(uKPetsciPrev_Seq,uKPetsci_Seq));
    _call(VecCopy(vKPetsciPrev_Seq,vKPetsci_Seq));
    _call(VecCopy(aKPetsciPrev_Seq,aKPetsci_Seq));
#endif
}; 

void DynamicNonlinearSystem::makeACheckPoint()
{

    if (!_tmpStateAllocated)
    {
        _call(VecDuplicate(vKPetsci,&vKPetsciTmp));
        _call(VecDuplicate(aKPetsci,&aKPetsciTmp));
#ifdef PARALLEL  
        _call(VecDuplicate(vKPetsci_Seq,&vKPetsciTmp_Seq));
        _call(VecDuplicate(aKPetsci_Seq,&aKPetsciTmp_Seq));
#endif
    };
    
    _call(VecCopy(vKPetsci,vKPetsciTmp));
    _call(VecCopy(aKPetsci,aKPetsciTmp));
#ifdef PARALLEL  
    _call(VecCopy(vKPetsci_Seq,vKPetsciTmp_Seq));
    _call(VecCopy(aKPetsci_Seq,aKPetsciTmp_Seq));
#endif
    ImplicitStaticNonlinearSystem::makeACheckPoint();
    
}

void DynamicNonlinearSystem::restoreACheckPoint()
{
    if (!_tmpStateAllocated)
    {
        ERROR("No checkpoint exists!");
        exit(-1);
    }
     _call(VecCopy(vKPetsciTmp,vKPetsci));
     _call(VecCopy(aKPetsciTmp,aKPetsci));
#ifdef PARALLEL  
    _call(VecCopy(vKPetsciTmp_Seq,vKPetsci_Seq));
    _call(VecCopy(aKPetsciTmp_Seq,aKPetsci_Seq));
#endif
    ImplicitStaticNonlinearSystem::restoreACheckPoint();
}

ImplicitNonlinearSystem::ImplicitNonlinearSystem() : DynamicNonlinearSystem()
{
    _gamma = 1.0 / 2.0;// Newmark Gamma
    _beta = (_gamma + 1.0 / 2.0) * (_gamma + 1.0 / 2.0) / 4.0;
}
ImplicitNonlinearSystem::~ImplicitNonlinearSystem(){};


bool ImplicitNonlinearSystem::systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon)
{
    // call linear system to find deltaTotalDofPetsci
    int value = _matrixModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _matrixModified = (value==1);
    
    if (not(_kspIsAllocated) or  _matrixModified)
    {
        if (_kspIsAllocated)
        {
            _call(KSPDestroy(&_ksp));
        }
        _kspIsAllocated = true;
        PC pc;           /* preconditioner context */
        Mat A, M;
        _stiffnessMatrix->getPetscMat(A);
        _massMatrix->getPetscMat(M);
        
        double factor = 1. / (dt * dt * _beta);
        _call(MatAXPY(A, factor, M, SAME_NONZERO_PATTERN)); // stiffness and mass have the same pattern
        
        _call(KSPCreate(PETSC_COMM_WORLD, &_ksp));
        _call(KSPSetOperators(_ksp, A, A));
        _call(KSPGetPC(_ksp, &pc));
        /// Assigning the solver
        if ((petsc_solver == "lu") or (petsc_solver == "cholesky")) // Direct solvers
        {
            _call(KSPSetType(_ksp,KSPPREONLY));
            _call(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
            _call(PCSetType(pc, petsc_solver.c_str()));
        } 
        else
        {
            _call(KSPSetType(_ksp, petsc_solver.c_str())); // Iterative solvers
            _call(PCSetType(pc, petsc_precon.c_str()));
        }
        _call(KSPSetTolerances(_ksp, 1E-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
        _call(KSPSetReusePreconditioner(_ksp,PETSC_TRUE));
        _call(PCSetReusePreconditioner(pc,PETSC_TRUE)); 
	_call(KSPSetFromOptions(_ksp));                
    };
    
    _call(VecAssemblyBegin(rhsTotalPetsci));
    _call(VecAssemblyEnd(rhsTotalPetsci));
    _call(KSPSolve(_ksp, rhsTotalPetsci, deltaTotalDofPetsci));
    
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    _call(VecAssemblyBegin(vKPetsci));
    _call(VecAssemblyEnd(vKPetsci));
    _call(VecAssemblyBegin(aKPetsci));
    _call(VecAssemblyEnd(aKPetsci));
    
    KSPConvergedReason reason;
    _call(KSPGetConvergedReason(_ksp, &reason));
    if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "PETSC has failed\n");
        PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
        ERROR("The linear system of equations did not converge");
        return false;
    }
    // update
    _call(VecScale(deltaTotalDofPetsci,-1.));
    double* localData;
    _call(VecGetArray(deltaTotalDofPetsci,&localData));
    int sizeIS = ISMapLocal.size();
    _call(VecSetValues(uKPetsci, sizeIS, &ISMapLocal[0], localData, ADD_VALUES));
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    // velocity and acceleration
    
#ifdef PARALLEL
    int istart, iend;
    _call(VecGetOwnershipRange(deltaTotalDofPetsci, &istart, &iend));
    int sizeISMeca = 0;
    int activeMechDofs = _activeDispNumDofs+_activeStochasticNumDofs;
    if (istart < activeMechDofs)
    {
        sizeISMeca = std::min(iend,activeMechDofs)-istart;
    }
    if (sizeISMeca >0)
    {
        // velocity
        vector<double> vIncrement(localData,localData+sizeISMeca);
        scale(vIncrement,(_gamma / (dt * _beta)));
        _call(VecSetValues(vKPetsci, sizeISMeca, &ISMapLocal[0], &vIncrement[0], ADD_VALUES));
        // 
        // acceleration 
        vector<double> aIncrement(localData,localData+sizeISMeca);
        scale(aIncrement,1./(dt * dt * _beta));
        _call(VecSetValues(aKPetsci, sizeISMeca, &ISMapLocal[0], &aIncrement[0], ADD_VALUES));
    }
    _call(VecAssemblyBegin(vKPetsci));
    _call(VecAssemblyEnd(vKPetsci));
    _call(VecAssemblyBegin(aKPetsci));
    _call(VecAssemblyEnd(aKPetsci));
#else
    double*v, *a;
    _call(VecGetArray(vKPetsci, &v));
    _call(VecGetArray(aKPetsci, &a));
    int activeMechDofs = _activeDispNumDofs+_activeStochasticNumDofs;
    for (int i = 0; i < activeMechDofs; i++) 
    {
        int ind = ISMapLocal[i];
        v[ind] +=  (_gamma / (dt * _beta)) * localData[i];
        a[ind] +=  localData[i] / (dt * dt * _beta);
    };
    _call(VecRestoreArray(vKPetsci, &v));
    _call(VecRestoreArray(aKPetsci, &a));
#endif          
    _call(VecRestoreArray(deltaTotalDofPetsci, &localData));
    _velocitySolutionModified = true;
    _accelerationSolutionModified = true;
    _solutionModified = true;
    _matrixModified = false;
    return true;
};

void ImplicitNonlinearSystem::dynamicPrediction(double timeRun, double dt, const vector<int>& scatter_uK) 
{
    
#ifdef PARALLEL
    // get local data
    int value = _solutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _solutionModified = (value==1);
    if (_solutionModified)
    {
        _call(VecAssemblyBegin(uKPetsci));
        _call(VecAssemblyEnd(uKPetsci));
        _call(VecScatterBegin(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _solutionModified = false;
    }
    
    value = _velocitySolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _velocitySolutionModified = (value==1);
    if (_velocitySolutionModified)
    {
        _call(VecAssemblyBegin(vKPetsci));
        _call(VecAssemblyEnd(vKPetsci));
        _call(VecScatterBegin(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _velocitySolutionModified = false;
    }
    
    value = _accelerationSolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _accelerationSolutionModified = (value==1);
    if (_accelerationSolutionModified)
    {
        _call(VecAssemblyBegin(aKPetsci));
        _call(VecAssemblyEnd(aKPetsci));
        _call(VecScatterBegin(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _accelerationSolutionModified = false;
    }
    // update displacement
    _call(VecAXPY(uKPetsci_Seq,dt,vKPetsci_Seq));
    double fact = dt * dt * (0.5 - _beta);
    _call(VecAXPY(uKPetsci_Seq,fact,aKPetsci_Seq));
    // update velocity
    fact = dt * (1 - _gamma);
    _call(VecAXPY(vKPetsci_Seq,fact,aKPetsci_Seq));
    // update acceleration
    _call(VecZeroEntries(aKPetsci_Seq));
    _call(VecZeroEntries(aKPetsci));
    
    double *u, *v, *uprev, *vprev;
    _call(VecGetArray(uKPetsci_Seq, &u));
    _call(VecGetArray(vKPetsci_Seq, &v));
    _call(VecGetArray(uKPetsciPrev_Seq, &uprev));
    _call(VecGetArray(vKPetsciPrev_Seq, &vprev));
    for ( int i=0; i< scatter_uK.size(); i++)
    {
        int rowR = _activeDofMapLocal[i];
        if (rowR <0)
        {
            // dirichlet nodes
            u[i] = uprev[i];
            v[i] = vprev[i];
        };
    };
    _call(VecSetValues(uKPetsci, scatter_uK.size(), &scatter_uK[0], u, INSERT_VALUES));
    _call(VecSetValues(vKPetsci, scatter_uK.size(), &scatter_uK[0], v, INSERT_VALUES));
    _call(VecRestoreArray(uKPetsci_Seq, &u));
    _call(VecRestoreArray(vKPetsci_Seq, &v));
    _call(VecRestoreArray(uKPetsciPrev_Seq, &uprev));
    _call(VecRestoreArray(vKPetsciPrev_Seq, &vprev));
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    _call(VecAssemblyBegin(vKPetsci));
    _call(VecAssemblyEnd(vKPetsci));

    _solutionModified = false;
    _velocitySolutionModified = false;
    _accelerationSolutionModified = false;
#else
    int activeMechDofs = _activeDispNumDofs+_activeStochasticNumDofs;
    double* u, *v, *a;
    _call(VecGetArray(uKPetsci, &u));
    _call(VecGetArray(vKPetsci, &v));
    _call(VecGetArray(aKPetsci, &a));
    
    for (int j = 0; j < activeMechDofs; j++) 
    {
        int ind = ISMapLocal[j];
        u[ind] += dt * v[ind] + dt * dt * (0.5 - _beta) * a[ind];
        v[ind] += dt * (1 - _gamma) * a[ind];
        a[ind] = 0; 
    }
    _call(VecRestoreArray(uKPetsci, &u));
    _call(VecRestoreArray(vKPetsci, &v));
    _call(VecRestoreArray(aKPetsci, &a));
#endif
}


ExplicitNonlinearSystem::ExplicitNonlinearSystem(): DynamicNonlinearSystem(),_gamma(0.5){}
ExplicitNonlinearSystem::~ExplicitNonlinearSystem(){};

void ExplicitNonlinearSystem::allocateSystem(const vector<int>& activeDofMapLocal,
                             const vector<int>& ISactiveDof,
                             int activeDispNumDofs,
                             int activeStochasticNumDofs,
                             int activeExtraNumDofs)
{
    DynamicNonlinearSystem::allocateSystem(activeDofMapLocal,ISactiveDof,
                    activeDispNumDofs,activeStochasticNumDofs,activeExtraNumDofs);
    _call(VecDuplicate(rhsTotalPetsci, &_massLumped));
    _call(VecZeroEntries(_massLumped));
};

void ExplicitNonlinearSystem::zeroLumpedMassVector() 
{
    _call(VecAssemblyBegin(_massLumped));
    _call(VecAssemblyEnd(_massLumped));
    _call(VecZeroEntries(_massLumped));
}

void ExplicitNonlinearSystem::addToLumpedMassVector(int row, double val)
{
    // no lumped mass
    _call(VecSetValue(_massLumped,row,val,ADD_VALUES));
};

bool ExplicitNonlinearSystem::systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon)
{
    _call(VecAssemblyBegin(_massLumped));
    _call(VecAssemblyEnd(_massLumped));
        
    _call(VecAssemblyBegin(rhsTotalPetsci));
    _call(VecAssemblyEnd(rhsTotalPetsci));
    
    // deltaTotalDofPetsci is used for the local acceleration
    // deltaTotalDofPetsci = rhsTotalPetsci/_massLumped
    _call(VecPointwiseDivide(deltaTotalDofPetsci,rhsTotalPetsci,_massLumped));
    //
    double* localData;
    _call(VecGetArray(deltaTotalDofPetsci,&localData));
    int activeMechDofs = _activeDispNumDofs+_activeStochasticNumDofs;
    int totalActiveDofs =  activeMechDofs+ _activeExtraNumDofs;     
#ifdef PARALLEL
    // mechanics 
    // update all to acceleration vector
    _call(VecSetValues(aKPetsci, ISMapLocal.size(), &ISMapLocal[0], &localData[0], INSERT_VALUES));
    _call(VecAssemblyBegin(aKPetsci));
    _call(VecAssemblyEnd(aKPetsci));
    _call(VecScatterBegin(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
    _call(VecScatterEnd(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
    // modify 
    double* ulocal,* vlocal, *alocal, *alocalprev;
    _call(VecGetArray(uKPetsci_Seq, &ulocal));
    _call(VecGetArray(vKPetsci_Seq, &vlocal));
    _call(VecGetArray(aKPetsci_Seq, &alocal));
    _call(VecGetArray(aKPetsciPrev_Seq, &alocalprev));
    
    int localSize;
    _call(VecGetSize(uKPetsci_Seq,&localSize));
    for (int i = 0; i < localSize; i++)
    {
        int row = _activeDofMapLocal[i];
        if ( row > -1 and row < activeMechDofs)
        {
            // mech dof
            vlocal[i] +=  dt*((1.-_gamma)*alocalprev[i] + _gamma*alocal[i]);
        }
        else if (row >=activeMechDofs)
        {
            // extraDof
            ulocal[i] = dt * alocal[i];
            alocal[i] = 0.;
        }
    }
    // propagate to global vector
    _call(VecSetValues(uKPetsci, _scatterIS.size(), &_scatterIS[0], ulocal, INSERT_VALUES));
    _call(VecSetValues(vKPetsci, _scatterIS.size(), &_scatterIS[0], vlocal, INSERT_VALUES));
    _call(VecSetValues(aKPetsci, _scatterIS.size(), &_scatterIS[0], alocal, INSERT_VALUES));
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    _call(VecAssemblyBegin(vKPetsci));
    _call(VecAssemblyEnd(vKPetsci));
    _call(VecAssemblyBegin(aKPetsci));
    _call(VecAssemblyEnd(aKPetsci));
    _call(VecRestoreArray(uKPetsci_Seq, &ulocal));
    _call(VecRestoreArray(vKPetsci_Seq, &vlocal));
    _call(VecRestoreArray(aKPetsci_Seq, &alocal));
    _call(VecRestoreArray(aKPetsciPrev_Seq, &alocalprev));
#else
    double *u, *v, *a, *aprev;
    _call(VecGetArray(uKPetsci, &u));
    _call(VecGetArray(vKPetsci, &v));
    _call(VecGetArray(aKPetsci, &a));
    _call(VecGetArray(aKPetsciPrev, &aprev));
    for (int i = 0; i < activeMechDofs; i++) 
    {
        int ind = ISMapLocal[i];
        a[ind]  = localData[i];
        v[ind] +=  dt*((1.-_gamma)*aprev[ind] + _gamma*a[ind]);
    };
    for (int i=activeMechDofs; i<totalActiveDofs; i++)
    {
        int ind = ISMapLocal[i];
        u[ind] = dt * localData[i];
    };
    _call(VecRestoreArray(uKPetsci, &u));
    _call(VecRestoreArray(vKPetsci, &v));
    _call(VecRestoreArray(aKPetsci, &a));
    _call(VecRestoreArray(aKPetsci, &aprev));
#endif 
    _call(VecRestoreArray(deltaTotalDofPetsci, &localData));
    _velocitySolutionModified = false;
    _accelerationSolutionModified = false;
    _solutionModified = false;
    return true;
};


void ExplicitNonlinearSystem::dynamicPrediction(double timeRun, double dt, const vector<int>& scatter_uK) 
{
#ifdef PARALLEL
    // get local data
    int value = _solutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _solutionModified = (value==1);
    if (_solutionModified)
    {
        _call(VecAssemblyBegin(uKPetsci));
        _call(VecAssemblyEnd(uKPetsci));
        _call(VecScatterBegin(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_uKScatter, uKPetsci, uKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _solutionModified = false;
    }
    
    value = _velocitySolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _velocitySolutionModified = (value==1);
    if (_velocitySolutionModified)
    {
        _call(VecAssemblyBegin(vKPetsci));
        _call(VecAssemblyEnd(vKPetsci));
        _call(VecScatterBegin(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_vKScatter, vKPetsci, vKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _velocitySolutionModified = false;
    }
    
    value = _accelerationSolutionModified ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    _accelerationSolutionModified = (value==1);
    if (_accelerationSolutionModified)
    {
        _call(VecAssemblyBegin(aKPetsci));
        _call(VecAssemblyEnd(aKPetsci));
        _call(VecScatterBegin(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _call(VecScatterEnd(_aKScatter, aKPetsci, aKPetsci_Seq, INSERT_VALUES, SCATTER_FORWARD));
        _accelerationSolutionModified = false;
    }
    // update displacement
    _call(VecAXPY(uKPetsci_Seq,dt,vKPetsci_Seq));
    double fact = dt * dt * 0.5;
    _call(VecAXPY(uKPetsci_Seq,fact,aKPetsci_Seq));
    double *u, *uprev;
    _call(VecGetArray(uKPetsci_Seq, &u));
    _call(VecGetArray(uKPetsciPrev_Seq, &uprev));
    for ( int i=0; i< scatter_uK.size(); i++)
    {
        int rowR = _activeDofMapLocal[i];
        if (rowR <0)
        {
            u[i] = uprev[i];
        };
    };
    _call(VecSetValues(uKPetsci, scatter_uK.size(), &scatter_uK[0], u, INSERT_VALUES));
    _call(VecRestoreArray(uKPetsci_Seq, &u));
    _call(VecRestoreArray(uKPetsciPrev_Seq, &uprev));
    _call(VecAssemblyBegin(uKPetsci));
    _call(VecAssemblyEnd(uKPetsci));
    _solutionModified = false;
#else
    int activeMechDofs = _activeDispNumDofs+_activeStochasticNumDofs;
    double* u, *v, *a;
    _call(VecGetArray(uKPetsci, &u));
    _call(VecGetArray(vKPetsci, &v));
    _call(VecGetArray(aKPetsci, &a));
    
    for (int j = 0; j < activeMechDofs; j++) 
    {
        int ind = ISMapLocal[j];
        u[ind] += dt * v[ind] + 0.5* dt * dt  * a[ind];
    }
    _call(VecRestoreArray(uKPetsci, &u));
    _call(VecRestoreArray(vKPetsci, &v));
    _call(VecRestoreArray(aKPetsci, &a));
#endif
}
