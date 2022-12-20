//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file nonlinearSystem.h
  \brief This file contains all functions related to the system
*/


#ifndef _nonlinearSystem_H_
#define _nonlinearSystem_H_

#include "configuration.h"
#include "petsc.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"


class nonlinearSystem
{
    public:
        virtual ~nonlinearSystem(){}
        virtual bool isAllocated() const = 0;
        virtual void allocateSolution(int totalSize, const vector<int>& scatterIS) = 0;
        virtual void clearSolution() = 0;
        virtual void allocateSystem(const vector<int>& activeDofMapLocal,
                                     const vector<int>& ISactiveDof,
                                     int activeDispNumDofs,
                                     int activeStochasticNumDofs,
                                     int activeExtraNumDofs) = 0;
        virtual void clearSystem() = 0;
        virtual void zeroRHS() = 0;
        virtual void addToRHS(int row, double val) = 0;
        virtual void getNormInfRHS(vector<double>& norm)= 0; // return norm inf of each field
        virtual void zeroMatrix() = 0;
        virtual void addToMatrix(int row, int col, double val) = 0;
        virtual void zeroLumpedMassVector() = 0;
        virtual void addToLumpedMassVector(int row, double val) = 0;
        virtual bool systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon) = 0;
        virtual void setToSolution(int row, double val) = 0;
        virtual void setToSolution(const vector<int>& row, const vector<double>& val) = 0;
        virtual void addToSolution(int row, double val) = 0;
        virtual void setToVelocitySolution(int row, double val) = 0;
        virtual void setToVelocitySolution(const vector<int>& row, const vector<double>& val) = 0;
        virtual void addToVelocitySolution(int row, double val) = 0;
        virtual void setToAccelerationSolution(int row, double val) = 0;
        virtual void setToAccelerationSolution(const vector<int>& row, const vector<double>& val) = 0;
        virtual void addToAccelerationSolution(int row, double val) = 0;
        virtual void getFromSolution(const vector<int>& rows, vector<double>& uK) = 0;
        virtual void getFromSolutionByRange(int startIndex, int endIndex, vector<double>& uK) = 0;
        virtual void getFromVelocitySolution(const vector<int>& rows, vector<double>& uK) = 0;
        virtual void getFromVelocitySolutionByRange(int startIndex, int endIndex, vector<double>& uK) = 0;
        virtual void getFromAccelerationSolution(const vector<int>& rows, vector<double>& uK) = 0;
        virtual void getFromAccelerationSolutionByRange(int startIndex, int endIndex, vector<double>& uK) = 0;
        virtual void nextStep() = 0;
        virtual void resetToPreviousSolution() = 0;
        virtual classMatrix* getMatrix() = 0;
        virtual classMatrix* getMassMatrix() = 0;
        virtual bool withStiffnessMatrix() const = 0;
        virtual bool withMassMatrix() const = 0;
        virtual void dynamicPrediction(double timeRun, double dt, const vector<int>& scatter_uK) = 0;
        
        virtual void makeACheckPoint() = 0;
        virtual void restoreACheckPoint() = 0;
        
};

class ImplicitStaticNonlinearSystem : public nonlinearSystem
{
    protected:
        void _call(PetscErrorCode ierr) const {CHKERRV(ierr);}
    
    protected:
        // for solution
        bool _solutionAllocated;
        Vec uKPetsci, uKPetsciPrev, uKPetsciTmp;         //Displacement Vector (all DOF)  
        bool _solutionModified;
        bool _tmpStateAllocated;
        
#ifdef PARALLEL
        vector<int> _scatterIS;
        VecScatter _uKScatter;
        Vec uKPetsci_Seq, uKPetsciTmp_Seq;  //Sequental vector to gather the displacements of my processor
#endif
      
        // for system
        bool _systemAllocated;
        Vec rhsTotalPetsci;        //Right Hand Side Vector (fext-fint) (only the active DOF)
        Vec deltaTotalDofPetsci;   // local increment
        classMatrix* _stiffnessMatrix; // stiffness matrix
        //
        Vec rhsDispPetsci;
        Vec rhsStochasticPetsci;
        Vec rhsExtraPetsci;
        
        vector<int> _activeDofMapLocal;
        int _activeDispNumDofs; // 
        int _activeStochasticNumDofs; //
        int _activeExtraNumDofs;  //
        
        // Is local to all Dofs - local proc to uKPetsci
        vector<int> ISMapLocal; 
        bool _matrixModified;
                
        KSP _ksp;
        bool _kspIsAllocated;
        
    public:
        ImplicitStaticNonlinearSystem();
        virtual ~ImplicitStaticNonlinearSystem();
        virtual bool isAllocated() const {return _systemAllocated and _solutionAllocated;};
        virtual void allocateSolution(int totalDofs, const vector<int>& scatterIS);
        virtual void clearSolution();
        virtual void allocateSystem(const vector<int>& activeDofMapLocal,
                                     const vector<int>& ISactiveDof,
                                     int activeDispNumDofs,
                                     int activeStochasticNumDofs,
                                     int activeExtraNumDofs);
        
        virtual void clearSystem(); 
        virtual void zeroRHS();
        virtual void addToRHS(int row, double val);
        virtual void getNormInfRHS(vector<double>& norm);
        virtual void zeroMatrix();
        virtual void addToMatrix(int row, int col, double val);
        virtual void zeroLumpedMassVector() {} // no lumped mass
        virtual void addToLumpedMassVector(int row, double val){}// no lumped mass
        virtual bool systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon);
        virtual void setToSolution(int row, double val);
        virtual void setToSolution(const vector<int>& row, const vector<double>& val);
        virtual void addToSolution(int row, double val);
        
        virtual void setToVelocitySolution(int row, double val){};
        virtual void setToVelocitySolution(const vector<int>& row, const vector<double>& val){};
        virtual void addToVelocitySolution(int row, double val){};
        
        virtual void setToAccelerationSolution(int row, double val) {};
        virtual void setToAccelerationSolution(const vector<int>& row, const vector<double>& val) {};
        virtual void addToAccelerationSolution(int row, double val) {};
        
        virtual void getFromSolution(const vector<int>& rows, vector<double>& uK);
        virtual void getFromSolutionByRange(int startIndex, int endIndex, vector<double>& uK);
        virtual void getFromVelocitySolution(const vector<int>& rows, vector<double>& uK){};
        virtual void getFromVelocitySolutionByRange(int startIndex, int endIndex, vector<double>& uK){};
        virtual void getFromAccelerationSolution(const vector<int>& rows, vector<double>& uK){};
        virtual void getFromAccelerationSolutionByRange(int startIndex, int endIndex, vector<double>& uK){};
        
        virtual void nextStep();
        virtual void resetToPreviousSolution();
        
        virtual classMatrix* getMatrix() {return _stiffnessMatrix;};
        virtual classMatrix* getMassMatrix() {return NULL;};
        virtual bool withStiffnessMatrix() const {return true;};
        virtual bool withMassMatrix() const {return false;};
        virtual void dynamicPrediction(double timeRun, double dt, const vector<int>& scatter_uK){};
        virtual void makeACheckPoint();
        virtual void restoreACheckPoint();
        
};

class DynamicNonlinearSystem : public ImplicitStaticNonlinearSystem
{
    protected: 
        Vec vKPetsci, vKPetsciPrev, vKPetsciTmp;         //velocity Vector (all DOF)  
        Vec aKPetsci, aKPetsciPrev, aKPetsciTmp;         //velocity Vector (all DOF)
#ifdef PARALLEL
        VecScatter _vKScatter, _aKScatter;
        Vec vKPetsci_Seq, aKPetsci_Seq, vKPetsciTmp_Seq, aKPetsciTmp_Seq;  //Sequental vector to gather the displacements of my processor
        Vec uKPetsciPrev_Seq, vKPetsciPrev_Seq, aKPetsciPrev_Seq;
#endif
        classMatrix* _massMatrix;
        bool _velocitySolutionModified;
        bool _accelerationSolutionModified;
            
    public:
        DynamicNonlinearSystem();
        virtual ~DynamicNonlinearSystem();
        
        virtual void allocateSolution(int totalDofs, const vector<int>& scatterIS);
        virtual void clearSolution();
        
         virtual void allocateSystem(const vector<int>& activeDofMapLocal,
                                     const vector<int>& ISactiveDof,
                                     int activeDispNumDofs,
                                     int activeStochasticNumDofs,
                                     int activeExtraNumDofs);
        
        virtual void clearSystem();
        virtual void setToVelocitySolution(int row, double val);
        virtual void setToVelocitySolution(const vector<int>& row, const vector<double>& val);
        virtual void addToVelocitySolution(int row, double val);
        
        virtual void setToAccelerationSolution(int row, double val);
        virtual void setToAccelerationSolution(const vector<int>& row, const vector<double>& val);
        virtual void addToAccelerationSolution(int row, double val);
        virtual void zeroLumpedMassVector()  =0; // no lumped mass
        virtual void addToLumpedMassVector(int row, double val) =0; // virtual function
        
        virtual bool systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon)=0; // virtual function
        
        virtual void getFromVelocitySolution(const vector<int>& rows, vector<double>& uK);
        virtual void getFromVelocitySolutionByRange(int startIndex, int endIndex, vector<double>& uK);
        
        virtual void getFromAccelerationSolution(const vector<int>& rows, vector<double>& uK);
        virtual void getFromAccelerationSolutionByRange(int startIndex, int endIndex, vector<double>& uK);
        
        virtual void nextStep();
        virtual void resetToPreviousSolution();
        virtual classMatrix* getMassMatrix() {return _massMatrix;};
        virtual bool withStiffnessMatrix() const {return true;};
        virtual bool withMassMatrix() const {return true;};
        virtual void dynamicPrediction(double timeRun, double dt, const vector<int>& scatter_uK) =0; // virtual function
        
        virtual void makeACheckPoint();
        virtual void restoreACheckPoint();
};


class ImplicitNonlinearSystem : public DynamicNonlinearSystem
{
    protected: 
        double _gamma, _beta; 
            
    public:
        ImplicitNonlinearSystem();
        virtual ~ImplicitNonlinearSystem();
        
        virtual void zeroLumpedMassVector()  {}; // no lumped mass
        virtual void addToLumpedMassVector(int row, double val){};
        
        virtual bool systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon);
        virtual bool withStiffnessMatrix() const {return true;};
        virtual bool withMassMatrix() const {return true;};
        virtual void dynamicPrediction(double timeRun, double dt, const vector<int>& scatter_uK);
};


class ExplicitNonlinearSystem : public DynamicNonlinearSystem
{
    public:
        Vec _massLumped;
        double _gamma;
        
    public:
        ExplicitNonlinearSystem();
        virtual ~ExplicitNonlinearSystem();
        
        virtual void allocateSystem(const vector<int>& activeDofMapLocal,
                                     const vector<int>& ISactiveDof,
                                     int activeDispNumDofs,
                                     int activeStochasticNumDofs,
                                     int activeExtraNumDofs);
        
        virtual void zeroLumpedMassVector() ;
        virtual void addToLumpedMassVector(int row, double val);        
        virtual bool systemSolve(double timeRun, double dt, string petsc_solver, string petsc_precon);
        
        virtual classMatrix* getMassMatrix() {return NULL;};
        virtual classMatrix* getMatrix() {return NULL;}
        virtual bool withStiffnessMatrix() const {return false;};
        virtual bool withMassMatrix() const {return false;};
        
        virtual void dynamicPrediction(double timeRun, double dt, const vector<int>& scatter_uK);
    
};

#endif //_nonlinearSystem_H_