//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _solverExplicit_H_
#define _solverExplicit_H_

#include "configuration.h"
#include "solvers.h"
#include "maths.h"

/*! \brief This file contains all functions related the explicit solver (Newmark Algorithm). Go to the time stepping documentation of MuPhiSim for further information*/
class classSolverExplicit : public solvers {
private:
    int _flagAV;
    double _c1;
    double _cl;

public:
    classSolverExplicit(vector<classNodes *> &nod, int ndim, int numDof, vector<classDirichlet *> &DNod,
                        double scaleFactor, classOutputs *&outputs, double t0, double tf, string solverType,
                        const vector<int> &activeDofMapLocal, const vector<int>& nActiveDof, const vector<int> &ISactiveDof, int flag_AV,
                        const vector<double>& para_AV, classPrintForces *&forceOutputs, bool flagActivation,
                        string &petsc_solver, string &petsc_precon, bool flagFSI);

    virtual void calculate(vector<classNeumannBCs *> &NeumannBCs, vector<classGPs *> &GPs, int sim,
                           vector<classNodes *> &nod, vector<classElements *> &elem, const vector<int> &nodFEM,
                           const vector<int> &nodMM, const vector<int> &nodShare, const vector<int> &elemFEM, const vector<int> &elemMM,
                           vector<classElements *> &geoShape, vector<constitutiveModels *> &consMod,
                           const vector<int> &nodLocal, int nNodGlobal, int nNodLocal);
    
    virtual void computeInternalForce(const vector<classGPs *>& GPs, vector<double>& fint, double timeRun, double dt) const;
    // for explicit scheme
    virtual void computeLumpedMassVector(const vector<classGPs *>& GPs, nonlinearSystem& nlsys) const;
    
    void artificialViscosity(double elementSize, const vector<double> &F1, const vector<double> &F0, vector<double> &Piola_vis,
                             double rho, double dt) const;

    void bulkViscosity(double elementSize1, const vector<double>& F1, const vector<double>& F0, vector<double> &Piola_vis, double rho,
                       double dt) const;
     
    virtual void computeRightHandSide(nonlinearSystem* nlsys, vector<double>&inerForce, vector<double>&intForce, vector<double>&extForce,
                            int totalNumberDofsPerNode, const vector<double> &allField, const vector<double>& allAcc, 
                            vector<classGPs *> &GPs, vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, double timeRun,
                            double dt, double FSIstepRatio);
    
    virtual void setLastSolutionsToSystem(nonlinearSystem* nlsys, 
                                        const vector<int>& scatter_uK,
                                        const vector<int>& scatter_vvK,
                                        const vector<double>& uK,
                                        const vector<double>& vK,
                                        const vector<double>& aK,
                                        const vector<double>& vvK) const;
    virtual void applyDisplacementBC(nonlinearSystem* nlsys, int numMechDofs, double timeRun, double dt, 
                                         const vector<int> &nodLocal, int nNodGlobal, int nNodLocal);
    virtual void getSolutionsFormSystem(nonlinearSystem* nlsys, 
                                        int numMechDofs,  int numExtraDof,
                                        vector<double>& uK,
                                        vector<double>& vK,
                                        vector<double>& aK,
                                        vector<double>& vvK) const;
    ~classSolverExplicit() {};
};


#endif

