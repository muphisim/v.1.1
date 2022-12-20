//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _solverImplicit_H_
#define _solverImplicit_H_

#include "configuration.h"
#include "solvers.h"
#include "nonlinearSystem.h"
#include "maths.h"

/*! \brief This file contains all functions related the implicit solver (quasi-static and Newmark Algorithm). Go to the time stepping documentation of MuPhiSim for further information*/
class classSolverImplicitBase : public solvers {
 
    public:
        classSolverImplicitBase(vector<classNodes *> &nod, int ndim, int numDof, vector<classDirichlet *> &DNod,
                            double scaleFactor, classOutputs *&outputs, double t0, double tf, string solverType,
                            const vector<int> &activeDofMapLocal, const vector<int>& nActiveDof, const vector<int> &ISactiveDof,
                            classPrintForces *&forceOutputs, bool flagActivation,
                            string &petsc_solver, string &petsc_precon, bool flagFSI);
        virtual ~classSolverImplicitBase() {};
        virtual nonlinearSystem* createSystem() const =0;
        virtual void setLastSolutionsToSystem(nonlinearSystem* nlsys, 
                                        const vector<int>& scatter_uK,
                                        const vector<int>& scatter_vvK,
                                        const vector<double>& uK,
                                        const vector<double>& vK,
                                        const vector<double>& aK,
                                        const vector<double>& vvK) const =0;
        virtual void getSolutionsFormSystem(nonlinearSystem* nlsys, 
                                        int numMechDofs,  int numExtraDof,
                                        vector<double>& uK,
                                        vector<double>& vK,
                                        vector<double>& aK,
                                        vector<double>& vvK) const =0;
        virtual void computeRightHandSide(nonlinearSystem* nlsys, vector<double>&inerForce, vector<double>&intForce, vector<double>&extForce,
                            int totalNumberDofsPerNode, const vector<double> &allField, const vector<double>& allAcc, 
                            vector<classGPs *> &GPs, vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, double timeRun,
                            double dt, std::map<pair<int,int>,double>& stiffBC, bool tangentEstimation, double FSIstepRatio) =0;
        
        virtual void applyDisplacementBC(nonlinearSystem* nlsys, int numMechDofs, double timeRun, double dt, 
                                         const vector<int> &nodLocal, int nNodGlobal, int nNodLocal) =0;
        
        virtual void computeStiffness(double timeRun, double dt, vector<double>& uK, 
                                        const vector<classGPs *>& GPs, nonlinearSystem* nlsys) const;
        virtual void computeMassMatrix(const vector<classGPs *>& GPs, classMatrix& M) const;
        virtual void computeInternalForce(const vector<classGPs *>& GPs, vector<double>& fint, double timeRun, double dt) const;
        virtual void calculate(vector<classNeumannBCs *> &NeumannBCs, vector<classGPs *> &GPs, int sim,
                                vector<classNodes *> &nod, vector<classElements *> &elem, const vector<int> &nodFEM,
                                const vector<int> &nodMM, const vector<int> &nodShare, const vector<int> &elemFEM, const vector<int> &elemMM,
                                vector<classElements *> &geoShape, vector<constitutiveModels *> &consMod,
                                const vector<int> &nodLocal, int nNodGlobal, int nNodLocal);

       
};


/*! \brief This file contains all functions related the static implicit solver. Go to the time stepping documentation of MuPhiSim for further information*/
class classSolverImplicitStatic : public classSolverImplicitBase {

    public:
        classSolverImplicitStatic(vector<classNodes *> &nod, int ndim, int numDof, vector<classDirichlet *> &DNod,
                              double scaleFactor, classOutputs *&outputs, double t0, double tf, string solverType,
                              const vector<int> &activeDofMapLocal, const vector<int>& nActiveDof, const vector<int> &ISactiveDof,
                              classPrintForces *&forceOutputs, bool flagActivation, string &petsc_solver, string &petsc_precon,
                              bool flagFSI);
        virtual ~classSolverImplicitStatic(){}
        
        virtual nonlinearSystem* createSystem()  const {return new ImplicitStaticNonlinearSystem();};
        virtual void setLastSolutionsToSystem(nonlinearSystem* nlsys, 
                                        const vector<int>& scatter_uK,
                                        const vector<int>& scatter_vvK,
                                        const vector<double>& uK,
                                        const vector<double>& vK,
                                        const vector<double>& aK,
                                        const vector<double>& vvK) const;
        virtual void getSolutionsFormSystem(nonlinearSystem* nlsys, 
                                        int numMechDofs,  int numExtraDof,
                                        vector<double>& uK,
                                        vector<double>& vK,
                                        vector<double>& aK,
                                        vector<double>& vvK) const;
        virtual void computeRightHandSide(nonlinearSystem* nlsys, vector<double>&inerForce, vector<double>&intForce, vector<double>&extForce,
                            int totalNumberDofsPerNode, const vector<double> &allField, const vector<double>& allAcc, 
                            vector<classGPs *> &GPs, vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, double timeRun,
                            double dt, map<pair<int,int>,double>& stiffBC, bool tangentEstimation, double FSIstepRatio);
        
        virtual void applyDisplacementBC(nonlinearSystem* nlsys, int numMechDofs, double timeRun, double dt, 
                                         const vector<int> &nodLocal, int nNodGlobal, int nNodLocal);
                                         
        inline virtual void initValues(const vector<double> &u, const vector<double> &v, const vector<double> &a, const vector<double> &vv, const vector<double> &extF) {
            _u = u;
            _v.resize(_u.size());
            setAll(_v,0);
            _a.resize(_u.size());
            setAll(_a,0);
            _vv = vv;
            _extF = extF;
        };
};

/*! \brief This file contains all functions related the dynamic implicit solver (Newmark Algorithm). Go to the time stepping documentation of MuPhiSim for further information*/
class classSolverImplicit : public classSolverImplicitBase {

    public:
        classSolverImplicit(vector<classNodes *> &nod, int ndim, int numDof, vector<classDirichlet *> &DNod,
                              double scaleFactor, classOutputs *&outputs, double t0, double tf, string solverType,
                              const vector<int> &activeDofMapLocal, const vector<int>& nActiveDof, const vector<int> &ISactiveDof,
                              classPrintForces *&forceOutputs, bool flagActivation, string &petsc_solver, string &petsc_precon,
                              bool flagFSI);
        virtual ~classSolverImplicit(){}
        virtual nonlinearSystem* createSystem()  const {return new ImplicitNonlinearSystem();};
        virtual void setLastSolutionsToSystem(nonlinearSystem* nlsys, 
                                        const vector<int>& scatter_uK,
                                        const vector<int>& scatter_vvK,
                                        const vector<double>& uK,
                                        const vector<double>& vK,
                                        const vector<double>& aK,
                                        const vector<double>& vvK) const;
        virtual void getSolutionsFormSystem(nonlinearSystem* nlsys, 
                                        int numMechDofs,  int numExtraDof,
                                        vector<double>& uK,
                                        vector<double>& vK,
                                        vector<double>& aK,
                                        vector<double>& vvK) const;
        virtual void computeRightHandSide(nonlinearSystem* nlsys, vector<double>&inerForce, vector<double>&intForce, vector<double>&extForce,
                            int totalNumberDofsPerNode, const vector<double> &allField, const vector<double>& allAcc, 
                            vector<classGPs *> &GPs, vector<classNodes *> &nod, vector<classElements *> &elem,
                            vector<classNeumannBCs *> &NeumannBCs, int &ndim, double timeRun,
                            double dt, map<pair<int,int>,double>& stiffBC, bool tangentEstimation, double FSIstepRatio);
        
        virtual void applyDisplacementBC(nonlinearSystem* nlsys, int numMechDofs, double timeRun, double dt, 
                                         const vector<int> &nodLocal, int nNodGlobal, int nNodLocal);
};

#endif
