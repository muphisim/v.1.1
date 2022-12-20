//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _solvers_H_
#define _solvers_H_

#include "configuration.h"
#include <petscsnes.h>
#include <petscksp.h>
#include "geometry.h"
#include "timeStepping.h"
#include "nonlinearSystem.h"
#include "FSIAdapter.h"
#include "errorEstimation.h"


/*! \brief  A function to get  the CPU time of day in seconds.*/
double getCPUTime();
/*! \brief  A function to get  the current time of day in seconds.*/
double getTimeOfDay();

/*! \brief  This is the wrapper for all solvers. This base class implements general and virtual functions. ALL VIRTUAL FUNCTIONS MUST BE IMPLEMENTE BY THE USER WHEN ADDING NEW SOLVERS*/
class solvers {
public:
    enum StiffnessMatrixType{None=0,Current=1, Initial=2, Iteration, Interval, Beginning};
protected:

    int _numNodes; /*!< \brief Number of nodes in the domain*/
    int _ndim; /*!< \brief Dimension of the problem*/
    int _numDof; /*!< \brief Total number of degrees of freedom*/
    double _scaleFactor;
    int _numSteps; /*!< \brief useful in implicit scheme, if _numSteps=0, _scaleFactor is used to determin time step*/
    double _absTol, _relTol; /*!< \brief tolerance in implict scheme or quasi-static scheme*/
    int _maxIterations; /*!< \brief maximal number of iterations in each time step*/
    int _numIterationsIncreaseTimeStep; /*!< \brief when the time step is reduced, it is increase if the number of iterations larger than this value*/
    AdaptiveTimeStep* _adaptiveTimeStep;
    double _maxTimeStep;   /*!< \brief maximal time step in simulation*/
    classOutputs *_outputs;
    classPrintForces *_forceOutputs;
    double _tf, _t0;
    vector<classDirichlet *> _DNod;
    string _solverType;
    vector<double> _u, _v, _a, _vv, _extF;
    vector<int> _ISactiveDof;
    vector<int> _nActiveDof; /*!< \brief Number of active degrees of freedom*/ //Only for the mechanical part
    vector<int> _activeDofMapLocalInitial;
    vector<int> _activeDofMapLocal;
    double _elementSize;
    double _c0;
    bool _flagActivation;
    bool _flagFSI;
    string _petsc_solver;
    string _petsc_precon;
    StiffnessMatrixType _stiffnessMatrixType;
    int _stiffnessMatrixFixedIteration; // the stiffness matrix does not change afther this value, it is applicable only _stiffnessMatrixType = Step    
    TimeStepping _timeSteppingPlan;
    int _comRank, _comSize;
    int _maxFailedSteps;
    bool _tangentByPerturbation;
    double _perturbationTol;
    trueDisplacementField* _trueDispField;
    
#ifdef FSI
    FSIAdapter* _adapter;
#endif 
    
public:

    solvers(vector<classNodes *> &nod, int ndim, int numDof, vector<classDirichlet *> &DNod, double scaleFactor,
            classOutputs *&outputs, double t0, double tf, string solverType,const vector<int> &activeDofMapLocal,
            const vector<int>& nActiveDof, const vector<int> &ISactiveDof, classPrintForces *&forceOutputs, bool flagActivation,
            string &petsc_solver, string &petsc_precon, bool flagFSI);
    virtual ~solvers();
    
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
  @param[in] timeRunGlobal Global time of the simulation (important for multisolvers simulations)
*/
    virtual void calculate(vector<classNeumannBCs *> &NeumannBCs, vector<classGPs *> &GPs, int sim, vector<classNodes *> &nod,
              vector<classElements *> &elem, const vector<int> &nodFEM, const vector<int> &nodMM, const vector<int> &nodShare,
              const vector<int> &elemFEM, const vector<int> &elemMM, vector<classElements *> &geoShape,
              vector<constitutiveModels *> &consMod, const vector<int> &nodLocal, int nNodGlobal, int nNodLocal) = 0;
              
              
    void setTangentByPerturbation(bool fl, double tol);
    void setNumberOfSteps(int nstep);
    void setAbsoluteTolerance(double absTol);
    void setRelativeTolerance(double relTol);
    void setMaximalNumberOfIterations(int maxIt);
    void setMaximalTimeStep(double dtMax);
    void setAllowableFailedSteps(int nstep);
    void setTimeSteppingPlan(const TimeStepping& plan);
    void setStiffnessMatrixType(solvers::StiffnessMatrixType type, int iteration=-1);
    void setNumIterationsIncreaseTimeStep(int nIter);
    void allocateTrueDisplacementFieldForErrorAnalysis(const char what[], const std::vector<double>& data);
    void setAdaptiveTimeStep(const char what[], const vector<double>& data);
    void updateActivation(vector<classGPs *> &GPs, double timeRun, vector<int> &_activeDofMapLocal,
                          vector<int> &_nActiveDof, vector<int> &_ISactiveDof, int totalDofs, int numExtraDof,
                          const vector<int>& nodLocal, int nNodGlobal, int nNodLocal, int nprcs, int rank);
                          
    void nextStep(vector<classGPs *> &GPs, nonlinearSystem* lsys, bool resetInitialGP) const;
    void resetToPreviousStep(vector<classGPs *> &GPs, nonlinearSystem* lsys) const;
    void makeACheckPoint(vector<classGPs *> &GPs, nonlinearSystem* lsys) const;
    void restoreACheckPoint(vector<classGPs *> &GPs, nonlinearSystem* lsys) const;

    double calculateCriticalTimeStep(vector<classNodes *> &nod, vector<classElements *> &elem);

    double calculateCriticalTimeStep(vector<classNodes *> &nod, vector<classElements *> &elem, const vector<double> &uK);

    void preAllocateMatrix(const vector<classGPs *>& GPs, classMatrix &KK_total) const;
    void computeKinematicVariables(vector<classGPs *>& GPs, const vector<double> &uKK) const;
    void computeState(vector<classGPs *>& GPs, const vector<double> &uKK, double dt, double timeRun, bool withStiffness) const;
    void computeInertialForce(const vector<classGPs *>& GPs, const vector<double> &aK, vector<double>& finer) const;
    
    
    inline void setNumNodes(int numNodes) {
        _numNodes = numNodes;
    };

    inline int getNumNodes() const {
        return _numNodes;
    };

    inline int getNumMulti() const{
        return _nActiveDof.size() - 1;
    };

    inline string getType() const{
        return _solverType;
    };

    inline double getSf() const{
        return _scaleFactor;
    };

    inline classOutputs *getOutputs() {
        return _outputs;
    };

    inline double getTf() const{
        return _tf;
    };

    inline double getT0() const{
        return _t0;
    };

    inline vector<classDirichlet *>& getDirBCs() {
        return _DNod;
    };

    inline void
    getFinalValues(vector<double> &u, vector<double> &v, vector<double> &a, vector<double> &vv, vector<double> &extF) const {
        u = _u;
        v = _v;
        a = _a;
        vv = _vv;
        extF = _extF;
    };

    inline virtual void
    initValues(const vector<double> &u, const vector<double> &v, const vector<double> &a, const vector<double> &vv, const vector<double> &extF) {
        _u = u;
        _v = v;
        _a = a;
        _vv = vv;
        _extF = extF;
    };

};

#include "solversList.h"

#endif
