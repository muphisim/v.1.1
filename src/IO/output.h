//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _output_H_
#define _output_H_

#include "configuration.h"
#include "classGPs.h"
#include "classNodes.h"
#include "classElements.h"

/*! \brief  Outputs class to define the outputs variables*/
class classOutputs {

protected:

    int _numFiles; /*!< \brief Number of files to be printed (approximation) */
    vector<int> _outVars;

public:

    classOutputs();

    inline void setNumFiles(int data) {
        _numFiles = data;
    };

    inline void setOutVars(const vector<int>& data) {
        _outVars = data;
    };

    inline int getNumFiles() {
        return _numFiles;
    };

    inline const vector<int>& getOutVars() {
        return _outVars;
    };

    ~classOutputs() {
        _outVars.clear();
    };
};

/*! \brief  classPrintForces to define the nodes and directions of forces to be printed*/
class classPrintForces {

protected:

    vector<int> _nodesToWrite;
    vector<int> _getDirectionsToWrite;
    double _period;
    int _iterationsCompleted;
    bool _firstTimeForcesOutput;

public:

    classPrintForces();

    inline void setDirectionsToWrite(vector<int> data) {
        this->_getDirectionsToWrite = data;
    };

    inline void setNodesToWrite(const vector<int>& data) {
        this->_nodesToWrite = data;
    };

    inline void setPeriod(double data) {
        this->_period = data;
    };

    inline const vector<int>& getDirectionsToWrite() const {
        return this->_getDirectionsToWrite;
    };

    inline const vector<int>& getNodesToWrite() const {
        return this->_nodesToWrite;
    };

    inline double getPeriod() const {
        return this->_period;
    };

    inline int getIterationsCompleted() const { return this->_iterationsCompleted; };

    inline void setIterationCompleted() { this->_iterationsCompleted++; };

    inline void firstTimeComplete() { this->_firstTimeForcesOutput = false; };

    inline bool isFirstTime() const { return this->_firstTimeForcesOutput; };

    ~classPrintForces() {
        _getDirectionsToWrite.clear();
        _nodesToWrite.clear();
    };
};

extern void
outputsManagement(vector<classNodes *> &nodes, vector<classElements *> &elements, int &ndim, const vector<int> &nodesFEM,
                  const vector<int> &nodesMM, const vector<int> &nodesShare, const vector<int> &elementsFEM, const vector<int> &elementsMM,
                  const vector<double> &nonRealDisp, const vector<double> &vK, const vector<double> &aK, int sim, int numDof,
                  vector<classGPs *>& GPs, classOutputs *outputs, double dt, double timeRun, int rank);

extern void outputsManagement(vector<classNodes *> &nod, vector<classElements *> &elem, int &ndim, const vector<int> &nodFEM,
                              const vector<int> &nodMM, const vector<int> &nodShare, const vector<int> &elemFEM, const vector<int> &elemMM,
                              const vector<double> &uK, const vector<double> &vK, const vector<double> &aK, const vector<double> &vv,
                              const vector<double> &Fv, int sim, int nDof, vector<classGPs *>& GPs, classOutputs *outputs,
                              double dt, double timeRun, int rank);

extern void writeExternalForces(double timeRun, const vector<int> &outputNodes, const vector<double> &extForces, int ndim,
                                const vector<int> &direction, bool firsttime, double t0, double tf);

extern void printForcesOut(const vector<double> &intForce, int totalDofs, double dt, const vector<double> &inertiaForceFull,
                           int type, classPrintForces *&forceOutputs, double _t0, double _tf, double timeRun, int _ndim,
                           int rank, int _numDof, int nNodGlobal, int nNodLocal, const vector<int> &nodLocal, int nprcs);

extern void realDisplacementManagement(vector<classNodes *> &nod, int &ndim, const vector<int> &nodMM, const vector<int> &nodFEM,
                                       const vector<int> &nodShare, const vector<double> &ficDisp, int sim, int nNod, int nNMM,
                                       int nNFEM, vector<double> &U);

extern void realDisplacements(const vector<double> &uKK, vector<double> &realU, int numNodes, vector<classNodes *> &nodes,
                              const vector<int> &nodesMM, int ndim);

#endif
