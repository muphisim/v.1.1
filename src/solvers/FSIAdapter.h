//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file FSIAdapter.h
  \brief This file contains all functions related the FSI adapter using preCICE
*/


#ifndef _FSIAdapter_H_
#define _FSIAdapter_H_

#ifdef FSI
#include "SolverInterface.hpp"
#include <string>
#include <vector>
#include "classNeumannBCs.h"
#include "classNodes.h"

class FSIAdapter
{
  protected:
    precice::SolverInterface* _precice;
    // data
    std::string _meshName, _dataWriteName, _dataReadName;
    std::vector<int> _boundaryNodes;
    std::vector<int> _vertexIds; // list of boundary nodes
    int _vertexSize; // number of boundary nodes
    std::vector<double> _forces;
    std::vector<double> _forcesPrev;

    // data checkpoint of the solver
    double _dtIterFSI, _timeIterFSI;

  
  public:
    FSIAdapter(const std::string config_file_name="../precice-config.xml",
               const std::string solver_name = "Solid", 
               const std::string mesh_name = "Solid-Mesh",
               const std::string data_write_name = "Displacement",
               const std::string data_read_name ="Force");
    ~FSIAdapter();
    double initialize(std::vector<classNeumannBCs*> &NeumannBCs);
    void updateDisplacementInterface(double dt, const std::vector<double>& disp, int numNodes);
    void readBlockVectorData();
    void updateForceInterface(std::vector<double>& force, int numNodes, double fact);
    double advance(double dt);
    bool isCouplingOngoing();
    void finalize();
    
    const std::string& actionReadIterationCheckpoint();
    const std::string& actionWriteIterationCheckpoint();
    
    bool isActionRequired(const std::string& action);
    void markActionFulfilled(const std::string& action);
    
                      
    void makeACheckPoint(double timeRun, double dt);
    void restoreACheckPoint(double& timeRun, double& dt);
};

#endif // FSI

#endif //_FSIAdapter_H_
