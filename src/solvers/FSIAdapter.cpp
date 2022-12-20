//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file FSIAdapter.cpp
  \brief This file contains all functions related the FSI adapter using preCICE
*/


#include "FSIAdapter.h"
#include "classGPs.h"

#ifdef FSI

/*! \brief constructor
  @param[in] config_file_name path to the configure file
  @param[in] solver_name name of the solver
  @param[in] mesh_name name of the mesh
  @param[in] data_write_name name of data for writing
  @param[in] data_read_name name of data for reading
*/
FSIAdapter::FSIAdapter(const std::string config_file_name,
               const std::string solver_name, 
               const std::string mesh_name,
               const std::string data_write_name,
               const std::string data_read_name):
               _meshName(mesh_name), _dataWriteName(data_write_name), _dataReadName(data_read_name), 
               _vertexSize(0), _dtIterFSI(0.), _timeIterFSI(0.)
{
  int rank = 0;
  int nprcs = 1;
#ifdef PARALLEL 
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nprcs);
#endif
  _precice = new precice::SolverInterface(solver_name,config_file_name,rank,nprcs);
};

FSIAdapter::~FSIAdapter()
{
  delete _precice;
}

/*! \brief initialize function
  @param[in] NeumannBCs vector contains all the Neumann BCs in which the FSI interface can be loaded
*/
double FSIAdapter::initialize(std::vector<classNeumannBCs*> &NeumannBCs)
{
  classFSI *input_FSI = NULL;
  for (int i = 0; i < NeumannBCs.size(); i++) 
  {
    if (NeumannBCs[i]->getType() == "FSI") 
    {
       input_FSI = static_cast<classFSI *>(NeumannBCs[i]);
       break;
    }
  }
  if (input_FSI == NULL)
  {
    ERROR("classFSI must be used in NeumannBCs to make FSI simulation");
    exit(-1);
  }
  
  _vertexSize = input_FSI->getNumberOfBoundaryNodes();
  _boundaryNodes = input_FSI->getListOfBoundaryNodes();
  std::vector<std::vector<double> >& bcCoordinates = input_FSI->getListOfBoundaryCoordinates();

  int dim = _precice->getDimensions();
  std::vector<double> boundary_coords(_vertexSize*dim); // coords of coupling vertices
  for (int j=0; j <_vertexSize; j++)
  {
    std::vector<double>& xyz = bcCoordinates[j];
    if (dim == 2) 
    { //2D
      boundary_coords[dim*j] = xyz[0];
      boundary_coords[dim*j + 1] = xyz[1];
    } else if (dim == 3) 
    { //3D
      boundary_coords[dim*j] = xyz[0];
      boundary_coords[dim*j + 1] = xyz[1];
      boundary_coords[dim*j + 2] = xyz[2];
    }
  };
  
  int mesh_id  = _precice->getMeshID(_meshName);
  _vertexIds.resize(_vertexSize);
  _precice->setMeshVertices(mesh_id,_vertexSize, boundary_coords.data(), _vertexIds.data());
  _forcesPrev.resize(_vertexSize*dim,0.);
  _forces.resize(_vertexSize*dim,0.);
  return _precice->initialize();
};

/*! \brief update the displacement at interface and send to other solver
  @param[in] dt time step
  @param[in] disp vector of all degrees of freedom of unknown field
  @param[in] numNodes number of nodes of the domain
*/
void FSIAdapter::updateDisplacementInterface(double dt, const std::vector<double>& disp, int numNodes)
{
  if (_precice->isWriteDataRequired(dt))
  {
    int dim = _precice->getDimensions();
    std::vector<double> displacements(_vertexSize*dim,0.);
    // displacements : vector containing the value of the displacement on the boundary nodes
    for (int j=0; j <_vertexSize; j++)
    {
      int i = _boundaryNodes[j];
      if (dim == 2) 
      { //2D
        displacements[dim*j] = disp[i];
        displacements[dim*j + 1] = disp[i + numNodes];
      } 
      else if (dim == 3) 
      { //3D
        displacements[dim*j] = disp[i];
        displacements[dim*j + 1] = disp[i + numNodes];
        displacements[dim*j + 2] = disp[i + 2*numNodes];
      }
    }
    int mesh_id  = _precice->getMeshID(_meshName);
    int write_id = _precice->getDataID(_dataWriteName, mesh_id);
    _precice->writeBlockVectorData(write_id, _vertexSize, _vertexIds.data(), displacements.data());
  };
}

/*! \brief read data from the other solver
*/
void FSIAdapter::readBlockVectorData()
{
  if (_precice->isReadDataAvailable()) 
  {
    int mesh_id  = _precice->getMeshID(_meshName);
    int read_id  = _precice->getDataID(_dataReadName, mesh_id);
    _forcesPrev = _forces;
    _precice->readBlockVectorData(read_id, _vertexSize, _vertexIds.data(), _forces.data());
  };
}

/*! \brief update the external force vector with the data from other solver
  @param[in] force external force vector
  @param[in] numNodes number of nodes of the domain
*/
void FSIAdapter::updateForceInterface(std::vector<double>& force, int numNodes, double fact)
{  
  int dim = _precice->getDimensions();
  for (int j=0; j <_vertexSize; j++)
  {
    int i = _boundaryNodes[j];
    if (dim == 2) 
    { //2D
      force[i] += (1.-fact)*_forcesPrev[dim * j] + fact*_forces[dim * j];
      force[i + numNodes] += (1.-fact)*_forcesPrev[dim * j + 1] + fact*_forces[dim * j + 1];
    }
    else if (dim == 3) 
    { //3D
      force[i] += (1.-fact)*_forcesPrev[dim * j] + fact*_forces[dim * j];
      force[i + numNodes] += (1.-fact)*_forcesPrev[dim * j + 1] + fact*_forces[dim * j + 1];
      force[i + 2 * numNodes] += (1.-fact)*_forcesPrev[dim * j + 2] + fact*_forces[dim * j + 2];
    }
  };
};

/*! \brief advance function in preCICE
  @param[in] dt time step
*/
double FSIAdapter::advance(double dt)
{
  return _precice->advance(dt);
};

bool FSIAdapter::isCouplingOngoing()
{
  return _precice->isCouplingOngoing();
}

/*! \brief finalize when finishing coupling
*/
void FSIAdapter::finalize()
{
    _precice->finalize();
};

const std::string& FSIAdapter::actionReadIterationCheckpoint()
{
  return precice::constants::actionReadIterationCheckpoint(); 
}
const std::string& FSIAdapter::actionWriteIterationCheckpoint()
{
  return precice::constants::actionWriteIterationCheckpoint();
}

bool FSIAdapter::isActionRequired(const std::string& action)
{
  return _precice->isActionRequired(action);
}
void FSIAdapter::markActionFulfilled(const std::string& action)
{
  _precice->markActionFulfilled(action);
}

/*! \brief make a checkpoint of the solver
  @param[in] timeRun time
  @param[in] dt time step
*/
void FSIAdapter::makeACheckPoint(double timeRun, double dt)
{
  _timeIterFSI = timeRun;
  _dtIterFSI = dt;
};

/*! \brief load a checkpoint of the solver
  @param[out] timeRun time
  @param[out] dt time step
*/
void FSIAdapter::restoreACheckPoint(double& timeRun, double& dt)
{
  timeRun = _timeIterFSI;
  dt = _dtIterFSI;
};

#endif //FSI