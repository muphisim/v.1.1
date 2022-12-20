//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classGPs.cpp
  \brief This file contains the definition of the functions used of the classGPs
*/

#include "classGPs.h"
#include "constitutiveManager.h"


/*! Constructor for classGPs
  @param[in] me Label for the GP
  @param[in] xyz Position of the node
  @param[in] weight Weight for the integration
  @param[in] Jacobian Jacobian for the integration
  @param[in] ndim Dimension of the domain 
*/
classGPs::classGPs(int me, const vector<double>& xyz, double weight, double Jacobian, int ndim) : 

    classSpatialPoint(xyz, ndim,me),  w(weight),J(abs(Jacobian)),
    _wJ(weight*abs(Jacobian)),_currentState(ndim),_previousState(ndim),_initialState(ndim),
    activate(1),flagFIni(0), _tmpState(NULL)
{

};

/*! This function allows to get the local Dofs at GP
  @param[in] numNodes Number of nodes in the domain
  @param[out] vdofs vector of Dofs 
*/
void classGPs::getLocalDofs(int numNodes, std::vector<int>& vdofs) const
{
    int numDofMechPerNode = _currentState.getNumMechanicalDofsPerNode();
    int numExtraDofPerNode = _currentState.getNumExtraDofsPerNode();
    //
    int numDofPerNode = _currentState.getTotalDofsPerNode();
    int numNeighbours = neighbours.size();
    vdofs.resize(numDofPerNode*numNeighbours,0);
    int nstoch = numDofMechPerNode/ndim;
    // mech dofs
    for (int i = 0; i < ndim; i++) 
    {
        for (int m = 0; m < nstoch; m++) 
        {
            for (int a=0; a < numNeighbours; a++)
            {
                int indi11 = a + numNeighbours * i + m * numNeighbours * ndim;
                vdofs[indi11] = neighbours[a] + numNodes * i + m * numNodes * ndim;
            }
        }
    }
    // extraDof
    for (int i=numDofMechPerNode; i < numDofPerNode; i++)
    {
        for (int a=0; a < numNeighbours; a++)
        {
            int indi11 = a + numNeighbours * i;
            vdofs[indi11] = neighbours[a] + numNodes * i;
        }
    }
};

void classGPs::setConsModel(vector<constitutiveModels *>& mlaws)
{
    if (mlaws.size() == 0)
    {
        //WARNING("no constitutive laws at this GP");
        return;
    }
    int numDofPerNodes = 0;
    for (int i=0; i< mlaws.size(); i++)
    {
        _allaws.addLaw(mlaws[i]); // law per field
        numDofPerNodes += mlaws[i]->getNbrDofConsMod(); // number dof per field
    };
    
    int numMech = mlaws[0]->getNbrDofConsMod();
    if (numMech > ndim)
    {
        _currentState.allocateMechanicalField(numMech);
    }
    int numExtra = numDofPerNodes - numMech;
    if (numExtra >0)
    {
        _currentState.allocateExtraField(numExtra);
    }
    _previousState = (_currentState);
    _initialState = (_currentState);
};

void classGPs::makeACheckPoint()
{
    if (_tmpState == NULL)
    {
        _tmpState = new classStates(_currentState);
    }
    (*_tmpState) = (_currentState);
}

void classGPs::restoreACheckPoint()
{
    if (_tmpState == NULL)
    {
        ERROR("no checkpoint exists");
    }
    _currentState = (*_tmpState);
}