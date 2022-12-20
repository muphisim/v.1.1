//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#include "constitutiveManager.h"
#include "classGPs.h"

constitutiveManager::constitutiveManager(){}
constitutiveManager::~constitutiveManager()
{
    
}

/*! \brief This function to clear all existing material laws.
*/
void constitutiveManager::clear()
{
    _allLaws.clear();
}
/*! \brief This function adds a material law to manager.
*/
void constitutiveManager::addLaw(constitutiveModels* law)
{
    // extraDof exists only from the second law
    if (_allLaws.size() > 0)
    {
        int dofMechPerNode = getNumMechDofPerNode();
        int nbDofPerNode = 0;
        for (int i=0; i< _allLaws.size(); i++)
        {
            nbDofPerNode += _allLaws[i]->getNbrDofConsMod();
        }
        law->clearExtraFieldIndexes();
        for (int i=nbDofPerNode; i < nbDofPerNode+law->getNbrDofConsMod(); i++)
        {
            law->addExtraField(i-dofMechPerNode);
        }
    }
    if (!law->isInitialised())
    {
        law->initLawFromOptions();
    }
    _allLaws.push_back(law);
}

/*! \brief This function return the total number of dofs of each node
 */
int constitutiveManager::getTotalNumDofPerNode() const
{
    int nbDofPerNode = 0;
    for (int i=0; i< _allLaws.size(); i++)
    {
        nbDofPerNode += _allLaws[i]->getNbrDofConsMod();
    }
    return nbDofPerNode;
}

/*! \brief This function return the total number of mechancal dofs of each node  
 */
int constitutiveManager::getNumMechDofPerNode() const
{
    return _allLaws[0]->getNbrDofConsMod();
}
/*! \brief This function return the total number of extra dofs of each node
 */
int constitutiveManager::getNumExtraDofPerNode() const
{
    int nbDofPerNode = 0;
    for (int i=1; i< _allLaws.size(); i++)
    {
        nbDofPerNode += _allLaws[i]->getNbrDofConsMod();
    }
    return nbDofPerNode;
}

/*! \brief This function initialises the internal state.
 */
void constitutiveManager::initIntVars(vector<double>& intVars) const
{
    intVars.clear();
    for (int i = 0; i < _allLaws.size(); i++) 
    {
        _allLaws[i]->initIntVars(intVars); //The loop takes the size of the internal variables vector of all the constitutive model used in the problem for each GP
    }
}

/*! \brief This function checks the activation at a GP.
  @param[inout] GaussPoint all data at Gauss point
  @param[in] timeRun Current time
*/
void constitutiveManager::checkActivation(classGPs *GaussPoint, double timeRun) const
{
    for (int i=0; i< _allLaws.size(); i++)
    {
        _allLaws[i]->checkActivation(GaussPoint,timeRun);
    }
}

/*! \brief This function updates the constitutive response.
  @param[inout] GaussPoint all data at Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void constitutiveManager::constitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) const
{
    // copy    
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    if (flagTanMod and (!currState.tangentAllocated()))
    {
        currState.allocateTangent();
    };
    currState.startFromOtherState(prevState);
    
    //add prediction first
    for (int i=0; i< _allLaws.size(); i++)
    {
        _allLaws[i]->predictIntVars(GaussPoint,dt,timeRun);
    }
    // constititive laws are estimated after checking all laws
    for (int i=0; i< _allLaws.size(); i++)
    {
        _allLaws[i]->updateConstitutive(GaussPoint,dt,timeRun,flagTanMod);
    }
};

/*! \brief This function computes kinematics avaiable.
  @param[inout] GaussPoint all data at Gauss point
  @param[in] uK unknown vector
  @param[in] numNodes number of nodes in the domain
*/
void constitutiveManager::computeKinematicVariables(classGPs *GaussPoint, const vector<double>& uK, int numNodes) const
{
    classStates& currState = GaussPoint->getCurrentState();
    // mechanics
    int ndim = currState.getDimension();
    int numMechDofs = currState.getNumMechanicalDofsPerNode();
    int numExtraDofs = currState.getNumExtraDofsPerNode();
    int numStochFields = numMechDofs/ndim -1;
    //
    const vector<double>& dphi = GaussPoint->getDphi();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    
    // mechanical
    vector<double>& F = currState.getDeformationGradient();
    setAll(F,0.);
    for (int i=0; i< ndim; i++)
    {
        F[i * ndim + i] += 1.;
        for (int a = 0; a < numNeighbours; a++) 
        {
            int uKloc = neighbours[a] + i * numNodes;
            for (int j = 0; j < ndim; j++) 
            {
                F[i * ndim + j] += dphi[a * ndim + j]* uK[uKloc];
            }
        }
    }
    
    if (numStochFields > 0)
    {
        for (int k=0; k < numStochFields; k++)
        {
            for (int a = 0; a < numNeighbours; a++) 
            {
                for (int i=0; i< ndim; i++)
                {
                    int uKloc =  neighbours[a] + numNodes * i + (k+1) * numNodes * ndim;
                    for (int j = 0; j < ndim; j++) 
                    {
                        F[i * ndim + j + (k+1)*ndim*ndim] += dphi[a * ndim + j]* uK[uKloc];
                    }
                }
            }
        }
    };
    // extraDof
    if (numExtraDofs > 0)
    {
        vector<double>& field = currState.getExtraFields();
        vector<classTensor1 >& gradField = currState.getExtraFieldGradients();
        const vector<double>& phi = GaussPoint->getPhi();
        for (int ifield = 0; ifield < numExtraDofs; ifield ++)
        {
            field[ifield] = 0;
            gradField[ifield].setAll(0);
            for (int a = 0; a < numNeighbours; a++) 
            {
                int uKloc = neighbours[a] + (numMechDofs+ ifield) * numNodes;
                field[ifield] += phi[a]*uK[uKloc];
                for (int j = 0; j < ndim; j++) 
                {
                    gradField[ifield](j) += dphi[a * ndim + j]* uK[uKloc];
                }
            }
        }
    };
};

/*! \brief This function computes the inertial force vector 
  @param[in] numNodes numbor of nodes in the domain
  @param[in] aK acceleration vector
  @param[in] GaussPoint Gauss Point
  @param[out] vDofs local dofs
  @param[out] vForce force vector
*/
void constitutiveManager::computeInertialForce(int numNodes, const vector<double>& aK, classGPs *GaussPoint, vector<int>& vDofs, mVector& vForce) const
{
    GaussPoint->getLocalDofs(numNodes,vDofs);
    int numDofs = vDofs.size();
    mVector acc(numDofs);
    for (int i=0; i< vDofs.size(); i++)
    {
        acc(i) = aK[vDofs[i]];
    }
    vForce.resize(numDofs,true);
    vector<int> vDofsPosition;
    mMatrix mass;
    mVector vForcePart, aPart;
    for (int i=0; i< _allLaws.size(); i++)
    {
        const Term* term = _allLaws[i]->getTerm();
        term->getDofsLocalPosition(GaussPoint, vDofsPosition);
        term->getMassMatrix(GaussPoint,mass);
        int partSize = vDofsPosition.size();
        aPart.resize(partSize,true);
        for (int j=0; j< partSize; j++)
        {
            aPart(j) =  acc(vDofsPosition[j]);
        }
        mass.mult(aPart,vForcePart);
        for (int j=0; j< partSize; j++)
        {
            vForce(vDofsPosition[j])+= vForcePart(j);
        }
    }
}

/*! \brief This function computes the internal force vector 
  @param[in] numNodes numbor of nodes in the domain
  @param[in] GaussPoint Gauss Point
  @param[out] vDofs local dofs
  @param[out] vForce force vector
*/
void constitutiveManager::computeInternalForce(int numNodes, classGPs *GaussPoint,vector<int>& vDofs, mVector& vForce) const
{
    GaussPoint->getLocalDofs(numNodes,vDofs);
    int numDofs = vDofs.size();
    vForce.resize(numDofs,true);
    vector<int> vDofsPosition;
    mVector vForcePart;
    for (int i=0; i< _allLaws.size(); i++)
    {
        const Term* term = _allLaws[i]->getTerm();
        term->getDofsLocalPosition(GaussPoint, vDofsPosition);
        term->getForceVector(GaussPoint,vForcePart);
        int partSize = vDofsPosition.size();
        for (int j=0; j< partSize; j++)
        {
            vForce(vDofsPosition[j])+= vForcePart(j);
        }
    }
};

/*! \brief This function computes the stiffness matrix
  @param[in] numNodes numbor of nodes in the domain
  @param[in] GaussPoint Gauss Point
  @param[out] vDofs local dofs
  @param[out] stiffness stiffness matrix
*/
void constitutiveManager::computeStiffnessMatrix(double timeRun, double dt, vector<double>& uK, int numNodes, 
                        classGPs *GaussPoint, vector<int>& vDofs, mMatrix& stiffness, bool byPerturbation, double tol) const
{
   
    
    if (byPerturbation)
    {
        mVector vForce;
        computeInternalForce(numNodes, GaussPoint, vDofs, vForce);
        int numDofs = vDofs.size();
        stiffness.resize(numDofs,numDofs,true);
        mVector vForcePlus;
        vector<int> vDofsTmp;
        // store current state
        classStates currState = GaussPoint->getCurrentState();
        for (int i=0; i< numDofs; i++)
        {
            double curVal = uK[vDofs[i]];
            uK[vDofs[i]] += tol;
            computeKinematicVariables(GaussPoint, uK, numNodes);
            constitutive(GaussPoint, dt, timeRun, false);
            computeInternalForce(numNodes, GaussPoint, vDofsTmp, vForcePlus);
            for (int j=0; j< numDofs; j++)
            {
                stiffness(j,i) = (vForcePlus(j)-vForce(j))/tol;
            }
            uK[vDofs[i]] = curVal;
        }
        // back to current satete
        GaussPoint->getCurrentState() = currState;
    }
    else
    {
        GaussPoint->getLocalDofs(numNodes,vDofs);
        int numDofs = vDofs.size();
        stiffness.resize(numDofs,numDofs,true);
        vector<int> vDofsPosition;
        mMatrix stiffPart;
        for (int i=0; i< _allLaws.size(); i++)
        {
            const Term* term = _allLaws[i]->getTerm();
            term->getDofsLocalPosition(GaussPoint,vDofsPosition);
            term->getStiffnessMatrix(GaussPoint,stiffPart);
            int partSize = vDofsPosition.size();
            if (term->withCoupling())
            {
                // 
                for (int j=0; j< partSize; j++)
                {
                    for (int k=0; k< numDofs; k++)
                    {
                        stiffness(vDofsPosition[j],k) +=  stiffPart(j,k);
                    }
                }
            }
            else
            {
                // square form
                for (int j=0; j< partSize; j++)
                {
                    for (int k=0; k< partSize; k++)
                    {
                        stiffness(vDofsPosition[j],vDofsPosition[k]) +=  stiffPart(j,k);
                    }
                }
                
            }
        }
    
    }
};

/*! \brief This function computes the mass matrix
  @param[in] numNodes numbor of nodes in the domain
  @param[in] GaussPoint Gauss Point
  @param[out] vDofs local dofs
  @param[out] mass mass matrix
*/
void constitutiveManager::computeMassMatrix(int numNodes, classGPs *GaussPoint, vector<int>& vDofs, mMatrix& mass) const
{
    GaussPoint->getLocalDofs(numNodes,vDofs);
    int numDofs = vDofs.size();
    mass.resize(numDofs,numDofs,true);
    vector<int> vDofsPosition;
    mMatrix stiffPart;
    for (int i=0; i< _allLaws.size(); i++)
    {
        const Term* term = _allLaws[i]->getTerm();
        term->getDofsLocalPosition(GaussPoint,vDofsPosition);
        term->getMassMatrix(GaussPoint,stiffPart);
        int partSize = vDofsPosition.size();
        for (int j=0; j< partSize; j++)
        {
            for (int k=0; k< partSize; k++)
            {
                mass(vDofsPosition[j],vDofsPosition[k]) +=  stiffPart(j,k);
            }
        }
    }
};

/*! \brief This function computes the lumped mass vector
  @param[in] numNodes numbor of nodes in the domain
  @param[in] GaussPoint Gauss Point
  @param[out] vDofs local dofs
  @param[out] mass mass matrix
*/
void constitutiveManager::computeLumpedMassVector(int numNodes, classGPs *GaussPoint, vector<int>& vDofs, mVector& mass) const
{
    GaussPoint->getLocalDofs(numNodes,vDofs);
    int numDofs = vDofs.size();
    mass.resize(numDofs,true);
    vector<int> vDofsPosition;
    mVector massPart;
    for (int i=0; i< _allLaws.size(); i++)
    {
        const Term* term = _allLaws[i]->getTerm();
        term->getDofsLocalPosition(GaussPoint,vDofsPosition);
        term->getLumpedMassVector(GaussPoint,massPart);
        int partSize = vDofsPosition.size();
        for (int j=0; j< partSize; j++)
        {
            mass(vDofsPosition[j]) +=  massPart(j);
        }
    };
};


/*! \brief This function computes the force vector for extraDof in explict 
  @param[in] numNodes numbor of nodes in the domain
  @param[in] GaussPoint Gauss Point
  @param[out] vDofs local dofs
  @param[out] vForce force vector
*/
void constitutiveManager::computeInternalForceExplicit(int numNodes, classGPs *GaussPoint,vector<int>& vDofs, mVector& vForce , double timeRun, double dt)  const
{
    GaussPoint->getLocalDofs(numNodes,vDofs);
    int numDofs = vDofs.size();
    vForce.resize(numDofs,true);
    vector<int> vDofsPosition;
    mVector vForcePart;
    // start from extraDof
    for (int i=0; i< _allLaws.size(); i++)
    {
        const Term* term = _allLaws[i]->getTerm();
        term->getDofsLocalPosition(GaussPoint, vDofsPosition);
        term->getForceVectorExplicitScheme(GaussPoint,vForcePart, timeRun,dt);
        int partSize = vDofsPosition.size();
        for (int j=0; j< partSize; j++)
        {
            vForce(vDofsPosition[j])+= vForcePart(j);
        }
    }
};
