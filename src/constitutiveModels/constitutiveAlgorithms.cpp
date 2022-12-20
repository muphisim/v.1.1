//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#include "constitutiveAlgorithms.h"
#include "classGPs.h"

/*! \brief This function to get the index of Dofs at the Gauss point in the local vector
  @param[in] GaussPoint Gauss Point
  @param[out] dofPosition Array of all indexes of the related dofs define by terms
*/
void MechanicalTerm::getDofsLocalPosition(const classGPs *GaussPoint, vector<int>& dofPosition) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    dofPosition.resize(numberMechDofPerNode*numNeighbours);
    for (int i=0; i< numberMechDofPerNode*numNeighbours; i++)
    {
        dofPosition[i] = i;
    }
}
/*! \brief This function to get the internal force
  @param[in] GaussPoint Gauss Point
  @param[out] vForce force vector
*/
void MechanicalTerm::getForceVector(const classGPs *GaussPoint, mVector& vForce) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    // stress
    const vector<double>& PK1 = currState.getFirstPiolaKirchhoffStress();
    // internal force
    vForce.resize(numNeighbours*numberMechDofPerNode,true);
    for (int a = 0; a < numNeighbours; a++)
    {
        for (int i = 0; i < ndim; i++) 
        {
            for (int J = 0; J < ndim; J++) 
            {
                vForce(a+i*numNeighbours) += wJ * PK1[i * ndim + J] * dphi[a * ndim + J];;
            }

        }
    };
}

/*! \brief This function to get the stiffness matrix
  @param[in] GaussPoint Gauss Point
  @param[out] stiff stiffness matrix
*/
void MechanicalTerm::getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    const classTensor4* tanModuli = currState.getTangent();
    int totalNumDofs = numberMechDofPerNode*numNeighbours;
    stiff.resize(totalNumDofs,totalNumDofs,true);
    // update
    for (int a = 0; a < numNeighbours; a++)
    {
        for (int b = 0; b < numNeighbours; b++) 
        {
            for (int ii = 0; ii < ndim; ii++) 
            {
                int indi1 = a + numNeighbours * ii;
                for (int kk = 0; kk < ndim; kk++) 
                {
                    int indi2 = b + numNeighbours * kk;
                    for (int JJ = 0; JJ < ndim; JJ++) 
                    {
                        for (int LL = 0; LL < ndim; LL++) 
                        {
                            stiff(indi1,indi2) += wJ *tanModuli->get(ii, JJ, kk, LL) * dphi[a * ndim + JJ] * dphi[b * ndim + LL];
                        }
                    }
                }
            }
        }
    }
};

/*! \brief This function is to get the mass matrix
  @param[in] GaussPoint Gauss Point
  @param[out] stiff stiffness matrix
*/
void MechanicalTerm::getMassMatrix(const classGPs *GaussPoint, mMatrix& mass) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    int totalNumDofs = numberMechDofPerNode*numNeighbours;
    mass.resize(totalNumDofs,totalNumDofs,true);
    double rho = GaussPoint->getConstitutiveManager().getConsModel()[0]->getRho();
    for (int a = 0; a < numNeighbours; a++) 
    {
        for (int b = 0; b < numNeighbours; b++) 
        {
            for (int ii = 0; ii < ndim; ii++) 
            {
                int indi1 = a + numNeighbours * ii;
                int indi2 = b + numNeighbours * ii;
                {
                    mass(indi1,indi2) += wJ * rho * phi[a] * phi[b];
                }
            }
        }
    }
};

/*! \brief This function to get the lumped mass vector
          @param[in] GaussPoint Gauss Point
          @param[out] mass mass vector
        */
void MechanicalTerm::getLumpedMassVector(const classGPs *GaussPoint, mVector& mass) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    string type = GaussPoint->getElementType();
    double rho = GaussPoint->getConstitutiveManager().getConsModel()[0]->getRho();
    const vector<int> neighbours = GaussPoint->getNeighbours(); // GaussPointoral list of nodes in which a given GPs lies
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    int totalNumDofs = numberMechDofPerNode*numNeighbours;
    mass.resize(totalNumDofs,true);    
    if (type == "CPS6") 
    {
        double Jacobian = GaussPoint->getJ();
        double nbGPperEl = 3;
        for (int a = 0; a < numNeighbours; a++) {
            for (int ii = 0; ii < ndim; ii++) {
                int indi1 = a + numNeighbours * ii;
                for (int kk = 0; kk < ndim; kk++) {
                    if (ii == kk) {
                        if (a < 3) {
                            mass(indi1) += 0.5 * Jacobian * rho * 1.0 /
                                            (12.0 * nbGPperEl); //area=isoparametric area * jacobian
                            
                        } else {
                            mass(indi1) += 0.5 * Jacobian * rho * 0.25 / nbGPperEl;
                        }

                    }
                }
            }
        }
    } 
    else if (type == "C3D10") 
    {
        double Jacobian = GaussPoint->getJ();
        double nbGPperEl = 4;
        for (int a = 0; a < numNeighbours; a++) {
            for (int ii = 0; ii < ndim; ii++) {
                int indi1 = neighbours[a] + numNeighbours * ii;
                for (int kk = 0; kk < ndim; kk++) {
                    if (ii == kk) {
                        if (a < 4) {
                            mass(indi1) += 1.0 / 6.0 * Jacobian * rho * 1.0 / (32.0 * nbGPperEl);
                        } else {
                            mass(indi1) += 1.0 / 6.0 * Jacobian * rho * 7.0 /
                                            (48.0 * nbGPperEl);//area = isoparametric area*jacobian determinant
                        }

                    }
                }
            }
        }
    } 
    else
    {
        for (int a = 0; a < numNeighbours; a++) {
            for (int b = 0; b < numNeighbours; b++) {
                for (int ii = 0; ii < ndim; ii++) {
                    int indi1 = a + numNeighbours * ii;
                    for (int kk = 0; kk < ndim; kk++) {
                        if (ii == kk) {
                            mass(indi1) += wJ * rho * phi[a] *  phi[b];
                        }
                    }
                }
            }
        }
    }
};

/*! \brief This function to get fNint in explicit scheme with extraDof
  @param[in] GaussPoint Gauss Point
  @param[out] vForce force vector
  @param[in] timeRun current time
  @param[in] dt time step
*/
void MechanicalTerm::getForceVectorExplicitScheme(const classGPs *GaussPoint, mVector& vForce, double timeRun, double dt) const
{
    // the same as forcevector
    getForceVector(GaussPoint, vForce);
};


/*! \brief This function to get the index of Dofs at the Gauss point in the local vector
  @param[in] GaussPoint Gauss Point
  @param[out] dofPosition Array of all indexes of the related dofs define by terms
*/
void MechanicalTermStochastic::getDofsLocalPosition(const classGPs *GaussPoint, vector<int>& dofPosition) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    dofPosition.resize(numberMechDofPerNode*numNeighbours);
    for (int i=0; i< numberMechDofPerNode*numNeighbours; i++)
    {
        dofPosition[i] = i;
    }
};

/*! \brief This function to get the internal force
  @param[in] GaussPoint Gauss Point
  @param[out] vForce force vector
*/
void MechanicalTermStochastic::getForceVector(const classGPs *GaussPoint, mVector& vForce) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    // stress
    const vector<double>& PK1 = currState.getFirstPiolaKirchhoffStress();
    
    int ndim = GaussPoint->getDimension();
    int nstoch = numberMechDofPerNode/ndim;
    
    vForce.resize(numNeighbours*numberMechDofPerNode,true);
    for (int a = 0; a < numNeighbours; a++) 
    {
        for (int i = 0; i < ndim; i++) 
        {
            for (int m = 0; m < nstoch; m++) 
            {
                int indi11 = a + numNeighbours * i + m * numNeighbours * ndim;
                for (int J = 0; J < ndim; J++) 
                {
                    vForce(indi11) += wJ * PK1[m * ndim * ndim + i * ndim + J] * dphi[a * ndim + J];
                }
            }
        }
    };
}

/*! \brief This function to get the stiffness matrix
  @param[in] GaussPoint Gauss Point
  @param[out] stiff stiffness matrix
*/
void MechanicalTermStochastic::getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    int totalNumDofs = numberMechDofPerNode*numNeighbours;
    const classTensor6* tanModuli = currState.getTangentStochastic();
    
    int ndim = GaussPoint->getDimension();
    int nstoch = numberMechDofPerNode/ndim;
    
    stiff.resize(totalNumDofs,totalNumDofs,true);
    for (int a = 0; a < numNeighbours; a++) {
        for (int b = 0; b < numNeighbours; b++) {
            for (int mm = 0; mm < nstoch; mm++) {
                for (int nn = 0; nn < nstoch; nn++) {
                    for (int ii = 0; ii < ndim; ii++) {
                        int indi1 = a + numNeighbours * ii + mm * numNeighbours * ndim;
                        for (int kk = 0; kk < ndim; kk++) {
                            int indi2 = b + numNeighbours * kk + nn * numNeighbours * ndim;
                            for (int JJ = 0; JJ < ndim; JJ++) {
                                for (int LL = 0; LL < ndim; LL++) {
                                    stiff(indi1,indi2) += wJ*(tanModuli->get(mm, nn, ii, JJ, kk, LL) *dphi[a*ndim + JJ] * dphi[b*ndim + LL]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/*! \brief This function is to get the mass matrix
  @param[in] GaussPoint Gauss Point
  @param[out] stiff stiffness matrix
*/
void MechanicalTermStochastic::getMassMatrix(const classGPs *GaussPoint, mMatrix& mass) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    int totalNumDofs = numberMechDofPerNode*numNeighbours;
    mass.resize(totalNumDofs,totalNumDofs,true);
    double rho = GaussPoint->getConstitutiveManager().getConsModel()[0]->getRho();
    for (int a = 0; a < numNeighbours; a++) 
    {
        for (int b = 0; b < numNeighbours; b++) 
        {
            for (int ii = 0; ii < ndim; ii++) 
            {
                int indi1 = a + numNeighbours * ii;
                int indi2 = b + numNeighbours * ii;
                mass(indi1,indi2) += wJ * rho * phi[a] * phi[b];
            }
        }
    }
};

/*! \brief This function to get the lumped mass vector
          @param[in] GaussPoint Gauss Point
          @param[out] mass mass vector
        */
void MechanicalTermStochastic::getLumpedMassVector(const classGPs *GaussPoint, mVector& mass) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    string type = GaussPoint->getElementType();
    double rho = GaussPoint->getConstitutiveManager().getConsModel()[0]->getRho();
    const vector<int> neighbours = GaussPoint->getNeighbours(); // GaussPointoral list of nodes in which a given GPs lies
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    int totalNumDofs = numberMechDofPerNode*numNeighbours;
    mass.resize(totalNumDofs,true);    
    
    int nstoch = numberMechDofPerNode/ndim;
    if (type == "CPS6")
    {
        double Jacobian = GaussPoint->getJ();
        double nbGPperEl = 3;
        for (int a = 0; a < numNeighbours; a++) {
            for (int ii = 0; ii < ndim * nstoch; ii++) {
                int indi1 = neighbours[a] + numNeighbours * ii;
                for (int kk = 0; kk < ndim * nstoch; kk++) {
                    if (ii == kk) {
                        if (a < 3) {
                            mass(indi1) += 0.5 * Jacobian * rho * 1.0 /
                                            (12.0 * nbGPperEl); //area=isoparametric area * jacobian
                        } else {
                            mass(indi1) += 0.5 * Jacobian * rho * 0.25 / nbGPperEl;
                        }

                    }
                }
            }
        }
    }
    else if (type == "C3D10") 
    {
        double Jacobian = GaussPoint->getJ();
        double nbGPperEl = 4;
        for (int a = 0; a < numNeighbours; a++) {
            for (int ii = 0; ii < ndim * nstoch; ii++) {
                int indi1 = neighbours[a] + numNeighbours * ii;
                for (int kk = 0; kk < ndim * nstoch; kk++) {
                    if (ii == kk) {
                        if (a < 4) {
                            mass(indi1) += 1.0 / 6.0 * Jacobian * rho * 1.0 / (32.0 * nbGPperEl);
                        } else {
                            mass(indi1) += 1.0 / 6.0 * Jacobian * rho * 7.0 /
                                            (48.0 * nbGPperEl);//area = isoparametric area*jacobian determinant
                        }

                    }
                }
            }
        }
    } 
    else 
    {
        for (int a = 0; a < numNeighbours; a++) {
            for (int b = 0; b < numNeighbours; b++) {
                double phia = phi[a];
                double phib = phi[b];
                for (int ii = 0; ii < ndim * nstoch; ii++) {
                    int indi1 = a + numNeighbours * ii;
                    for (int kk = 0; kk < ndim * nstoch; kk++) {
                        if (ii == kk) {
                            mass(indi1) += wJ * rho * phia * phib;
                        }
                    }
                }
            }
        }
    }
}

/*! \brief This function to get fNint in explicit scheme with extraDof
  @param[in] GaussPoint Gauss Point
  @param[out] vForce force vector
  @param[in] timeRun current time
  @param[in] dt time step
*/
void MechanicalTermStochastic::getForceVectorExplicitScheme(const classGPs *GaussPoint, mVector& vForce, double timeRun, double dt) const
{
    // the same as force vector
    getForceVector(GaussPoint, vForce);
};

/*! \brief constructor
  @param[in] fieldIndex index of the extraDof field
*/
ExtraDofTerm::ExtraDofTerm(int fieldIndex): _fieldIndex(fieldIndex){}

/*! \brief This function to get the index of Dofs at the Gauss point in the local vector
  @param[in] GaussPoint Gauss Point
  @param[out] dofPosition Array of all indexes of the related dofs define by terms
*/
void ExtraDofTerm::getDofsLocalPosition(const classGPs *GaussPoint, vector<int>& dofPosition) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    dofPosition.resize(numNeighbours);
    for (int i=0; i< numNeighbours; i++)
    {
        dofPosition[i] = (numberMechDofPerNode+_fieldIndex)*numNeighbours+i;
    };
    
};
/*! \brief This function to get the internal force
  @param[in] GaussPoint Gauss Point
  @param[out] vForce force vector
*/
void ExtraDofTerm::getForceVector(const classGPs *GaussPoint, mVector& vForce) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    
    double source = currState.getExtraFieldSources()[_fieldIndex];
    const classTensor1& flux = currState.getExtraFieldFluxes()[_fieldIndex];
    
    vForce.resize(numNeighbours,true);
    
    for (int a = 0; a < numNeighbours; a++)
    {
        vForce(a) += wJ*source*phi[a];
        for (int i=0; i< ndim; i++)
        {
            vForce(a) -= wJ * flux(i) * dphi[a * ndim + i];
        };
    };
}
/*! \brief This function to get the stiffness matrix
  @param[in] GaussPoint Gauss Point
  @param[out] stiff stiffness matrix
*/
void ExtraDofTerm::getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    
    double DsourceDExtraDof = (*currState.getDextraFieldSourcesDextraFields())[_fieldIndex][_fieldIndex];
    const classTensor1& DfluxDExtraDof = (*currState.getDextraFieldFluxesDextraFields())[_fieldIndex][_fieldIndex];
    const classTensor2& DfluxDgradExtraDof = (*currState.getDextraFieldFluxesDextraFieldGrads())[_fieldIndex][_fieldIndex];
    
    stiff.resize(numNeighbours,numNeighbours,true);
    // no coupling with other fields
    for (int a = 0; a < numNeighbours; a++)
    {
        for (int b =0; b < numNeighbours; b++)
        {
            stiff(a,b) += wJ*DsourceDExtraDof*phi[a]*phi[b];
            for (int i=0; i< ndim; i++)
            {
                stiff(a,b) -= wJ * DfluxDExtraDof(i) * dphi[a * ndim + i]*phi[b];
                for (int j=0; j< ndim; j++)
                {
                    stiff(a,b) -= wJ * DfluxDgradExtraDof(i,j) * dphi[a * ndim + i]*dphi[b * ndim + j];
                }
            };
            
        }
    };
};

/*! \brief This function to get the mass matrix
  @param[in] GaussPoint Gauss Point
  @param[out] mass mass matrix
*/
void ExtraDofTerm::getMassMatrix(const classGPs *GaussPoint, mMatrix& mass) const
{
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    mass.resize(numNeighbours,numNeighbours,true);
};

/*! \brief This function to get the lumped mass vector
          @param[in] GaussPoint Gauss Point
          @param[out] mass mass vector
        */
void ExtraDofTerm::getLumpedMassVector(const classGPs *GaussPoint, mVector& mass) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int> neighbours = GaussPoint->getNeighbours(); // GaussPointoral list of nodes in which a given GPs lies
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    mass.resize(numNeighbours,true);    
    
    double fieldDensity = currState.getExtraFieldDensities()[_fieldIndex];

    for (int a = 0; a < numNeighbours; a++)
    {
        for (int b = 0; b < numNeighbours; b++) 
        {
            mass(a) += wJ * fieldDensity * phi[a] * phi[b];
        }
    }
}

/*! \brief This function to get fNint in explicit scheme with extraDof
  @param[in] GaussPoint Gauss Point
  @param[out] vForce force vector
  @param[in] timeRun current time
  @param[in] dt time step
*/
void ExtraDofTerm::getForceVectorExplicitScheme(const classGPs *GaussPoint, mVector& vForce, double timeRun, double dt) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    
    double source = currState.getExtraFieldSources()[_fieldIndex];
    const classTensor1& flux = currState.getExtraFieldFluxes()[_fieldIndex];
    
    vForce.resize(numNeighbours,true);
    double T = currState.getExtraFields()[_fieldIndex];
    double Tprev = prevState.getExtraFields()[_fieldIndex];
    double fieldDensity = currState.getExtraFieldDensities()[_fieldIndex];
    
    double onlyFieldSource = (T-Tprev)*fieldDensity/dt;
    // source includes already the term onlyFieldSource
    for (int a = 0; a < numNeighbours; a++)
    {
        vForce(a) += wJ*(source-onlyFieldSource -Tprev*fieldDensity/dt)*phi[a];
        for (int i=0; i< ndim; i++)
        {
            vForce(a) -= wJ * flux(i) * dphi[a * ndim + i];
        };
    };
}


/*! \brief constructor
  @param[in] fieldIndexes indexes of the extraDof fields that the coupling occurs
*/
MechanicalExtraDofFullCouplingTerm::MechanicalExtraDofFullCouplingTerm(const vector<int>& fieldIndexes): 
            MechanicalTerm(), _fieldIndexes(fieldIndexes){}



/*! \brief This function to get the stiffness matrix
  @param[in] GaussPoint Gauss Point
  @param[out] stiff stiffness matrix
*/
void MechanicalExtraDofFullCouplingTerm::getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const
{
    
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    int numExtraDofs = currState.getNumExtraDofsPerNode();
    int totalDofs = (numberMechDofPerNode+numExtraDofs)*numNeighbours;
    stiff.resize(numberMechDofPerNode*numNeighbours,totalDofs,true);
    // add the mechanical part
    const classTensor4* tanModuli = currState.getTangent();
    int nbFields = _fieldIndexes.size();
    const vector<classTensor2>& DPK1DT = *(currState.getDfirstPiolaKirchhoffDextraFields());
    for (int a = 0; a < numNeighbours; a++)
    {
        for (int b = 0; b < numNeighbours; b++) 
        {
            for (int ii = 0; ii < ndim; ii++) 
            {
                int indi1 = a + numNeighbours * ii;
                // mecha-mecha
                for (int kk = 0; kk < ndim; kk++) 
                {
                    int indi2 = b + numNeighbours * kk;
                    for (int JJ = 0; JJ < ndim; JJ++) 
                    {
                        for (int LL = 0; LL < ndim; LL++) 
                        {
                            stiff(indi1,indi2) += wJ *tanModuli->get(ii, JJ, kk, LL) * dphi[a * ndim + JJ] * dphi[b * ndim + LL];
                        }
                    }
                    
                }
                // mecha-extraDof
                for (int kk =0; kk< nbFields; kk++)
                {
                    int field = _fieldIndexes[kk];
                    int indi2 = b + numNeighbours * (numberMechDofPerNode+field);
                    for (int JJ = 0; JJ < ndim; JJ++) 
                    {
                        stiff(indi1, indi2) += wJ * DPK1DT[field](ii,JJ) * dphi[a * ndim + JJ]*phi[b];
                    }
                }
                
            }
        }
    }
};

/*! \brief constructor
  @param[in] fieldIndex index of the extraDof field
*/
ExtraDofMechanicalFullCouplingTerm::ExtraDofMechanicalFullCouplingTerm(int fieldIndex): ExtraDofTerm(fieldIndex){}


/*! \brief This function to get the stiffness matrix
  @param[in] GaussPoint Gauss Point
  @param[out] stiff stiffness matrix
*/
void ExtraDofMechanicalFullCouplingTerm::getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<int>& neighbours = GaussPoint->getNeighbours();
    int numNeighbours = neighbours.size();
    const vector<double>& phi = GaussPoint->getPhi();
    const vector<double>& dphi = GaussPoint->getDphi();
    double wJ = GaussPoint->getWeightJ();
    int ndim = GaussPoint->getDimension();
    int numberMechDofPerNode = currState.getNumMechanicalDofsPerNode();
    int numExtraDofPerNode = currState.getNumExtraDofsPerNode();
    int totalDofs = (numberMechDofPerNode+numExtraDofPerNode)*numNeighbours;
    
    // couling with F and other fields
    const vector<double>& DsourceDExtraDof = (*currState.getDextraFieldSourcesDextraFields())[_fieldIndex];
    const classTensor2& DsourceDF = (*currState.getDextraFieldSourcesDdeformationGradient())[_fieldIndex];
    const vector<classTensor1>& DfluxDExtraDof = (*currState.getDextraFieldFluxesDextraFields())[_fieldIndex];
    const vector<classTensor2>& DfluxDgradExtraDof = (*currState.getDextraFieldFluxesDextraFieldGrads())[_fieldIndex];
    const classTensor3& DfluxDF = (*currState.getDextraFieldFluxesDdeformationGradient())[_fieldIndex];
    
    stiff.resize(numNeighbours,totalDofs,true);
    // no coupling with other fields
    for (int a = 0; a < numNeighbours; a++)
    {
        for (int b =0; b < numNeighbours; b++)
        {
            // extraDof- extraDof
            for (int kk =0; kk < numExtraDofPerNode; kk++)
            {
                int indi2 = b + (numberMechDofPerNode+kk)*numNeighbours;
            
                stiff(a,indi2) += wJ*DsourceDExtraDof[kk]*phi[a]*phi[b];
                for (int i=0; i< ndim; i++)
                {
                    stiff(a,indi2) -= wJ * DfluxDExtraDof[kk](i) * dphi[a * ndim + i]*phi[b];
                    for (int j=0; j< ndim; j++)
                    {
                        stiff(a,indi2) -= wJ * DfluxDgradExtraDof[kk](i,j) * dphi[a * ndim + i]*dphi[b * ndim + j];
                    }
                };
            }
            // extraDof meca
            for (int kk=0; kk< ndim; kk ++)
            {
                int indi2 = b + kk*numNeighbours;
                for (int j=0; j< ndim; j++)
                {
                    stiff(a,indi2) += wJ*DsourceDF(kk,j)*phi[a]*dphi[b * ndim + j];
                }
                for (int i=0; i< ndim; i++)
                {
                    for (int j=0; j< ndim; j++)
                    {
                        stiff(a,indi2) -= wJ * DfluxDF(i,kk,j) * dphi[a * ndim + i]*dphi[b * ndim + j];
                    }
                }
            }
            
        }
    };
};