//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file hyperElasticStVenantKirchhoffThermoMechanics.cpp
  \brief This constitutive model is thought just for Thermo-Mechanical analysis in small deformations. It is a wrapper to couple Hyper Elastic St-Venant-Kirchhoff model
  with a thermomechanical couling term. In that scenario, all stresses are equivalent and thus in PK1, the Cauchy stress will be returned
  Mechanical part of Fully-Coupled Thermo-Mechanical Analysis
*/
#include "hyperElasticStVenantKirchhoffThermoMechanics.h"
#include "maths.h"
#include "classGPs.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young's modulus
  @param[in] nu Poisson ratio
  @param[in] stressState 
  @param[in] alpha Thermal expansion coefficients
  @param[in] theta_ref Reference temperature at which the materialis defined as undeformed in the absence of mechanical loads
  @param[in] fieldIndex Index of the temperature field in the list of extra dofs (0 by default)
  
*/
classHyperElasticStVenantKirchhoffThermoMechanics::classHyperElasticStVenantKirchhoffThermoMechanics(int ndim,
                                                                                                     string name,
                                                                                                     string name_cons,
                                                                                                     double rho,
                                                                                                     double Young,
                                                                                                     double nu,
                                                                                                     int stressState,
                                                                                                     double alpha, double theta_ref,
                                                                                                     int fieldIndex) : 
        classHyperElasticStVenantKirchhoff(ndim, name, name_cons, rho, Young, nu, stressState),
        _alphatensor(ndim), _theta_ref(theta_ref)
{
     _alphatensor.addDiagonal(alpha);
    _tempratureFieldIndex = fieldIndex;
    
}

classHyperElasticStVenantKirchhoffThermoMechanics::~classHyperElasticStVenantKirchhoffThermoMechanics() 
{
};

void classHyperElasticStVenantKirchhoffThermoMechanics::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalExtraDofFullCouplingTerm(std::vector<int>(1,_tempratureFieldIndex));
};


/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classHyperElasticStVenantKirchhoffThermoMechanics::initIntVars(vector<double> &intVars) {
    _intVarsSize = intVars.size();
    intVars.push_back(0.0);  //activation time
    intVars.push_back(0.0);  //activation state
};


/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classHyperElasticStVenantKirchhoffThermoMechanics::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    double CurrTemperature = currState.getExtraFields()[_tempratureFieldIndex];
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    
    intVarsCurr[_intVarsSize] = intVarsPrev[_intVarsSize]; //If it's inside this function, that means the element is active
    intVarsCurr[_intVarsSize + 1] = 1; //If it's inside this function, that means the element is active
    

    //we need to get the elastic stiffness tensor to calculate the thermomechanical coupling tensor m
    vector < double > mtilde(ndim * ndim, 0);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                for (int l = 0; l < ndim; l++) {
                    mtilde[i * ndim + j] += _elasticStiffTensor(i, j, k, l) * _alphatensor(k,l);
                }
            }
        }
    }
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    
    vector < double >  m = mtilde; //Thermomechanical coupling tensor = elasticStiffTensor*alphatensor.
    multiSTensorScalar(mtilde, ndim, ndim, _theta_ref - CurrTemperature); //The stress is divided in a pure mechanical part and a coupling part (mtilde)

    // mechanical part without temperature effect
    vector<double> PK2(ndim *ndim, 0);
    vector<double> rightCauchy(ndim *ndim, 0);
    vector<double> FCurrinv(ndim *ndim, 0);


    classHyperElasticStVenantKirchhoff::stress(PK2, E);
    multTensorTensor3(FCurr, ndim, ndim, PK2, ndim, PK1);
    

    inverse(FCurr, ndim, FCurrinv);
    setAll(PK2,0.);
    multTensorTensor3(FCurrinv, ndim, ndim, PK1, ndim, PK2); //we have PK1 and we need PK2
    sumTensorTensor(PK2, mtilde, PK2, ndim, 1); //Sum of pure mechanical and coupling parts
    setAll(PK1,0.);
    multTensorTensor3(FCurr, ndim, ndim, PK2, ndim, PK1); //We get the new PK1, included the coupling term
    computeCauchy(FCurr, ndim, PK1,Cauchy);
    computeVMS(Cauchy, ndim, VMS);

    //Below, we get some terms of the stiffness matrices for this thermomechanical constitutive model
    // COUPLING TENSOR aiJ
    vector < double > atensor(ndim * ndim, 0);
    for (int i = 0; i < ndim; i++) {
        for (int J = 0; J < ndim; J++) {
            for (int P = 0; P < ndim; P++) {
                for (int Q = 0; Q < ndim; Q++) {
                    atensor[i * ndim + J] += 0.5 * m[P * ndim + Q] * (delta[P * ndim + J] * FCurr[i * ndim + Q] +

                                                                      FCurr[i * ndim + P] * delta[Q * ndim + J]);
                }
            }
        }
    }
    
    // compute source
    vector<double> dF(ndim * ndim, 0);
    for (int i = 0; i < ndim; i++) 
    {
        for (int J = 0; J < ndim; J++) 
        {
            dF[i * ndim + J] += (FCurr[i * ndim + J]-FPrev[i * ndim + J])/dt;
        }
    }
    double& mecaSource = currState.getExtraFieldSources()[_tempratureFieldIndex];
    mecaSource += CurrTemperature*dotProduct(atensor,dF,ndim * ndim);
    
    
    if (flagTanMod)
    {
        classTensor4 *tanModuli = currState.getTangent();
        for (int i = 0; i < ndim; i++) {
            for (int J = 0; J < ndim; J++) {
                for (int k = 0; k < ndim; k++) {
                    for (int L = 0; L < ndim; L++) {
                        double val = 0;
                        for (int K = 0; K < ndim; K++) {
                            for (int I = 0; I < ndim; I++) {
                                val += _elasticStiffTensor(I, J, K, L) * FCurr[i * ndim + I] * FCurr[k * ndim + K];
                            }
                        }
                        if (i == k) {
                            val += PK2[J * ndim + L]+mtilde[J * ndim + L];
                        }
                        tanModuli->setValues(i, J, k, L, val);
                    }
                }
            }
        }
        
        classTensor2& dPdtheta =  (*currState.getDfirstPiolaKirchhoffDextraFields())[_tempratureFieldIndex];
        //dPiJ/dtheta
        for (int i = 0; i < ndim; i++) {
            for (int J = 0; J < ndim; J++) {
                for (int K = 0; K < ndim; K++) 
                {
                    dPdtheta(i,J) += -1 * FCurr[i * ndim + K] * m[K * ndim + J];
                }
            }
        }
        
        //daiJ/dFsT
        classTensor4 dadF(ndim);
        for (int i = 0; i < ndim; i++) {
            for (int J = 0; J < ndim; J++) {
                for (int P = 0; P < ndim; P++) {
                    for (int Q = 0; Q < ndim; Q++) {
                        for (int s = 0; s < ndim; s++) {
                            for (int T = 0; T < ndim; T++) {
                                dadF(i,J,s,T) += 0.5 * m[P * ndim + Q] *
                                                    (delta[P * ndim + J] *
                                                     delta[i * ndim + s] *
                                                     delta[Q * ndim + T]
                                                     + delta[i * ndim + s] *
                                                       delta[P * ndim + T] *
                                                       delta[Q * ndim + J]);
                            }
                        }
                    }
                }
            }
        }
        
        classTensor2& DsourceDF = (*currState.getDextraFieldSourcesDdeformationGradient())[_tempratureFieldIndex];
        double& DsourceDT = (*currState.getDextraFieldSourcesDextraFields())[_tempratureFieldIndex][_tempratureFieldIndex];
        
        DsourceDT += dotProduct(atensor,dF,ndim * ndim);
        for (int i = 0; i < ndim; i++) 
        {
            for (int J = 0; J < ndim; J++) 
            {
                DsourceDF(i,J) += CurrTemperature*atensor[i * ndim + J]/dt;
                for (int s = 0; s < ndim; s++) 
                {
                    for (int T = 0; T < ndim; T++) 
                    {
                        DsourceDF(i,J) += CurrTemperature*dadF(s,T,i,J)*dF[s* ndim + T];
                    }
                }
            }
        }
    }
};
