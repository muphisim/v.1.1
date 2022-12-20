//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file isoMorphogeneticGrowth.cpp
  \brief This file contains all functions related to isotropic growth constitutive model. There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/
#include "isoMorphogeneticGrowth.h"
#include "maths.h"
#include "classGPs.h"
/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
  @param[in] Gc Growth multiplier
*/
classIsoMorphogeneticGrowth::classIsoMorphogeneticGrowth(int ndim, string name, string name_cons, double rho,
                                                         double Young, double nu, double Gc) : 
                                                         
            constitutiveModels(ndim, name,  name_cons, rho),
            _Young(Young),_nu(nu),_lambda((Young * nu) / ((1. + nu) * (1. - 2. * nu))),
            _mu(0.5 * Young / (1. + nu)),_Gc(Gc)
{

    
};

classIsoMorphogeneticGrowth::~classIsoMorphogeneticGrowth()
{
 
}

void classIsoMorphogeneticGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term !=NULL) delete _term;
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classIsoMorphogeneticGrowth::soundSpeed() const 
{
    double rho = this->getRho();
    double sound = sqrt(_Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classIsoMorphogeneticGrowth::initIntVars(vector<double> &intVars) {
    intVars.push_back(1.0);

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classIsoMorphogeneticGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();

    double theta = intVarsPrev[0];
    if (timeRun > 1E-16) {
        theta += _Gc * dt;
    } else {
        theta = intVarsPrev[0];
    }
    intVarsCurr[0] = theta;
    double Je = determinantTensor(FCurr, ndim) / theta;
    vector<double> Finv(ndim *ndim, 0);
    vector<double> FinvT(ndim *ndim, 0);
    inverse(FCurr, ndim, Finv);
    transpose(Finv, ndim, FinvT);
    double logJe = log(Je);
    double ndime = double(ndim + 1E-10);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            PK1[i * ndim + j] =
                    _mu * FCurr[i * ndim + j] / pow(theta, 2.0 / ndime) + (_lambda * logJe - _mu) * FinvT[i * ndim + j];
        }
    }
    computeCauchy(FCurr, ndim, PK1,Cauchy);
    computeVMS(Cauchy, ndim, VMS);
    if (flagTanMod) {
        vector<double> delta(ndim *ndim, 0);
        for (int i=0; i< ndim; i++)
        {
            delta[i * ndim + i]  = 1;
        }
        
        classTensor4 *tanModuli = currState.getTangent();
        double temp = 0;
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
                for (int k = 0; k < ndim; k++) {
                    for (int l = 0; l < ndim; l++) {
                        temp = 0;
                        temp = _lambda * Finv[j * ndim + i] * Finv[l * ndim + k] -
                               (_lambda * logJe - _mu) * Finv[l * ndim + i] * Finv[j * ndim + k] +
                               _mu * delta[i * ndim + k] * delta[j * ndim + l] / pow(theta, 2.0 / ndime);
                        tanModuli->setValues(i, j, k, l, temp);
                    }
                }
            }
        }
    } 
}


/*! \brief This function provides the growth tensor. It would used for external calls, as the growth part is unstressed in the deformation gradient decomposition.
  @param[in] FCurr Current deformation gradient tensor
  @param[in] FPrev Current deformation gradient tensor
  @param[in] intVarsPrev Array of the previous state of the internal variables
  @param[out] intVarsCurr Array of the current state of the internal variables
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
  @param[in] PK1 First Piola Kirchhoff stress tensor
  @param[in] tanModuli Tangent moduli
  @param[out] Fg Growth tensor
*/
void classIsoMorphogeneticGrowth::getFg(const vector<double> &FCurr, const vector<double> &FPrev, const vector<double> &intVarsPrev,
                                        vector<double> &intVarsCurr, double dt, double timeRun, bool flagTanMod,
                                        const vector<double> &PK1, classTensor4 *tanModuli, vector<double> &Fg) const {
    int ndim = getDim();
    double theta = intVarsPrev[0];
    if (timeRun > 1E-16) {
        theta += _Gc * dt;
    } else {
        theta = intVarsPrev[0];
    }
    intVarsCurr[0] = theta; // update internal var

    double ndime = double(ndim + 1E-10);
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fg[i * ndim + j] = pow(theta, 1.0 / ndime) * delta[i * ndim + j];
        }
    }
}

/*! \brief This function provides the elastic Green-Lagrangian strain tensor. Please not that this E is not the same than the one provided in the results of MuPhiSim (that one is calculated with the whole deformation gradient) 
  @param[in] ndim Dimension of the domain
  @param[out] E Green-Lagrangian strain tensor
  @param[in] FCurr Current deformation gradient tensor
  @param[in] FPrev Previous deformation gradient tensor
  @param[in] intVarsPrev Array of the internal variables
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] nStep Current number of time steps
*/
void classIsoMorphogeneticGrowth::getE(int ndim, vector<double> &E, const vector<double> &FCurr, const vector<double> &FPrev,
                                       const vector<double> &intVarsPrev, double dt, double timeRun, int nStep) const{
    double theta = intVarsPrev[0];
    vector<double> C(ndim * ndim, 0);
    vector<double> Cinv(ndim * ndim, 0);
    vector<double> sPiola(ndim * ndim, 0);
    vector<double> Esmall(ndim * ndim, 0);
    vector<double> Fe(ndim * ndim, 0);
    vector<double> Cauchy(ndim * ndim, 0);
    double ndime = double(ndim + 1E-10);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fe[i * ndim + j] = FCurr[i * ndim + j] / pow(theta, 1.0 / ndime);
        }
    }
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    multSTensor3FirstTranspose(Fe, ndim, ndim, Fe, ndim, C); // C= F^t F
    sumTensorTensor(C, delta, E, ndim, -1);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            E[i * ndim + j] = 0.5 * E[i * ndim + j];
        }
    }
};
