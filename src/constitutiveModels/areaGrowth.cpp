//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file areaGrowth.cpp
  \brief This file contains all functions related to a area growth constitutive model. There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/
#include "areaGrowth.h"
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
classAreaGrowth::classAreaGrowth(int ndim, string name, string name_cons, double rho, double Young, double nu,
                                 double Gc) : constitutiveModels(ndim, name, name_cons, rho),
            _Young(Young),_nu(nu), _lambda((Young * nu) / ((1. + nu) * (1. - 2. * nu))),
            _mu(0.5 * Young / (1. + nu)),_Gc(Gc), _n0(ndim,0)
 {

    _n0[0] = 0;
    _n0[1] = 1;
    if (ndim == 3) {
        _n0[0] = 0;
        _n0[1] = 0;
        _n0[2] = 1;
    }
};

classAreaGrowth::~classAreaGrowth()
{

}

/*! \brief This function initialises this law from options
*/
void classAreaGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classAreaGrowth::soundSpeed() const 
{
    double rho = this->getRho();
    double sound = sqrt(_Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classAreaGrowth::initIntVars(vector<double> &intVars)
{
    intVars.push_back(1.0);
}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classAreaGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
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
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);

    double theta = intVarsPrev[0];
    if (timeRun > 1E-16) {
        theta += _Gc * dt;
    } else {
        theta = intVarsPrev[0];
    }
    intVarsCurr[0] = theta; // update internal var
    

    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    vector<double> PK1e(ndim *ndim, 0);
    vector<double> Fe(ndim *ndim, 0);
    vector<double> Feinv(ndim *ndim, 0);
    vector<double> FeinvT(ndim *ndim, 0);
    //double Je;
    vector<double> Fg(ndim *ndim, 0);
    vector<double> Fginv(ndim *ndim, 0);
    vector<double> Finv(ndim *ndim, 0);
    vector<double> FinvT(ndim *ndim, 0);
    inverse(FCurr, ndim, Finv);
    transpose(Finv, ndim, FinvT);
    vector<double> n(ndim, 0);
    vector<double> n0Xn0(ndim *ndim, 0); // Tensor produc
    vector<double> nXn0(ndim *ndim, 0); // Tensor produc
    vector<double> nXn(ndim *ndim, 0); // Tensor produc
    // Deformed surface normal#
    multTensorVector(FCurr, ndim, ndim, _n0, ndim, n);
    dyadicProduct(_n0, ndim, _n0, ndim, n0Xn0);
    dyadicProduct(n, ndim, _n0, ndim, nXn0);
    dyadicProduct(n, ndim, n, ndim, nXn);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fg[i * ndim + j] = sqrt(theta) * delta[i * ndim + j] + (1 - sqrt(theta)) * n0Xn0[i * ndim + j];
        }
    }
    inverse(Fg, ndim, Fginv);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fe[i * ndim + j] = 1 / sqrt(theta) * FCurr[i * ndim + j] + (1 - 1 / sqrt(theta)) * nXn0[i * ndim + j];
        }
    }
    inverse(Fe, ndim, Feinv);
    transpose(Feinv, ndim, FeinvT);

    double logJe = log(determinantTensor(Fe, ndim));
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            PK1e[i * ndim + j] = _mu * Fe[i * ndim + j] + (_lambda * logJe - _mu) * FeinvT[i * ndim + j];
        }
    }
    multSTensor3SecondTranspose(PK1e, ndim, ndim, Fginv, ndim, PK1);
    vector<double> FgFg(ndim *ndim, 0);
    multTensorTensor3(Fginv, ndim, ndim, Fginv, ndim, FgFg);
    computeCauchy(FCurr, ndim, PK1,Cauchy);
    computeVMS(Cauchy, ndim, VMS);
    // tangent moduli
    if (flagTanMod) 
    {
        classTensor4 *tanModuli = currState.getTangent();
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
                for (int k = 0; k < ndim; k++) {
                    for (int l = 0; l < ndim; l++) {
                        double temp = _mu * delta[i * ndim + k] * FgFg[j * ndim + l] +
                               _lambda * Finv[j * ndim + i] * Finv[l * ndim + k] -
                               (_lambda * logJe - _mu) * Finv[l * ndim + i] * Finv[j * ndim + k];
                        tanModuli->setValues(i, j, k, l, temp);
                    }
                }
            }
        }
    }
};

