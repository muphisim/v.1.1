//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file volMechDrivenGrowth.cpp
  \brief This file contains all functions related to growth as a consequence of variations in the deformation gradient (mechanical growth). There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/
#include "volMechDrivenGrowth.h"
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
classVolMechDrivenGrowth::classVolMechDrivenGrowth(int ndim, string name, string name_cons, double rho, double Young,
                                                   double nu, double Gc) : 
                                                   constitutiveModels(ndim, name, name_cons, rho),
        _Young(Young),_nu(nu),
        _lambda((Young * nu) / ((1. + nu) * (1. - 2. * nu))), _mu(0.5 * Young / (1. + nu)),
        _Gc(Gc), _n0(ndim,0)

{
    if (ndim == 3) {
        _n0[0] = 0.0;
        _n0[1] = 0.0;
        _n0[2] = 1.0; //3D
    } else {
        _n0[0] = 1.0;
        _n0[1] = 0.0;
    }
}

classVolMechDrivenGrowth::~classVolMechDrivenGrowth()
{
}

void classVolMechDrivenGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classVolMechDrivenGrowth::soundSpeed() const
{
    double rho = this->getRho();
    double sound = sqrt(_Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classVolMechDrivenGrowth::initIntVars(vector<double> &intVars) {
    intVars.push_back(1.0);
}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classVolMechDrivenGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    // double theta=intVarsPrev[0];
    // if(timeRun>1E-16){
    //   theta+=Gc*dt;
    // }else{
    //   theta=intVarsPrev[0];
    // }
    // intVarsCurr[0] = theta;

    // double Je  = determinantTensor(FCurr,ndim)/theta;
    // vector<double> Finv(ndim*ndim,0);
    // vector<double> FinvT(ndim*ndim,0);
    // inverse(FCurr,ndim,Finv);
    // transpose(Finv,ndim,FinvT);
    // double logJe=log(Je);
    // double ndime=double(ndim+1E-10);

    // for(int i=0;i<ndim;i++){
    //   for(int j=0;j<ndim;j++){
    //     PK1[i*ndim+j]=mu * FCurr[i*ndim+j]/pow(theta,2.0/ndime) + (lambda * logJe - mu) * FinvT[i*ndim+j];
    //   }
    // }

    // if(flagTanMod){
    //   double temp=0;
    //   for(int i=0; i<ndim;i++){
    //     for(int j=0; j<ndim;j++){
    // 	for(int k=0; k<ndim;k++){
    // 	  for(int l=0; l<ndim;l++){
    // 	    temp=0;
    // 	    temp=lambda*Finv[j*ndim+i]*Finv[l*ndim+k]-(lambda*logJe-mu)*Finv[l*ndim+i]*Finv[j*ndim+k]+mu*delta[i*ndim+k]*delta[j*ndim+l]/pow(theta,2.0/ndime);
    // 	    tanModuli->setValues(i,j,k,l,temp);
    // 	  }
    // 	}
    //     }
    //   }
    // }
    double theta = intVarsPrev[0];
    if (timeRun > 1E-16) {
        theta += _Gc * dt;
    } else {
        theta = intVarsPrev[0];
    }
    intVarsCurr[0] = theta;
    vector<double> delta(ndim *ndim,0);
    for (int i=0; i< ndim; i++)
        delta[i*ndim+i] = 1;
    vector<double> Fe(ndim *ndim,0);
    vector<double> Fe2(ndim *ndim,0);
    vector<double> Fg(ndim *ndim,0);
    vector<double> Fginv(ndim *ndim,0);
    vector<double> Finv(ndim *ndim,0);
    vector<double> FinvT(ndim *ndim,0);
    inverse(FCurr, ndim, Finv);
    transpose(Finv, ndim, FinvT);
    vector<double> n(ndim, 0);
    vector<double> nXn0(ndim *ndim,0);
    vector<double> nXn(ndim *ndim,0);
    vector<double> n0Xn0(ndim *ndim,0);

    multTensorVector(FCurr, ndim, ndim, _n0, ndim, n);
    dyadicProduct(_n0, ndim, _n0, ndim, n0Xn0);
    dyadicProduct(n, ndim, _n0, ndim, nXn0);
    dyadicProduct(n, ndim, n, ndim, nXn);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fg[i * ndim + j] = delta[i * ndim + j] + (theta - 1) * n0Xn0[i * ndim + j];
        }
    }
    inverse(Fg, ndim, Fginv);
    multTensorTensor3(FCurr, ndim, ndim, Fginv, ndim, Fe2);
    vector<double> be(ndim *ndim, 0);
    vector<double> kirch(ndim *ndim,  0);
    multSTensor3SecondTranspose(Fe2, ndim, ndim, Fe2, ndim, be);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fe[i * ndim + j] = FCurr[i * ndim + j] + (1 / theta - 1) * nXn0[i * ndim + j];
        }
    }

    double logJe = log(determinantTensor(Fe, ndim));

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            PK1[i * ndim + j] = _mu * FCurr[i * ndim + j] + (_lambda * logJe - _mu) * FinvT[i * ndim + j] +
                                (1 / (theta * theta) - 1) * nXn[i * ndim + j] * FinvT[i * ndim + j];
        }
    }
    computeCauchy(FCurr, ndim, PK1,Cauchy);
    computeVMS(Cauchy, ndim, VMS);
    if (flagTanMod) {
        classTensor4 *tanModuli = currState.getTangent();
        double temp = 0;
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
                for (int k = 0; k < ndim; k++) {
                    for (int l = 0; l < ndim; l++) {
                        temp = 0;
                        temp = _lambda * Finv[j * ndim + i] * Finv[l * ndim + k] -
                               (_lambda * logJe - _mu) * Finv[l * ndim + i] * Finv[j * ndim + k] +
                               _mu * delta[i * ndim + k] * delta[j * ndim + l] -
                               _mu * (1 / (theta * theta) - 1) * Finv[l * ndim + i] * Finv[j * ndim + k] *
                               nXn[i * ndim + j];
                        tanModuli->setValues(i, j, k, l, temp);
                    }
                }
            }
        }
    } 
}
