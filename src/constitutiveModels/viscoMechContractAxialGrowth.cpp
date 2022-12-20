//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file viscoMechContractAxialGrowth.cpp
  \brief This file contains all functions related to viscoElastic growth constitutive model with contractility. It is a wrapper the viscoElastic model and the growth model. The growth model should be defined in this file, and there will not be any call to any growth constitutive model. The viscoElastic part is the original formulation proposed by Simo et al (nonlinear viscoelasticity, Simo J.C., Computational Inelasticity, Springer). There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/

#include "viscoMechContractAxialGrowth.h"
#include "maths.h"
#include "viscoElastic.h"
#include "classGPs.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] K Bulk modulus
  @param[in] N Number of internal variables. Number of branches in the generalised Maxwell model
  @param[in] mus Array of mus for N branches. Its size is N+1, because the initial position mus[0] is the for the branch with the spring alone (Einfi)
  @param[in] etas Array of etas for N branches. Its size is N 
  @param[in] Gc Growth multiplier
  @param[in] n0 Direction of growth in the reference configuration
  @param[in] Kp Polymarisation rate of the F-actin
  @param[in] Kd Depolymerisation rate of the F-actin
  @param[in] rhoMC Density of F-actin
  @param[in] c0 Contractility rate
  @param[in] S0 Stress threshold for contractility
*/
classViscoMechContractAxialGrowth::classViscoMechContractAxialGrowth(int ndim, string name, string name_cons,
                                                                     double rho, double K, int N, const vector<double>& mus,
                                                                     const vector<double>& etas, double Gc, const vector<double>& n0,
                                                                     double Kp, double Kd, double rhoMC, double c0,
                                                                     double S0) : 
                                        constitutiveModels(ndim, name, name_cons, rho) {

    this->viscoElastic = new classViscoElastic(ndim, name, name_cons, rho, K, N, mus, etas);
    vector<double> _mu = mus;
    vector<double> _eta = etas;
    this->_sum_mu = 0;
    for (int i = 0; i < _mu.size(); ++i) {
        _sum_mu += _mu[i];
    }


    this->Young = 0;
    this->nu = 0;
    this->lambda = 0;
    this->mu = 0;
    this->n0 = n0;
    // Initialization
    this->lambda = (Young * nu) / ((1. + nu) * (1. - 2. * nu));
    this->mu = 0.5 * Young / (1. + nu);
    this->Gc = Gc; // in this constitutive model this is the initial value of the growth multiplier
    this->_N = N;
    this->Kp = Kp;
    this->Kd = Kd;
    this->rhoMC = rhoMC;
    this->c0 = c0;
    this->S0 = S0;

}

classViscoMechContractAxialGrowth::~classViscoMechContractAxialGrowth()
{
    delete viscoElastic;
}

void classViscoMechContractAxialGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classViscoMechContractAxialGrowth::soundSpeed() const {
    // double rho=this->rho;
    // double Young= this->Young;
    // double sound=sqrt(Young/rho);

    double sound = sqrt((_K + 4 / 3 * _sum_mu) / _rho);

    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classViscoMechContractAxialGrowth::initIntVars(vector<double> &intVars) {
    /* This finite viscous constitutive model has as many internal variables as maxwell branches. Then, for each branch there are two tensors Hn and Sn. Then the internal variables array will be built as follows:
       ndim=3;
       TotalSize=N*2*ndim*ndim; That reads: N internal variables with two ndim order tensors
    */
    int ndim = getDim();
    vector<double> Hn(ndim *ndim, 0.);
    vector<double> Sn(ndim *ndim, 0.);
    vector<double> Fcn(ndim *ndim, 0.);// Previous contractile tensor
    for (int i = 0; i < ndim; i++) {
        Fcn[i * ndim + i] = 1.0;
    }
    int finalSize = _N * (Hn.size() + Sn.size());

    this->_finalSize = finalSize;
// It should be initialised to zero, these loops are just to illustrate how to access to the internal variables;
    for (int i = 0; i < _finalSize; i++) {
        intVars.push_back(0.);
    }

    _finalSize = _finalSize + 1;
    intVars.push_back(0); // Growth rate initialised to zero
    // and just at the end, I add the Gc int var of the growth part

    _finalSize = _finalSize + 1;
    intVars.push_back(this->Gc); // This internal variable is initialised to 1


    for (int i = 0; i < ndim * ndim; i++) {
        intVars.push_back(Fcn[i]);
        _finalSize++;
    }


}

                                                
 /*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classViscoMechContractAxialGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    const vector<double> &Fn = currState.getDeformationGradient();
    const vector<double> &F0 = prevState.getDeformationGradient();
    vector<double> &P = currState.getFirstPiolaKirchhoffStress();
    const vector<double> &PK1Prev = prevState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    computeGL(Fn, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    classTensor4 *tanModuli = currState.getTangent();
    classTensor4 *viscousTanModuli;
    viscousTanModuli = new classTensor4(ndim);
    classTensor4 *tanPart12;
    vector<double> I(ndim *ndim, 0.);
    tanPart12 = new classTensor4(ndim);
    classTensor4 *tanPart123;
    tanPart123 = new classTensor4(ndim);
    classTensor4 *tanPart122;
    tanPart122 = new classTensor4(ndim);
    for (int i = 0; i < ndim; i++) {
        I[i * ndim + i] = 1.; // defining the kronecker delta
    }

    //  double rhoMCReal=rhoMC*rho;
    double rhoMCReal = rhoMC;
    vector<double> n0Xn0(ndim *ndim,
    0);
    vector<double> Fa(ndim *ndim,
    0.);
    vector<double> Fainv(ndim *ndim,
    0.);
    vector<double> Fg(ndim *ndim,
    0.);
    vector<double> F0inv(ndim *ndim,
    0.);
    vector<double> Fv(ndim *ndim,
    0.);
    vector<double> Fginv(ndim *ndim,
    0.);
    vector<double> intVarsPrevGrowth(1, 0.), intVarsCurrGrowth(1, 0.);
    int numIntVarContract = ndim * ndim;
    vector<double> Fc(ndim *ndim,
    0.);
    vector<double> Fcn(ndim *ndim,
    0.);
    vector<double> intVarsPrevVisco(_finalSize - numIntVarContract - 2, 0.), intVarsCurrVisco(
            _finalSize - numIntVarContract - 2, 0.);

    // Organisation of the internal variables
    // The second last position of the whole array finalSize
    intVarsPrevGrowth[0] = intVarsPrev[_finalSize - 1 - numIntVarContract];
    for (int i = 0; i < _finalSize - numIntVarContract - 2; i++) {
        intVarsPrevVisco[i] = intVarsPrev[i];
        intVarsCurrVisco[i] = intVarsCurr[i];
    }
    Fcn.assign(&intVarsPrev[_finalSize - numIntVarContract],
               &intVarsPrev[_finalSize - numIntVarContract] + numIntVarContract);

    /////////////////////////////////////////////////////////////////////////////////
    // Calculation of Fg due to the axial growth
    double theta = intVarsPrevGrowth[0];
    dyadicProduct(n0, ndim, n0, ndim, n0Xn0);
    double Jn1 = determinantTensor(Fn, ndim);
    // Jn1=1;
    double prevTheta = 0;
    if (timeRun > 1E-16) {
        prevTheta = theta;
        //theta=(1-Kd*dt)*prevTheta+Jn1/rhoMCReal*Kp*dt;
        theta = (prevTheta + Jn1 * Kp * dt / rhoMCReal) / (1 + Kd * dt);

    } else {
        theta = intVarsPrevGrowth[0];
    }
    intVarsCurrGrowth[0] = theta; // update internal var
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fg[i * ndim + j] = I[i * ndim + j] + (theta - 1) * n0Xn0[i * ndim + j];
            //cout<<Fg[i*ndim+j]<<" ";
        }
        //    cout<<endl;
    }
    //////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // Contractility
    //  inverse(F0,ndim,F0inv);
    //multTensorTensor3(F0inv,ndim,ndim,PK1Prev,ndim, Sn);
    multSTensor3SecondTranspose(PK1Prev, ndim, ndim, F0, ndim, Cauchy);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Cauchy[i * ndim + j] = Cauchy[i * ndim + j] / (determinantTensor(F0, ndim));
        }
    }
    double Stmp = 0;
    computeVMS(Cauchy, ndim, VMS);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Stmp += Cauchy[i * ndim + j] * n0Xn0[i * ndim + j];
            // cout<<Fcn[i*ndim+j]<<" ";
        }
        //   //  cout<<endl;
    }
    // cout<< "Stress sc "<< Stmp<<endl;
    double c = 0;
    if (Stmp > S0 || Stmp <= 0) {
        c = 0;
        // cout<< "Stress sc "<< Stmp<< " S0 "<<S0<<endl;
    } else {
        c = c0 * (1 - Stmp / S0);
    }
    if (c < 0) {
        c = 0;
    }
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fc[i * ndim + j] = Fcn[i * ndim + j] - c * dt * n0Xn0[i * ndim + j];
            //  cout<<Fc[i*ndim+j]<<" ";
        }
        // cout<<endl;
    }
    if (ndim == 2) {
        Fc[ndim * ndim - 1] =
                1 / Fc[0];// I have to write a function to make it incompressible. This is a temporal solution
    } else {
        Fc[8] = sqrt(1 / Fc[0]);// I have to write a function to make it incompressible. This is a temporal solution
        Fc[4] = sqrt(1 / Fc[0]);// I have to write a function to make it incompressible. This is a temporal solution
    }
    ////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////
    // ViscoElasticity
    multTensorTensor3(Fc, ndim, ndim, Fg, ndim, Fa);
    inverse(Fa, ndim, Fainv);
    multTensorTensor3(Fn, ndim, ndim, Fainv, ndim, Fv);

    // This one will me the viscous tangent moduli
    static_cast<classViscoElastic *>(this->viscoElastic)->updateConstitutive(Fv, F0, intVarsPrevVisco, intVarsCurrVisco,
                                                                             dt, timeRun, flagTanMod, P, P,
                                                                             viscousTanModuli);

    // // // If Fg(F)
    tensorProductSecondTensorSecondTensor(Fainv, ndim, ndim, I, ndim, ndim, tanPart122);
    vector<double> tempFA(ndim *ndim,
    0.);
    vector<double> Fninv(ndim *ndim,
    0.);
    vector<double> FninvT(ndim *ndim,
    0.);
    double a = Jn1 * Kp * dt / (rhoMCReal * theta * (prevTheta + Jn1 * Kp * dt / rhoMCReal));
    multTensorTensor3(n0Xn0, ndim, ndim, Fn, ndim, tempFA);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            tempFA[i * ndim + j] = tempFA[i * ndim + j] * a;
        }
    }
    inverse(Fn, ndim, Fninv);
    transpose(Fninv, ndim, FninvT);
    tensorProductSecondTensorSecondTensor(tempFA, ndim, ndim, FninvT, ndim, ndim, tanPart123);
    sumFourthTensorFourthTensor(ndim, tanPart122, tanPart123, tanPart12);

    // if not Fg(F)
    //tensorProductSecondTensorSecondTensor(Fainv, ndim, ndim, I, ndim, ndim, tanPart12);

    multFourthTensorFourthTensor(ndim, viscousTanModuli, tanPart12, tanModuli);

    delete viscousTanModuli;
    delete tanPart12;
    delete tanPart122;
    delete tanPart123;

    ///////////////////////////////////////////////////////////////////////////////
    // Updating internal variables
    for (int i = 0; i < _finalSize - numIntVarContract - 2; i++) {
        intVarsCurr[i] = intVarsCurrVisco[i];
    }
    for (int i = 0; i < ndim * ndim; i++) {
        intVarsCurr[_finalSize - numIntVarContract + i] = Fc[i];
    }
    //  intVarsCurr[_finalSize-3]=intVarsPrev[_finalSize-3]+(theta-prevTheta)/dt;
    intVarsCurr[_finalSize - 2 - numIntVarContract] = (theta - prevTheta) / dt;
    intVarsCurr[_finalSize - 1 - numIntVarContract] = intVarsCurrGrowth[0];
}
