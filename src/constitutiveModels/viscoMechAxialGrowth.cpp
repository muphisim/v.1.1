//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file viscoMechAxialGrowth.cpp
  \brief This file contains all functions related to viscoElastic growth constitutive model. It is a wrapper the viscoElastic model and the growth model. The growth model should be defined in this file, and there will not be any call to any growth constitutive model. The viscoElastic part is the original formulation proposed by Simo et al (nonlinear viscoelasticity, Simo J.C., Computational Inelasticity, Springer). There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/

#include "viscoMechAxialGrowth.h"
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
  @param[in] Gc Initial growth multiplier
  @param[in] n0 Direction of growth in the reference configuration
  @param[in] Kp Polymarisation rate of the microtubules
  @param[in] Kd Depolymerisation rate of the microtubules
  @param[in] rhoMicro Density of microtubules 
*/
classViscoMechAxialGrowth::classViscoMechAxialGrowth(int ndim, string name, string name_cons, double rho, double K,
                                                     int N, vector<double> mus, vector<double> etas, double Gc,
                                                     vector<double> n0, double Kp, double Kd, double rhoMicro)
        : constitutiveModels(ndim, name, name_cons, rho) {

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
    this->Young = Young;
    this->nu = nu;
    this->lambda = (Young * nu) / ((1. + nu) * (1. - 2. * nu));
    this->mu = 0.5 * Young / (1. + nu);
    this->Gc0 = Gc; // in this constitutive model this is the initial value of the growth multiplier
    this->_N = N;
    this->Kp = Kp;
    this->Kd = Kd;
    this->_K = K;
    this->rhoMicro = rhoMicro;
}

classViscoMechAxialGrowth::~classViscoMechAxialGrowth()
{
    delete viscoElastic;
};

void classViscoMechAxialGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classViscoMechAxialGrowth::soundSpeed() const{
    //double rho=this->rho;
    // double Young= this->Young;
    //  double sound=sqrt(Young/rho);
    //double factornu = (1.-nu)/((1.+nu)*(1.-2.*nu));
    //return sqrt(Young*factornu/rho);

    /* Lecture Notes in Applied Mechanics
       Volume 2
       Series Editor
       Prof. Dr.-Ing. Friedrich Pfeiffer*/
    // Fast shear wave:
    double sound = sqrt((_K + 4 / 3 * _sum_mu) / _rho);
    // Slow shear wave:
    //double sound=sqrt(_sum_mu/rho);
    //cout<<"hello "<< sound<<endl;
    // The fastest one will be always more critical
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classViscoMechAxialGrowth::initIntVars(vector<double> &intVars) {
    /* This finite viscous constitutive model has as many internal variables as maxwell branches. Then, for each branch there are two tensors Hn and Sn. Then the internal variables array will be built as follows:
       ndim=3;
       TotalSize=N*2*ndim*ndim; That reads: N internal variables with two ndim order tensors
    */
    int ndim = getDim();
    vector<double> Hn(ndim *ndim, 0.);
    vector<double> Sn(ndim *ndim, 0.);
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
    intVars.push_back(this->Gc0); // This internal variable is initialised to 1
}


/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classViscoMechAxialGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &Fn = currState.getDeformationGradient();
    const vector<double> &F0 = prevState.getDeformationGradient();
    vector<double> &P = currState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    computeGL(Fn, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    classTensor4 *tanModuli = currState.getTangent();
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();

    classTensor4 *viscousTanModuli;
    viscousTanModuli = new classTensor4(ndim);
    classTensor4 *tanPart12;
    tanPart12 = new classTensor4(ndim);
    classTensor4 *tanPart123;
    tanPart123 = new classTensor4(ndim);
    classTensor4 *tanPart122;
    tanPart122 = new classTensor4(ndim);
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }

    //double rhoMicroReal=rhoMicro*rho;
    double rhoMicroReal = rhoMicro;
    vector<double> n0Xn0(ndim *ndim, 0);
    vector<double> Fg(ndim *ndim, 0.);
    vector<double> Fv(ndim *ndim, 0.);
    vector<double> Fginv(ndim *ndim, 0.);
    vector<double> intVarsPrevGrowth(1, 0.), intVarsCurrGrowth(1, 0.);
    vector<double> intVarsPrevVisco(_finalSize - 2, 0.), intVarsCurrVisco(_finalSize - 2, 0.);
    // The last position of the whole array finalSize
    intVarsPrevGrowth[0] = intVarsPrev[_finalSize - 1];


    for (int i = 0; i < _finalSize - 2; i++) {
        intVarsPrevVisco[i] = intVarsPrev[i];
        intVarsCurrVisco[i] = intVarsCurr[i];
    }

    /////////////////////////////////////////////////////////////////////////////////
    // Calculation of Fg due to the axial growth
    double theta = intVarsPrevGrowth[0];
    dyadicProduct(n0, ndim, n0, ndim, n0Xn0);
    //double Gc=1;
    //Gc=0;
    double Jn1 = determinantTensor(Fn, ndim);
    //cout<< "determinant Jn"<< Jn1<<endl;
    //Jn1=1;
    double prevTheta = 0;
    if (timeRun > 1E-16) {
        prevTheta = theta;
        //theta=(1-Kd*dt)*prevTheta+Jn1*Kp*dt/rhoMicroReal;
        theta = (prevTheta + Jn1 * Kp * dt / rhoMicroReal) / (1 + Kd * dt);
        //theta+=Gc*dt;
    } else {
        theta = intVarsPrevGrowth[0];
    }

    intVarsCurrGrowth[0] = theta; // update internal var
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fg[i * ndim + j] = delta[i * ndim + j] + (theta - 1) * n0Xn0[i * ndim + j];
            //  cout<<Fg[i*ndim+j]<<" ";
        }
        //  cout<<endl;
    }
    //////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////
    // ViscoElasticity
    inverse(Fg, ndim, Fginv);
    multTensorTensor3(Fn, ndim, ndim, Fginv, ndim, Fv);

    // This one will me the viscous tangent moduli
    static_cast<classViscoElastic *>(this->viscoElastic)->updateConstitutive(Fv, F0, intVarsPrevVisco, intVarsCurrVisco,
                                                                             dt, timeRun, flagTanMod, P, P,
                                                                             viscousTanModuli);
    Cauchy = P;
    computeCauchy(Fn, ndim, P,Cauchy);
    computeVMS(Cauchy, ndim, VMS);
    // // If Fg(F)
    tensorProductSecondTensorSecondTensor(Fginv, ndim, ndim, delta, ndim, ndim, tanPart122);
    vector<double> tempFA(ndim *ndim, 0.);
    vector<double> Fninv(ndim *ndim, 0.);
    vector<double> FninvT(ndim *ndim, 0.);
    double a = Jn1 * Kp * dt / (rhoMicroReal * theta * (prevTheta + Jn1 * Kp * dt / rhoMicroReal));
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
    // tensorProductSecondTensorSecondTensor(Fginv, ndim, ndim, I, ndim, ndim, tanPart12);

    multFourthTensorFourthTensor(ndim, viscousTanModuli, tanPart12, tanModuli);

    delete viscousTanModuli;
    delete tanPart12;
    delete tanPart122;
    delete tanPart123;

    for (int i = 0; i < _finalSize - 2; i++) {
        intVarsCurr[i] = intVarsCurrVisco[i];
    }
    //    intVarsCurr[_finalSize-3]=intVarsPrev[_finalSize-3]+(theta-prevTheta)/dt;
    intVarsCurr[_finalSize - 2] = (theta - prevTheta) / dt;
    intVarsCurr[_finalSize - 1] = intVarsCurrGrowth[0];
}
