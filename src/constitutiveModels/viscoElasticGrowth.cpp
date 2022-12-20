//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file viscoElasticGrowth.cpp
  \brief This file contains all functions related to viscoElastic growth constitutive model. It is a wrapper to couple a growth model with a viscoElastic model. The growth model should be defined in this file. The viscoElastic part is the original formulation proposed by Simo et al (nonlinear viscoelasticity, Simo J.C., Computational Inelasticity, Springer). There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/

#include "viscoElasticGrowth.h"
#include "maths.h"
#include "viscoElastic.h"
#include "isoMorphogeneticGrowth.h"
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
*/
classViscoElasticGrowth::classViscoElasticGrowth(int ndim, string name, string name_cons, double rho, double K, int N,
                                                 vector<double> mus, vector<double> etas, double Gc)
        : constitutiveModels(ndim, name, name_cons, rho) 
{
    this->viscoElastic = new classViscoElastic(ndim, name, name_cons, rho, K, N, mus, etas);
    vector<double> _mu = mus;
    vector<double> _eta = etas;
    double _K = K;
    double sum_mu = 0.;
    for (int i = 0; i < _mu.size(); ++i) {
        sum_mu += _mu[i];
    }
    _Young = 9. * _K * sum_mu / (3. * _K + sum_mu);
    _nu = (3. * _K - 2. * sum_mu) / 2. / (3. * _K + sum_mu);
    //this->growth=new classFiberGrowth(ndim, "FiberGrowthModel", rho, _Young, _nu, Gc);
    this->growth = new classIsoMorphogeneticGrowth(ndim, name, name_cons, rho, _Young, _nu, Gc);
    this->_N = N;
}

classViscoElasticGrowth::~classViscoElasticGrowth()
{
    delete viscoElastic;
    delete growth;
}

void classViscoElasticGrowth::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classViscoElasticGrowth::soundSpeed() const
{
    double rho = this->getRho();
    double sound = sqrt(_Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classViscoElasticGrowth::initIntVars(vector<double> &intVars) {
    // This constitutive model has as many internal variables as maxwell branches. Then, for each branch there are two tensors Hn and Sn. Then the internal variables array will be built as follows:
    /*
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
    // and just at the end, I add the Gc int var of the growth part
    _finalSize = _finalSize + 1;
    intVars.push_back(1.0); // This internal variable is initialised to 1

};


/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classViscoElasticGrowth::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    if (flagTanMod) {
        ERROR("Tangent modulus is not available for the viscoelastic growth model");
        exit(0);
    }
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
    
    vector<double> Fa(ndim *ndim, 0.); // Contractille
    vector<double> Fainv(ndim *ndim, 0.);
    vector<double> Fg(ndim *ndim, 0.);
    vector<double> Fv(ndim *ndim, 0.);
    vector<double> Fginv(ndim *ndim, 0.);
    vector<double> intVarsPrevGrowth(1, 0.), intVarsCurrGrowth(1, 0.);
    vector<double> intVarsPrevVisco(_finalSize - 1, 0.), intVarsCurrVisco(_finalSize - 1, 0.);
    // The last position of the whole array finalSize
    intVarsPrevGrowth[0] = intVarsPrev[_finalSize - 1];

    for (int i = 0; i < _finalSize - 1; i++) {
        intVarsPrevVisco[i] = intVarsPrev[i];
        intVarsCurrVisco[i] = intVarsCurr[i];
    }
    //  static_cast<classFiberGrowth*>(this->growth)->getFg(Fn, F0, intVarsPrevGrowth, intVarsCurrGrowth, dt, timeRun, flagTanMod, P,  tanModuli, Fg);// It does not depend on anything, just on the time (linear growth)
    // static_cast<classAreaGrowth*>(this->growth)->getFg(Fn, F0, intVarsPrevGrowth, intVarsCurrGrowth, dt, timeRun, flagTanMod, P,  tanModuli, Fg);// It does not depend on anything, just on the time (linear growth)
    static_cast<classIsoMorphogeneticGrowth *>(this->growth)->getFg(Fn, F0, intVarsPrevGrowth, intVarsCurrGrowth, dt,
                                                                    timeRun, flagTanMod, P, tanModuli,
                                                                    Fg);// It does not depend on anything, just on the time (linear growth)
    inverse(Fg, ndim, Fginv);
    multTensorTensor3(Fn, ndim, ndim, Fginv, ndim, Fv);
    static_cast<classViscoElastic *>(this->viscoElastic)->updateConstitutive(Fv, F0, intVarsPrevVisco, intVarsCurrVisco,
                                                                             dt, timeRun, flagTanMod, P, P, tanModuli);

    for (int i = 0; i < _finalSize - 1; i++) {
        intVarsCurr[i] = intVarsCurrVisco[i];
    }
    intVarsCurr[_finalSize - 1] = intVarsCurrGrowth[0];
};
