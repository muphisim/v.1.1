//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file linearElastic.cpp
  \brief This constitutive model is thought only for testing and just for small deformations. In that scenario, all stresses are equivalent and thus in PK1, the Cauchy stress will be returned
*/
#include "linearElastic.h"
#include "maths.h"
#include "classGPs.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] _Young's modulus
  @param[in] _nu Poisson ratio
  @param[in] stressState 
*/
linearElastic::linearElastic(int ndim, string name, string name_cons, 
                    double rho, double Young, double nu, int stressState) : 
                    constitutiveModels(ndim, name, name_cons, rho),
                    _Young(Young),_nu(nu), _elasticStiffTensor(ndim),
                    _stressState(stressState)
{
    fillElasticStiffTensor(ndim,_elasticStiffTensor);
    
}

linearElastic::~linearElastic()
{
     
}

void linearElastic::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double linearElastic::soundSpeed() const
{
    double rho = this->getRho();
    double sound = sqrt(_Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void linearElastic::initIntVars(vector<double> &intVars) {
    intVars.push_back(0.0);  //activation
}

void linearElastic::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            E[i * ndim + j] = 0.5 * (FCurr[i * ndim + j]+FCurr[j * ndim + i] - 2.*delta[i * ndim + j]);
        }
    }
    stress(PK1, E);
    computeCauchy(FCurr, ndim, PK1, Cauchy);
    computeVMS(Cauchy, ndim, VMS);

    if (flagTanMod) 
    {
        classTensor4 *tanModuli = currState.getTangent();
        (*tanModuli)=(_elasticStiffTensor);
    } 
}

double linearElastic::defoEnergy(classGPs *GaussPoint) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress(); 
    int ndim = getDim(); 
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    return 0.5*(doubleContractionSecondTensors(PK1,FCurr,ndim)-doubleContractionSecondTensors(PK1,delta,ndim));
};

void linearElastic::stress(vector<double> &PK2, const vector<double> &E) const {
    int ndim = getDim();
    int stressState = this->_stressState;
    PK2 = E;
    double mu2 = _Young / (1. + _nu);
    double lam = _Young * _nu / ((1. + _nu) * (1. - (2 * _nu)));
    double traceE = 0;
    for (int i = 0; i < ndim; i++) {
        traceE += E[i * ndim + i];
    }

    if (ndim == 2 && stressState != 1) {
        if (stressState == 0)//plane stress
        {
            double fac = _Young / (1 - _nu * _nu);
            multiSTensorScalar(PK2, ndim, ndim, fac * (1 - _nu));//PK2=E/(1-v^2)*((1-v)E+I*trace(E)*v)
            for (int i = 0; i < ndim; i++) {
                PK2[i * ndim + i] += traceE * fac * _nu;
            }
        } else if (stressState == 3)//isotropic
        {
            ERROR("Not yet implemented");
            exit(-1);
        } else {
            INFO("Invalid stress state!");
            exit(-1);
        }

    } else//plane strain in 2D by definition or isotropic in 3d
    {
        multiSTensorScalar(PK2, ndim, ndim, mu2);
        for (int i = 0; i < ndim; i++) {
            PK2[i * ndim + i] += lam * traceE;    // PK2 = 2*mu*E + lam*tr(E)*I
        }
    }

}

/*! \brief This function defines the stiffness tensor. It will be the tanget moduli as well.
  @param[in] ndim Dimension of the domain
*/
void linearElastic::fillElasticStiffTensor(int ndim, classTensor4 &elasticStiffness) {
    double mu = 0.5 * _Young / (1. + _nu);
    double lambda = (_Young * _nu) / ((1. + _nu) * (1. - 2. * _nu));
    double twicemu = mu + mu;
    elasticStiffness.setAll(0);
    if (ndim == 2) {
        int stressState = this->_stressState;
        if (stressState == 0) {
            INFO("Plane stress");
            elasticStiffness.setValues(0, 0, 0, 0, _Young / (1 - _nu * _nu));
            elasticStiffness.setValues(1, 1, 0, 0, _nu * _Young / (1 - _nu * _nu));
            elasticStiffness.setValues(0, 0, 1, 1, _nu * _Young / (1 - _nu * _nu));
            elasticStiffness.setValues(1, 1, 1, 1, _Young / (1 - _nu * _nu));
            elasticStiffness.setValues(0, 1, 0, 1, (_Young / (1 - _nu * _nu)) * (1 - _nu) / 2);
            elasticStiffness.setValues(1, 0, 1, 0, (_Young / (1 - _nu * _nu)) * (1 - _nu) / 2);
            elasticStiffness.setValues(1, 0, 0, 1, (_Young / (1 - _nu * _nu)) * (1 - _nu) / 2);
            elasticStiffness.setValues(0, 1, 1, 0, (_Young / (1 - _nu * _nu)) * (1 - _nu) / 2);

        } else if (stressState == 1) {
            INFO("Plane strain");

            double factor = _Young / ((1 + _nu) * (1 - 2 * _nu));
            elasticStiffness.setValues(0, 0, 0, 0, factor * (1 - _nu));
            elasticStiffness.setValues(1, 1, 0, 0, factor * _nu);
            elasticStiffness.setValues(0, 0, 1, 1, factor * _nu);
            elasticStiffness.setValues(1, 1, 1, 1, factor * (1 - _nu));
            elasticStiffness.setValues(0, 1, 0, 1, factor * (1 - 2 * _nu) / 2);
            elasticStiffness.setValues(1, 0, 1, 0, factor * (1 - 2 * _nu) / 2);
            elasticStiffness.setValues(1, 0, 0, 1, factor * (1 - 2 * _nu) / 2);
            elasticStiffness.setValues(0, 1, 1, 0, factor * (1 - 2 * _nu) / 2);


        } else {
            INFO("Isotropic");
            elasticStiffness.setValues(0, 0, 0, 0, lambda + twicemu);
            elasticStiffness.setValues(0, 0, 1, 1, lambda);
            elasticStiffness.setValues(1, 1, 0, 0, lambda);
            elasticStiffness.setValues(1, 1, 1, 1, lambda + twicemu);
            elasticStiffness.setValues(0, 1, 0, 1, mu);
            elasticStiffness.setValues(1, 0, 1, 0, mu);
            elasticStiffness.setValues(1, 0, 0, 1, mu);
            elasticStiffness.setValues(0, 1, 1, 0, mu);
        }

    } else {

        INFO("Isotropic 3D");
        // isotropic materials, major and minor symmetries
        elasticStiffness.setValues(0, 0, 0, 0, lambda + twicemu);
        elasticStiffness.setValues(1, 1, 1, 1, lambda + twicemu);
        elasticStiffness.setValues(2, 2, 2, 2, lambda + twicemu);
        elasticStiffness.setValues(0, 0, 1, 1, lambda);
        elasticStiffness.setValues(0, 0, 2, 2, lambda);
        elasticStiffness.setValues(1, 1, 0, 0, lambda);
        elasticStiffness.setValues(1, 1, 2, 2, lambda);
        elasticStiffness.setValues(2, 2, 0, 0, lambda);
        elasticStiffness.setValues(2, 2, 1, 1, lambda);
        elasticStiffness.setValues(1, 2, 1, 2, mu);
        elasticStiffness.setValues(1, 2, 2, 1, mu);
        elasticStiffness.setValues(2, 1, 2, 1, mu);
        elasticStiffness.setValues(2, 1, 1, 2, mu);
        elasticStiffness.setValues(2, 0, 2, 0, mu);
        elasticStiffness.setValues(2, 0, 0, 2, mu);
        elasticStiffness.setValues(0, 2, 0, 2, mu);
        elasticStiffness.setValues(0, 2, 2, 0, mu);
        elasticStiffness.setValues(0, 1, 0, 1, mu);
        elasticStiffness.setValues(0, 1, 1, 0, mu);
        elasticStiffness.setValues(1, 0, 1, 0, mu);
        elasticStiffness.setValues(1, 0, 0, 1, mu);

    }
}


/*! \brief This function actives the elements if the activation time is lower than timeRun
  @param[inout] activation state
*/
void linearElastic::checkActivation(classGPs *GaussPoint, double timeRun) const
{
    const vector<double>& intVarActivation = GaussPoint->getCurrentState().getInternalVariables();
    if (intVarActivation[0] <= timeRun) {
        GaussPoint->setActivate(1);
    } else {
        GaussPoint->setActivate(0);
    }
}
