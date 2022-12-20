//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file hyperElasticStVenantKirchhoff.cpp
  \brief This constitutive model is thought only for testing and just for small deformations. In that scenario, all stresses are equivalent and thus in PK1, the Cauchy stress will be returned
*/
#include "hyperElasticStVenantKirchhoff.h"
#include "maths.h"
#include "classGPs.h"
#include "TensorOperations.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young's modulus
  @param[in] nu Poisson ratio
  @param[in] stressState 
*/
classHyperElasticStVenantKirchhoff::classHyperElasticStVenantKirchhoff(int ndim, string name, string name_cons,
                                                                       double rho, double Young, double nu,
                                                                       int stressState) : 
                    constitutiveModels(ndim, name, name_cons, rho),
                    _elasticStiffTensor(ndim),_Young(Young),_nu(nu),
                    _stressState (stressState)
                    
{
    this->fillElasticStiffTensor(ndim,_elasticStiffTensor);
}

classHyperElasticStVenantKirchhoff::~classHyperElasticStVenantKirchhoff()
{
}

void classHyperElasticStVenantKirchhoff::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classHyperElasticStVenantKirchhoff::soundSpeed()  const
{
//  double factornu = (1.-_nu)/((1.+_nu)*(1.-2.*_nu));
    double sound = sqrt(_Young / _rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classHyperElasticStVenantKirchhoff::initIntVars(vector<double> &intVars) {
    intVars.push_back(0.0);  //activation
}

void classHyperElasticStVenantKirchhoff::stress(vector<double> &PK2, const vector<double> &E) const
{
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
void classHyperElasticStVenantKirchhoff::fillElasticStiffTensor(int ndim, classTensor4& elasticStiffTensor) const {
    double Young = _Young;
    double nu = _nu;
    double mu = 0.5 * Young / (1. + nu);
    double lambda = (Young * nu) / ((1. + nu) * (1. - 2. * nu));
    double twicemu = mu + mu;
    elasticStiffTensor.setAll(0);
    if (ndim == 2) {
        int stressState = this->_stressState;
        if (stressState == 0) {
            INFO("Plane stress");
            elasticStiffTensor.setValues(0, 0, 0, 0, Young / (1 - nu * nu));
            elasticStiffTensor.setValues(1, 1, 0, 0, nu * Young / (1 - nu * nu));
            elasticStiffTensor.setValues(0, 0, 1, 1, nu * Young / (1 - nu * nu));
            elasticStiffTensor.setValues(1, 1, 1, 1, Young / (1 - nu * nu));
            elasticStiffTensor.setValues(0, 1, 0, 1, (Young / (1 - nu * nu)) * (1 - nu) / 2);
            elasticStiffTensor.setValues(1, 0, 1, 0, (Young / (1 - nu * nu)) * (1 - nu) / 2);
            elasticStiffTensor.setValues(1, 0, 0, 1, (Young / (1 - nu * nu)) * (1 - nu) / 2);
            elasticStiffTensor.setValues(0, 1, 1, 0, (Young / (1 - nu * nu)) * (1 - nu) / 2);

        } else if (stressState == 1) {
            INFO("Plane strain");

            double factor = Young / ((1 + nu) * (1 - 2 * nu));
            elasticStiffTensor.setValues(0, 0, 0, 0, factor * (1 - nu));
            elasticStiffTensor.setValues(1, 1, 0, 0, factor * nu);
            elasticStiffTensor.setValues(0, 0, 1, 1, factor * nu);
            elasticStiffTensor.setValues(1, 1, 1, 1, factor * (1 - nu));
            elasticStiffTensor.setValues(0, 1, 0, 1, factor * (1 - 2 * nu) / 2);
            elasticStiffTensor.setValues(1, 0, 1, 0, factor * (1 - 2 * nu) / 2);
            elasticStiffTensor.setValues(1, 0, 0, 1, factor * (1 - 2 * nu) / 2);
            elasticStiffTensor.setValues(0, 1, 1, 0, factor * (1 - 2 * nu) / 2);


        } else {
            INFO("Isotropic");
            elasticStiffTensor.setValues(0, 0, 0, 0, lambda + twicemu);
            elasticStiffTensor.setValues(0, 0, 1, 1, lambda);
            elasticStiffTensor.setValues(1, 1, 0, 0, lambda);
            elasticStiffTensor.setValues(1, 1, 1, 1, lambda + twicemu);
            elasticStiffTensor.setValues(0, 1, 0, 1, mu);
            elasticStiffTensor.setValues(1, 0, 1, 0, mu);
            elasticStiffTensor.setValues(1, 0, 0, 1, mu);
            elasticStiffTensor.setValues(0, 1, 1, 0, mu);
        }

    } else {

        INFO("Isotropic 3D");
        // isotropic materials, major and minor symmetries
        elasticStiffTensor.setValues(0, 0, 0, 0, lambda + twicemu);
        elasticStiffTensor.setValues(1, 1, 1, 1, lambda + twicemu);
        elasticStiffTensor.setValues(2, 2, 2, 2, lambda + twicemu);
        elasticStiffTensor.setValues(0, 0, 1, 1, lambda);
        elasticStiffTensor.setValues(0, 0, 2, 2, lambda);
        elasticStiffTensor.setValues(1, 1, 0, 0, lambda);
        elasticStiffTensor.setValues(1, 1, 2, 2, lambda);
        elasticStiffTensor.setValues(2, 2, 0, 0, lambda);
        elasticStiffTensor.setValues(2, 2, 1, 1, lambda);
        elasticStiffTensor.setValues(1, 2, 1, 2, mu);
        elasticStiffTensor.setValues(1, 2, 2, 1, mu);
        elasticStiffTensor.setValues(2, 1, 2, 1, mu);
        elasticStiffTensor.setValues(2, 1, 1, 2, mu);
        elasticStiffTensor.setValues(2, 0, 2, 0, mu);
        elasticStiffTensor.setValues(2, 0, 0, 2, mu);
        elasticStiffTensor.setValues(0, 2, 0, 2, mu);
        elasticStiffTensor.setValues(0, 2, 2, 0, mu);
        elasticStiffTensor.setValues(0, 1, 0, 1, mu);
        elasticStiffTensor.setValues(0, 1, 1, 0, mu);
        elasticStiffTensor.setValues(1, 0, 1, 0, mu);
        elasticStiffTensor.setValues(1, 0, 0, 1, mu);

    }
}


/*! \brief This function actives the elements if the activation time is lower than timeRun
  @param[inout] activation state
*/
 void classHyperElasticStVenantKirchhoff::checkActivation(classGPs *GaussPoint, double timeRun) const
 {
    const vector<double>& intVarActivation = GaussPoint->getCurrentState().getInternalVariables();
    if (intVarActivation[0] <= timeRun) 
    {
        GaussPoint->setActivate(1);
    } else {
        GaussPoint->setActivate(0);
    }
}


/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classHyperElasticStVenantKirchhoff::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    int ndim = getDim();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();

    vector<double> PK2(ndim *ndim, 0);
    vector<double> Finv(ndim *ndim, 0);

    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    stress(PK2, E);
    multTensorTensor3(FCurr, ndim, ndim, PK2, ndim, PK1);
    computeCauchy(FCurr, ndim, PK1, Cauchy);
    computeVMS(Cauchy, ndim, VMS);
    if (flagTanMod) {
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
                            val += PK2[J * ndim + L];
                        }
                        tanModuli->setValues(i, J, k, L, val);
                    }
                }
            }
        }
    }
};

double classHyperElasticStVenantKirchhoff::defoEnergy(classGPs *GaussPoint) const
{
    const classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress(); 
    int ndim = getDim(); 
    vector<double> E(ndim *ndim, 0);
    vector<double> PK2(ndim *ndim, 0);
    vector<double> rightCauchy(ndim *ndim, 0);
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    multSTensor3FirstTranspose(FCurr, ndim, ndim, FCurr, ndim, rightCauchy);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            E[i * ndim + j] = 0.5 * (rightCauchy[i * ndim + j] - delta[i * ndim + j]);
        }
    }
    stress(PK2, E);
    return 0.5*(doubleContractionSecondTensors(PK2,E,ndim));
};
