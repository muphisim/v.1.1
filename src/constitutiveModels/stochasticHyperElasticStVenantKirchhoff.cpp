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
#include "stochasticHyperElasticStVenantKirchhoff.h"
#include "maths.h"
#include "classGPs.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
*/
classStochasticHyperElasticStVenantKirchhoff::classStochasticHyperElasticStVenantKirchhoff(int ndim, string name,
                                                                                           string name_cons, double rho,
                                                                                           double Young, double nu,
                                                                                           int stressState,
                                                                                           vector<double> Stochastic_numbers,
                                                                                           vector<double> Stochastic_parameters,
                                                                                           vector<double> Stochastic_function,
                                                                                           string approximation,
                                                                                           int order, int resolution, bool stochasticMapping, string distribution)
        : constitutiveModels(ndim, name, name_cons, rho) {
    

    int number_func = 0;
    for (int j = 0; j < Stochastic_numbers.size(); j++) {
        number_func = number_func + Stochastic_numbers[j];
    }
    int number_stochastic = 0;
    if (approximation == "Haar") {
        number_stochastic = binomialCoefficients(pow(2, order + 1) + number_func - 1, number_func);
    } else {
        if (approximation == "Mixed") {
        number_stochastic = binomialCoefficients((order + 1)*(2*pow(2,resolution)) + number_func - 1, number_func);
        } else {
            number_stochastic = binomialCoefficients(order + number_func, order);
        }
    }
    this->_stochasticMapping = stochasticMapping;
    vector<double> zeros(number_stochastic, 0);
    this->elasticStiffTensor = new classTensor6(ndim, number_stochastic);
    this->Young = zeros;
    this->nu = zeros;
    this->Young[0] = Young;
    this->nu[0] = nu;
    int offset = 1;
    for (int j = 0; j < Stochastic_parameters.size(); j++) {
        if (Stochastic_parameters[j] == 1) {
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->Young[offset + i] = Young * Stochastic_function[i];
                //count++;
            }
        }
        if (Stochastic_parameters[j] == 2) {
            if(Stochastic_parameters.size()==2){
                offset = offset + Stochastic_numbers[0];
            }
            for (int i = 0; i < Stochastic_numbers[j]; i++) {
                this->nu[offset + i] = nu * Stochastic_function[i + offset - 1];
                offset++;
            }
        }
    }
    this->_stressState = stressState; //0 for plane stress, 1 for plane strain
    vector<double> C(number_stochastic * number_stochastic * number_stochastic, 0);
    if (approximation == "Haar") {
        Build_C_Haar(order, number_func, C);
        ConvertToHaar(this->Young, this->nu, order, Stochastic_parameters, Stochastic_numbers, distribution);
    } else {
        if (approximation == "Mixed"){
            Build_C_Mixed(order, resolution, number_func, C, distribution);
        } else {
            build_ThirdOrderTensor_C(order, number_func, C, distribution);
        }
    }
    C_Stoch = C;
    this->fillElasticStiffTensor(ndim, number_stochastic);
    this->nstoch = number_stochastic;
}

classStochasticHyperElasticStVenantKirchhoff::~classStochasticHyperElasticStVenantKirchhoff()
{
}

void classStochasticHyperElasticStVenantKirchhoff::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
     _term = new MechanicalTermStochastic();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticHyperElasticStVenantKirchhoff::soundSpeed() const{
    double rho = this->getRho();
    double Young = this->Young[0];
    //  double factornu = (1.-nu)/((1.+nu)*(1.-2.*nu));
    double sound = sqrt(Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticHyperElasticStVenantKirchhoff::initIntVars(vector<double> &intVars) {
    // This constitutive model does not have any internal variable. The arrays is not modified
    // The last position is always for the number of remeshings. MuPhiSim Internal variable
    intVars.push_back(0.0); // This internal variable is initialised to 0
}


void classStochasticHyperElasticStVenantKirchhoff::stress(vector<double> &PK2, const vector<double> &E) {
    int number_stochastic = nstoch;
    int ndim = getDim();
    int stressState = this->_stressState;
    PK2 = E;
    vector<double> mu2;
    vector<double> lam(number_stochastic, 0);
    vector<double> traceE;
    vector<double> mu(number_stochastic, 0);
    vector<double> lambda(number_stochastic, 0);
    vector<double> twicemu(number_stochastic, 0);
    vector<double> oneminusnu_carre(number_stochastic, 0);
    vector<double> Young_times_nu(number_stochastic, 0);
    vector<double> fac(number_stochastic, 0);
    vector<double> Young_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> nutimesYoung_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> oneplusnu_times_oneminustwonu(number_stochastic, 0);
    vector<double> oneminusnu_over2(number_stochastic, 0);
    vector<double> oneminusnu(number_stochastic, 0);
    vector<double> oneminustwicenu(number_stochastic, 0);
    oneminustwicenu = nu;
    vector<double> oneminusmu = nu;
    vector<double> oneplusnu = nu;

    oneminusnu_over2 = nu;
    linearCombinationVariableStochastic(oneminusnu_over2, -0.5, 0.5);
    linearCombinationVariableStochastic(oneplusnu, 1, 1);
    linearCombinationVariableStochastic(oneminustwicenu, -2, 1);
    DivisionRandomVariableRandomVariableStochastic(Young, oneplusnu, C_Stoch, mu);
    twicemu = mu;
    mu2 = mu;
    MultiplicationRandomVariableScalar(mu, 0.5);
    multiRandomVariableRandomVariableStochastic(nu, nu, C_Stoch, oneminusnu_carre);
    linearCombinationVariableStochastic(oneminusnu_carre, -1, 1);
    DivisionRandomVariableRandomVariableStochastic(Young, oneplusnu, C_Stoch, fac);
    multiRandomVariableRandomVariableStochastic(Young, nu, C_Stoch, Young_times_nu);
    multiRandomVariableRandomVariableStochastic(oneplusnu, oneminustwicenu, C_Stoch, oneplusnu_times_oneminustwonu);
    DivisionRandomVariableRandomVariableStochastic(Young_times_nu, oneminusnu_carre, C_Stoch,
                                                   nutimesYoung_over_one_minusnu_carre);
    DivisionRandomVariableRandomVariableStochastic(Young_times_nu, oneplusnu_times_oneminustwonu, C_Stoch, lam);

    for (int m = 0; m < nstoch; m++) {
        double traceE_temp = 0;
        for (int i = 0; i < ndim; i++) {
            traceE_temp += E[m * ndim * ndim + i * ndim + i];
        }
        traceE.push_back(traceE_temp);
    }
    if (ndim == 2 && stressState != 1) {
        if (stressState == 0)//plane stress
        {

            vector<double> PK2_temp(ndim *ndim
            *nstoch, 0);
            multiSTensorRandomVariableStochastic(PK2, ndim, ndim, fac, nstoch, C_Stoch,
                                                 PK2_temp);//PK2=E/(1-v^2)*((1-v)E+I*trace(E)*v)
            PK2 = PK2_temp;
            for (int m = 0; m < nstoch; m++) {
                for (int u = 0; u < nstoch; u++) {
                    for (int v = 0; v < nstoch; v++) {
                        for (int i = 0; i < ndim; i++) {
                            PK2[m * ndim * ndim + i * ndim + i] +=
                                    C_Stoch[m * nstoch * nstoch + u * nstoch + v] * traceE[u] *
                                    nutimesYoung_over_one_minusnu_carre[v];
                        }
                    }
                }
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
        vector<double> PK2_temp(ndim *ndim
        *nstoch, 0);
        multiSTensorRandomVariableStochastic(PK2, ndim, ndim, mu2, nstoch, C_Stoch, PK2_temp);
        PK2 = PK2_temp;
        for (int m = 0; m < nstoch; m++) {
            for (int i = 0; i < ndim; i++) {
                for (int u = 0; u < nstoch; u++) {
                    for (int v = 0; v < nstoch; v++) {
                        PK2[m * ndim * ndim + i * ndim + i] += C_Stoch[m * nstoch * nstoch + u * nstoch + v] * lam[u] *
                                                               traceE[v];    // PK2 = 2*mu*E + lam*tr(E)*I

                    }
                }
            }
        }
    }

}

/*! \brief This function defines the stiffness tensor. It will be the tanget moduli as well.
  @param[in] ndim Dimension of the domain
*/
void classStochasticHyperElasticStVenantKirchhoff::fillElasticStiffTensor(int ndim, int number_stochastic) {

    vector<double> mu(number_stochastic, 0);
    vector<double> lambda(number_stochastic, 0);
    vector<double> twicemu(number_stochastic, 0);
    vector<double> oneminusnu_carre(number_stochastic, 0);
    vector<double> Young_times_nu(number_stochastic, 0);
    vector<double> Young_times_one_minusnuover2(number_stochastic, 0);
    vector<double> Young_times_one_minusnuover2_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> Young_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> nutimesYoung_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> oneminusnu_over2(number_stochastic, 0);
    vector<double> oneplusnu = nu;
    vector<double> oneminustwicemu = nu;
    vector<double> oneminusnu = nu;
    vector<double> Young_over_oneplusnu(number_stochastic, 0);
    vector<double> factor(number_stochastic, 0);
    vector<double> factor_times_nu(number_stochastic, 0);
    vector<double> factor_times_oneminusnu(number_stochastic, 0);
    oneminusnu_over2 = nu;
    linearCombinationVariableStochastic(oneminusnu_over2, -0.5, 0.5);
    linearCombinationVariableStochastic(oneplusnu, 1, 1);
    linearCombinationVariableStochastic(oneminustwicemu, -2, 1);
    linearCombinationVariableStochastic(oneminusnu, -1, 1);
    DivisionRandomVariableRandomVariableStochastic(Young, oneplusnu, C_Stoch, mu);
    twicemu = mu;
    Young_over_oneplusnu = mu;
    MultiplicationRandomVariableScalar(mu, 0.5);
    DivisionRandomVariableRandomVariableStochastic(Young_over_oneplusnu, oneminustwicemu, C_Stoch, factor);
    multiRandomVariableRandomVariableStochastic(factor, nu, C_Stoch, factor_times_nu);
    lambda = factor_times_nu;
    multiRandomVariableRandomVariableStochastic(factor, oneminusnu, C_Stoch, factor_times_oneminusnu);
    multiRandomVariableRandomVariableStochastic(nu, nu, C_Stoch, oneminusnu_carre);
    linearCombinationVariableStochastic(oneminusnu_carre, -1, 1);
    DivisionRandomVariableRandomVariableStochastic(Young, oneminusnu_carre, C_Stoch, Young_over_one_minusnu_carre);
    multiRandomVariableRandomVariableStochastic(Young, oneminusnu_over2, C_Stoch, Young_times_one_minusnuover2);
    multiRandomVariableRandomVariableStochastic(Young, nu, C_Stoch, Young_times_nu);
    DivisionRandomVariableRandomVariableStochastic(Young_times_nu, oneminusnu_carre, C_Stoch,
                                                   nutimesYoung_over_one_minusnu_carre);
    DivisionRandomVariableRandomVariableStochastic(Young_times_one_minusnuover2, oneminusnu_carre, C_Stoch,
                                                   Young_times_one_minusnuover2_over_one_minusnu_carre);

        if (ndim == 2) {
            int stressState = this->_stressState;
            if (stressState == 0) {
                INFO("Plane stress");
                for(int beta=0; beta<number_stochastic; beta++) {
                    for (int psi = 0; psi < number_stochastic; psi++) {
                        vector<double> temp_young_over_oneminusnu_carre(number_stochastic * number_stochastic, 0);
                        vector<double> temp_nu_times_young_over_oneminusnu_carre(number_stochastic * number_stochastic, 0);
                        vector<double> temp_Young_times_one_minusnuover2_over_one_minusnu_carre(number_stochastic * number_stochastic, 0);
                        for (int gamma = 0; gamma < number_stochastic; gamma++) {
                            temp_young_over_oneminusnu_carre[psi + beta * number_stochastic] += C_Stoch[gamma + psi * number_stochastic +
                                                                               beta * number_stochastic *
                                                                               number_stochastic] * Young_over_one_minusnu_carre[gamma];
                            temp_nu_times_young_over_oneminusnu_carre[psi + beta * number_stochastic] += C_Stoch[gamma + psi * number_stochastic +
                                                                                   beta * number_stochastic *
                                                                                   number_stochastic] * nutimesYoung_over_one_minusnu_carre[gamma];
                            temp_Young_times_one_minusnuover2_over_one_minusnu_carre[psi + beta * number_stochastic] += C_Stoch[gamma + psi * number_stochastic +
                                                                                                                 beta * number_stochastic *
                                                                                                                 number_stochastic] * Young_times_one_minusnuover2_over_one_minusnu_carre[gamma];
                        }
                        this->elasticStiffTensor->setValues(beta, psi, 0, 0, 0, 0,
                                                            temp_young_over_oneminusnu_carre[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 1, 1, 1,
                                                            temp_young_over_oneminusnu_carre[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 0, 0, 1, 1,
                                                            temp_nu_times_young_over_oneminusnu_carre[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 1, 0, 0,
                                                            temp_nu_times_young_over_oneminusnu_carre[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 0, 1, 0, 1,
                                                            temp_Young_times_one_minusnuover2_over_one_minusnu_carre[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 0, 1, 1, 0,
                                                            temp_Young_times_one_minusnuover2_over_one_minusnu_carre[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 0, 1, 0,
                                                            temp_Young_times_one_minusnuover2_over_one_minusnu_carre[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 0, 0, 1,
                                                            temp_Young_times_one_minusnuover2_over_one_minusnu_carre[psi + beta * number_stochastic]);
                    }
                }
            } else if (stressState == 1) {
                INFO("Plane strain");


            } else {
                INFO("Isotropic");
                for(int beta=0; beta<number_stochastic; beta++) {
                    for (int psi = 0; psi < number_stochastic; psi++) {
                        vector<double> temp_mu(number_stochastic * number_stochastic, 0);
                        vector<double> temp_lambda(number_stochastic * number_stochastic, 0);
                        for (int gamma = 0; gamma < number_stochastic; gamma++) {
                            temp_mu[psi + beta * number_stochastic] += C_Stoch[gamma + psi * number_stochastic +
                                                                               beta * number_stochastic *
                                                                               number_stochastic] * mu[gamma];
                            temp_lambda[psi + beta * number_stochastic] += C_Stoch[gamma + psi * number_stochastic +
                                                                                   beta * number_stochastic *
                                                                                   number_stochastic] * lambda[gamma];
                        }
                        this->elasticStiffTensor->setValues(beta, psi, 0, 0, 0, 0,
                                                            temp_lambda[psi + beta * number_stochastic] +
                                                            2 * temp_mu[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 1, 1, 1,
                                                            temp_lambda[psi + beta * number_stochastic] +
                                                            2 * temp_mu[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 0, 0, 1, 1,
                                                            temp_lambda[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 1, 0, 0,
                                                            temp_lambda[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 0, 1, 0, 1,
                                                            temp_mu[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 0, 1, 1, 0,
                                                            temp_mu[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 0, 1, 0,
                                                            temp_mu[psi + beta * number_stochastic]);
                        this->elasticStiffTensor->setValues(beta, psi, 1, 0, 0, 1,
                                                            temp_mu[psi + beta * number_stochastic]);
                    }
                }

            }

        } else {

            INFO("Isotropic 3D");
            // isotropic materials, major and minor symmetries

        for(int beta=0; beta<number_stochastic; beta++){
            for(int psi=0; psi<number_stochastic; psi++){
                vector<double> temp_mu(number_stochastic*number_stochastic, 0);
                vector<double> temp_lambda(number_stochastic*number_stochastic, 0);
                for(int gamma=0; gamma<number_stochastic; gamma++){
                    temp_mu[psi + beta *number_stochastic] += C_Stoch[gamma + psi*number_stochastic + beta *number_stochastic*number_stochastic]*mu[gamma];
                    temp_lambda[psi + beta *number_stochastic] += C_Stoch[gamma + psi*number_stochastic + beta *number_stochastic*number_stochastic]*lambda[gamma];
                }
                this->elasticStiffTensor->setValues(beta, psi,0, 0, 0, 0, temp_lambda[psi + beta *number_stochastic] + 2*temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,1, 1, 1, 1, temp_lambda[psi + beta *number_stochastic] + 2*temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,2, 2, 2, 2, temp_lambda[psi + beta *number_stochastic] + 2*temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 0, 0, 1, 1, temp_lambda[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 0, 0, 2, 2, temp_lambda[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 1, 1, 0, 0, temp_lambda[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 1, 1, 2, 2, temp_lambda[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 2, 2, 0, 0, temp_lambda[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 2, 2, 1, 1, temp_lambda[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 1, 2, 1, 2, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 1, 2, 2, 1, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,2, 1, 2, 1, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi, 2, 1, 1, 2, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,2, 0, 2, 0, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,2, 0, 0, 2, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,0, 2, 0, 2, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,0, 2, 2, 0, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,0, 1, 0, 1, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,0, 1, 1, 0, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,1, 0, 1, 0, temp_mu[psi + beta *number_stochastic]);
                this->elasticStiffTensor->setValues(beta, psi,1, 0, 0, 1, temp_mu[psi + beta *number_stochastic]);
            }

        }
    }
}

/*! \brief This function defines the product of C and F summed on one stochastic variable
  @param[in/out] FCurr Tensor of order 4 C_{beta,iota,K,L}
  @param[in] FCurr Deformation tensor
*/
void classStochasticHyperElasticStVenantKirchhoff::build_CF(vector<double> &CF, vector<double> &FCurr) {
    int ndim = getDim();
    for(int alpha=0; alpha<nstoch; alpha++){
        for(int beta=0; beta<nstoch; beta++){
            for(int iota=0; iota<nstoch; iota++){
                for(int i=0; i<ndim; i++){
                    for(int K=0; K<ndim; K++){
                        CF[beta*ndim*ndim*nstoch + iota*ndim*ndim + i*ndim + K] += C_Stoch[iota*nstoch*nstoch + beta*nstoch + alpha]*FCurr[alpha*ndim*ndim+ i*ndim + K];
                    }
                }
            }
        }
    }
}


/*! \brief This function defines the product of C and F summed on one stochastic variable
  @param[in/out] FCurr Tensor of order 4 C_{beta,iota,K,L}
  @param[in] FCurr Deformation tensor
*/
void classStochasticHyperElasticStVenantKirchhoff::build_dEdF(vector<double> &dEdF, vector<double> &FCurr) {
    int ndim = getDim();
    vector<double> C_times_F(nstoch * nstoch * ndim * ndim);
    for (int psi = 0; psi < nstoch; psi++) {
        for (int eps = 0; eps < nstoch; eps++) {
            for (int phi = 0; phi < nstoch; phi++) {
                for (int k = 0; k < ndim; k++) {
                    for (int S = 0; S < ndim; S++) {
                        C_times_F[psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + S] +=
                                C_Stoch[psi * nstoch * nstoch + phi * nstoch + eps] * FCurr[phi * ndim * ndim + k * ndim + S];

                    }
                }
            }
        }
    }
    if (ndim==3) {
        for (int psi = 0; psi < nstoch; psi++) {
            for (int eps = 0; eps < nstoch; eps++) {
                for (int k = 0; k < ndim; k++) {
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                         ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 0] = C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                         ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 0] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 2 +
                         ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 0] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 2];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                         ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 0] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                         ndim * ndim * nstoch * 2 + ndim * ndim * eps + ndim * k + 0] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 2];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                         ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 1] = C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                         ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 1] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                         ndim * ndim * nstoch * 2 + ndim * ndim * eps + ndim * k + 1] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 2];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                         ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 1] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 2 +
                         ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 1] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 2];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 2 +
                         ndim * ndim * nstoch * 2 + ndim * ndim * eps + ndim * k + 2] = C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 2];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                         ndim * ndim * nstoch * 2 + ndim * ndim * eps + ndim * k + 2] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                         ndim * ndim * nstoch * 2 + ndim * ndim * eps + ndim * k + 2] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 2 +
                         ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 2] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                    dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 2 +
                         ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 2] = 0.5 * C_times_F[
                            psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                }
            }
        }
    }else {
        if(ndim==2){
            for (int psi = 0; psi < nstoch; psi++) {
                for (int eps = 0; eps < nstoch; eps++) {
                    for (int k = 0; k < ndim; k++) {
                        dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                             ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 0] = C_times_F[
                                psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                        dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                             ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 0] = 0.5 * C_times_F[
                                psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                        dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                             ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 0] = 0.5 * C_times_F[
                                psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                        dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                             ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 1] = C_times_F[
                                psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 1];
                        dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 1 +
                             ndim * ndim * nstoch * 0 + ndim * ndim * eps + ndim * k + 1] = 0.5 * C_times_F[
                                psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                        dEdF[ndim * ndim * ndim * ndim * nstoch * psi + ndim * ndim * ndim * nstoch * 0 +
                             ndim * ndim * nstoch * 1 + ndim * ndim * eps + ndim * k + 1] = 0.5 * C_times_F[
                                psi * nstoch * ndim * ndim + eps * ndim * ndim + k * ndim + 0];
                    }
                }
            }
        }
    }

}

/*! \brief This function defines the product of C and F summed on one stochastic variable
  @param[in/out] dSdF Tensor of order 4 C_{beta,iota,K,L}
  @param[in] FCurr Deformation tensor
*/
void classStochasticHyperElasticStVenantKirchhoff::build_dSdF(vector<double> &dSdF,vector<double> &dEdF) {
    int ndim = getDim();
    for(int beta=0;beta<nstoch;beta++) {
        for (int i = 0; i < ndim; i++) {
            for (int J = 0; J < ndim; J++) {
                for (int eps = 0; eps < nstoch; eps++) {
                    for (int k = 0; k < ndim; k++) {
                        for (int M = 0; M < ndim; M++) {
                            for (int psi = 0; psi < nstoch; psi++) {
                                for (int R = 0; R < ndim; R++) {
                                    for (int S = 0; S < ndim; S++) {
                                           dSdF [beta*ndim*ndim*ndim*ndim*nstoch+ i*ndim*ndim*ndim*nstoch + J*ndim*ndim*nstoch+ eps*ndim*ndim +k*ndim + M] +=
                                                   elasticStiffTensor->get(beta, psi, i, J, R, S)*dEdF[psi*ndim*ndim*ndim*ndim*nstoch + R*ndim*ndim*ndim*nstoch + S*ndim*ndim*nstoch + eps*ndim*ndim +  k*ndim + M];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticHyperElasticStVenantKirchhoff::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    if(timeRun > 1){
        this->Young[0] = 1e7;
    }
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    classStates& initState = GaussPoint->getInitialState();
    vector<double> &FIni = initState.getDeformationGradient();
    vector<double> &FCurr = currState.getDeformationGradient();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    int ndim = getDim();
    vector<double>FIniInv(nstoch*ndim*ndim,0);
    vector<double>FCurrMecha(nstoch*ndim*ndim,0);
    vector<double>det_FIni(nstoch,0);
    if(_stochasticMapping) {
        InverseTensorStochastic(FIni, ndim, nstoch, C_Stoch, FIniInv);
        multTensorTensor3Stochastic(FCurr, ndim, ndim, FIniInv, ndim, FCurrMecha, nstoch, C_Stoch);
        DeterminantStochasticTensor(FIni, ndim, nstoch, det_FIni, C_Stoch);
    }else {
            FCurrMecha = FCurr;
    }

    computeGLstoch(FCurrMecha, ndim, nstoch, C_Stoch, E);
    //computeEquivalentStrainStochastic(E, ndim, nstoch, C_Stoch, equivalentStrain);
    computeVolumetricStrainStochastic(E, ndim, nstoch, volumetricStrain);
    vector<double> PK2(ndim *ndim *nstoch, 0);
    vector<double> Finv(ndim *ndim *nstoch, 0);
    vector<double> CF(nstoch*nstoch*ndim*ndim,0);
    vector<double> dEdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);
    vector<double> dSdF(nstoch*nstoch*ndim*ndim*ndim*ndim,0);


    stress(PK2, E);
    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);
    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
    }
        multTensorTensor3Stochastic(FCurr, ndim, ndim, PK2, ndim, PK1, nstoch, C_Stoch);




    double val=0;
    computeCauchyStoch(FCurr, ndim, nstoch, C_Stoch, PK1, Cauchy);
    //computeStochasticVMS(Cauchy, ndim, nstoch, VMS, C_Stoch);
    build_CF(CF, FCurr);
    build_dEdF(dEdF, FCurr);
    build_dSdF(dSdF, dEdF);
    if (flagTanMod) {
        classTensor6* tanModuli = currState.getTangentStochastic();
        for (int iota = 0; iota < nstoch; iota++) {
            for (int i = 0; i < ndim; i++) {
                for (int J = 0; J < ndim; J++) {
                    for (int eps = 0; eps < nstoch; eps++) {
                        for (int k = 0; k < ndim; k++) {
                            for (int M = 0; M < ndim; M++) {
                                val =0;
                                 if(i==k){
                                     for(int beta =0; beta < nstoch; beta++){
                                         val += C_Stoch[iota*nstoch*nstoch + beta*nstoch + eps]*PK2[beta*ndim*ndim + M*ndim +J];
                                     }
                                 }
                                for(int beta =0; beta < nstoch; beta++) {
                                    for (int K = 0; K < ndim; K++) {
                                        val += CF[beta*ndim*ndim*nstoch + iota*ndim*ndim +i*ndim + K]*dSdF[beta*ndim*ndim*ndim*ndim*nstoch + K*ndim*ndim*ndim*nstoch + J*ndim*ndim*nstoch + eps*ndim*ndim + ndim*k +M];
                                    }
                                }
                                tanModuli->setValues(iota, eps, i, J, k, M, val);

                            }
                        }
                    }
                }
            }
        }
    } 
};
