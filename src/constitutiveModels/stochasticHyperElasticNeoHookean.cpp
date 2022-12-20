//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file stochasitcHyperElasticNeoHookean.cpp
  \brief 
*/
#include "stochasticHyperElasticNeoHookean.h"
#include "maths.h"
#include "classGPs.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
*/
classStochasticHyperElasticNeoHookean::classStochasticHyperElasticNeoHookean(int ndim, string name,
                                                                                           string name_cons, double rho,
                                                                                           double Young, double nu,
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
    this->C_Stoch = C;
    this->nstoch = number_stochastic;
    vector<double> lambda(number_stochastic, 0);
    vector<double> mu(number_stochastic, 0);
    vector<double> K1(number_stochastic, 0);
    vector<double> twicemu(number_stochastic, 0);
    vector<double> oneminusnu_carre(number_stochastic, 0);
    vector<double> Young_times_nu(number_stochastic, 0);
    vector<double> fac(number_stochastic, 0);
    vector<double> Young_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> nutimesYoung_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> oneplusnu_times_oneminustwonu(number_stochastic, 0);
    vector<double> oneminustwicenu(number_stochastic, 0);
    vector<double> factor(number_stochastic, 0);
    vector<double> factor_times_nu(number_stochastic, 0);
    vector<double> young_over_oneminustwonu(number_stochastic, 0);
    oneminustwicenu = this->nu;
    vector<double> oneplusnu = this->nu;
    vector<double> _nu = this->nu;
    vector<double> _Young = this->Young; 
    linearCombinationVariableStochastic(oneplusnu, 1, 1);
    linearCombinationVariableStochastic(oneminustwicenu, -2, 1);


    DivisionRandomVariableRandomVariableStochastic(_Young, oneplusnu, C_Stoch, mu);
    vector<double> Young_over_oneplusnu = mu;
    MultiplicationRandomVariableScalar(mu, 0.5);

   
    DivisionRandomVariableRandomVariableStochastic(Young_over_oneplusnu, oneminustwicenu, C_Stoch, factor);
    multiRandomVariableRandomVariableStochastic(factor, _nu, C_Stoch, factor_times_nu);
    lambda = factor_times_nu;
    DivisionRandomVariableRandomVariableStochastic(_Young, oneminustwicenu, C_Stoch, young_over_oneminustwonu);
    K1 = young_over_oneminustwonu;
    MultiplicationRandomVariableScalar(K1, double(1./3));
    this->lambda = lambda;
    this->mu = mu;
    this->K1 = K1;
}

classStochasticHyperElasticNeoHookean::~classStochasticHyperElasticNeoHookean()
{
}

void classStochasticHyperElasticNeoHookean::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
     _term = new MechanicalTermStochastic();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticHyperElasticNeoHookean::soundSpeed() const{
    double rho = this->getRho();
    double Young = this->Young[0];
    //  double factornu = (1.-nu)/((1.+nu)*(1.-2.*nu));
    double sound = sqrt(Young / rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticHyperElasticNeoHookean::initIntVars(vector<double> &intVars) {
    // This constitutive model does not have any internal variable. The arrays is not modified
    // The last position is always for the number of remeshings. MuPhiSim Internal variable
    intVars.push_back(0.0); // This internal variable is initialised to 0
}


void classStochasticHyperElasticNeoHookean::stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr) {
    int number_stochastic = nstoch;
    int ndim = getDim();
    setAll(PK1,0);
    vector<double> FCurrT(ndim *ndim *number_stochastic,0);
    vector<double> FInv(ndim *ndim *number_stochastic, 0);
    vector<double> FInvT(ndim *ndim *number_stochastic, 0);
    vector<double> Jacobian(number_stochastic,0);
    vector<double> Jacobian_minusone(number_stochastic,0);
    vector<double> JJ_minusone(number_stochastic,0);
    vector<double> K1_JJ_minusone(number_stochastic,0);
    vector<double> traceB(number_stochastic,0);
    vector<double> B(ndim *ndim *number_stochastic,0);
    vector<double> factor1(ndim *ndim *number_stochastic,0);
    vector<double> factor2(ndim *ndim *number_stochastic,0);
    vector<double> factor3(ndim *ndim *number_stochastic,0);
    vector<double> firstPart(ndim *ndim *number_stochastic,0);
    vector<double> Jacobian_power(number_stochastic,0);
    vector<double> mu_over_jacobian_power(number_stochastic,0);
    
     

    DeterminantStochasticTensor(FCurr, ndim, number_stochastic, Jacobian, C_Stoch);
    Jacobian_minusone = Jacobian;
    Jacobian_minusone[0] = Jacobian_minusone[0] - 1.0;
    
    InverseTensorStochastic(FCurr, ndim, number_stochastic, C_Stoch, FInv);
    transposeStochastic(FInv, ndim, number_stochastic, FInvT);
    transposeStochastic(FCurr, ndim, number_stochastic, FCurrT);
    
    multSTensor3SecondTransposeStochastic(FCurr, ndim, ndim, FCurr, ndim, B, number_stochastic, C_Stoch);

    for (int m = 0; m < number_stochastic; m++) {
        for (int i = 0; i < ndim; i++) {
            traceB[m] += B[m * ndim * ndim + i * ndim + i];
        }
    }

    multiRandomVariableRandomVariableStochastic(Jacobian, Jacobian_minusone, C_Stoch,
                                                 JJ_minusone);
    multiRandomVariableRandomVariableStochastic(K1, JJ_minusone, C_Stoch,
                                                 K1_JJ_minusone);                                            
    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, K1_JJ_minusone, number_stochastic, C_Stoch,
                                                 factor3);
    multiSTensorRandomVariableStochastic(FInvT, ndim, ndim, traceB, number_stochastic, C_Stoch,
                                                 factor2);
    double power = ((5.0 - double(ndim))/ double(ndim));                                                                                                                                  
    nthPowerIntegral(Jacobian, Jacobian_power, C_Stoch, power);
    for (int i = 0; i < number_stochastic; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                factor1[i*ndim*ndim + j*ndim +k] = FCurr[i*ndim*ndim + j*ndim +k] - (1/(double(ndim)))*factor2[i*ndim*ndim + j*ndim +k];
            }
        }
    }
    
    DivisionRandomVariableRandomVariableStochastic(mu, Jacobian_power, C_Stoch,
                                                   mu_over_jacobian_power);
                                                  
    multiSTensorRandomVariableStochastic(factor1, ndim, ndim, mu_over_jacobian_power, number_stochastic, C_Stoch,
                                                firstPart);
    
    
    for (int i = 0; i < number_stochastic; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                PK1[i*ndim*ndim + j*ndim +k] = firstPart[i*ndim*ndim + j*ndim +k] + factor3[i*ndim*ndim + j*ndim +k];
            }
        }
    }

    multTensorTensor3Stochastic(FInv, ndim, ndim, PK1, ndim, PK2, nstoch, C_Stoch);
    

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticHyperElasticNeoHookean::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    //if(timeRun > 1){
    //    this->Young[0] = 1e7;
    //}
    
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

    stress(PK2, PK1, FCurr);
    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);
    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
    }
        



    computeCauchyStoch(FCurr, ndim, nstoch, C_Stoch, PK1, Cauchy);
    //computeStochasticVMS(Cauchy, ndim, nstoch, VMS, C_Stoch);
    double val =0;
    if (flagTanMod) {
        //To be implemented
    } 

}

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
*/
classStochasticHyperElastic3DCompressibleNeoHookean::classStochasticHyperElastic3DCompressibleNeoHookean(int ndim, string name,
                                                                                           string name_cons, double rho,
                                                                                           double Young, double nu,
                                                                                           vector<double> Stochastic_numbers,
                                                                                           vector<double> Stochastic_parameters,
                                                                                           vector<double> Stochastic_function,
                                                                                           string approximation,
                                                                                           int order, int resolution, bool stochasticMapping, string distribution) : 
                        constitutiveModels(ndim, name, name_cons, rho) {
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
    this->nstoch = number_stochastic;

    vector<double> lambda(number_stochastic, 0);
    vector<double> mu(number_stochastic, 0);
    vector<double> K1(number_stochastic, 0);
    vector<double> twicemu(number_stochastic, 0);
    vector<double> oneminusnu_carre(number_stochastic, 0);
    vector<double> Young_times_nu(number_stochastic, 0);
    vector<double> fac(number_stochastic, 0);
    vector<double> Young_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> nutimesYoung_over_one_minusnu_carre(number_stochastic, 0);
    vector<double> oneplusnu_times_oneminustwonu(number_stochastic, 0);
    vector<double> factor(number_stochastic, 0);
    vector<double> factor_times_nu(number_stochastic, 0);
    vector<double> young_over_oneminustwonu(number_stochastic, 0);
    vector<double> oneminustwicenu = this->nu;
    vector<double> oneplusnu = this->nu;


    linearCombinationVariableStochastic(oneplusnu, 1, 1);
    linearCombinationVariableStochastic(oneminustwicenu, -2, 1);
    DivisionRandomVariableRandomVariableStochastic(this->Young, oneplusnu, C_Stoch, mu);
    vector<double> Young_over_oneplusnu = mu;
    MultiplicationRandomVariableScalar(mu, 0.5);
    
    DivisionRandomVariableRandomVariableStochastic(Young_over_oneplusnu, oneminustwicenu, C_Stoch, factor);
    multiRandomVariableRandomVariableStochastic(factor, this->nu, C_Stoch, factor_times_nu);
    lambda = factor_times_nu;

    DivisionRandomVariableRandomVariableStochastic(this->Young, oneminustwicenu, C_Stoch, young_over_oneminustwonu);
    K1 = young_over_oneminustwonu;
    MultiplicationRandomVariableScalar(K1, 0.33);
    this->lambda = lambda;
    this->mu = mu;
    this->K1 = K1; 
}

classStochasticHyperElastic3DCompressibleNeoHookean::~classStochasticHyperElastic3DCompressibleNeoHookean()
{
   
}

void classStochasticHyperElastic3DCompressibleNeoHookean::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalTermStochastic();;
};

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classStochasticHyperElastic3DCompressibleNeoHookean::soundSpeed() const
{
    return sqrt((lambda[0] + 2 * mu[0]) / _rho);
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classStochasticHyperElastic3DCompressibleNeoHookean::initIntVars(vector<double> &intVars) {
    // This constitutive model does not have internal variables

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classStochasticHyperElastic3DCompressibleNeoHookean::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
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
    vector<double>Jacobian(nstoch,0);
    vector<double>logJacobian(nstoch,0);
    vector<double>invC(nstoch*ndim*ndim,0);
    vector<double>I_minus_invC(nstoch*ndim*ndim,0);
    vector<double>rightCauchy(nstoch*ndim*ndim,0);
    vector<double>firstPart(nstoch*ndim*ndim,0);
    vector<double>secondPart(nstoch*ndim*ndim,0);
    vector<double>lambda_logJ(nstoch,0);
    
    DeterminantStochasticTensor(FCurr, ndim, nstoch, Jacobian, C_Stoch);
    logIntegral(Jacobian, logJacobian, C_Stoch);
    multSTensor3FirstTranspose(FCurr, ndim, ndim, FCurr, ndim, rightCauchy);
    InverseTensorStochastic(rightCauchy, ndim, nstoch, C_Stoch, invC);

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

    for (int i = 0; i < nstoch; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                I_minus_invC[i*ndim*ndim + j*ndim +k] = -invC[i*ndim*ndim + j*ndim +k];
            }
        }
    }

    for(int i=0; i<ndim; i++){
        I_minus_invC[i*ndim +i] = I_minus_invC[i*ndim +i] +1;
    }
    
    multiSTensorRandomVariableStochastic(I_minus_invC, ndim, ndim, mu, nstoch, C_Stoch,
                                                firstPart);
    multiRandomVariableRandomVariableStochastic(lambda, logJacobian, C_Stoch,
                                                 lambda_logJ);                                            
    multiSTensorRandomVariableStochastic(invC, ndim, ndim, lambda_logJ, nstoch, C_Stoch,
                                                secondPart);

    for (int i = 0; i < nstoch; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                PK2[i*ndim*ndim + j*ndim +k] = firstPart[i*ndim*ndim + j*ndim +k] + secondPart[i*ndim*ndim + j*ndim +k];
            }
        }
    }

    vector<double>FIniInvT(nstoch*ndim*ndim,0);
    vector<double>PK2_FIniInvT(nstoch*ndim*ndim,0);
    vector<double>FIniInv_PK2_FIniInvT(nstoch*ndim*ndim,0);
    if(_stochasticMapping){
        transposeStochastic(FIniInv, ndim, nstoch, FIniInvT);
        multTensorTensor3Stochastic(PK2, ndim, ndim, FIniInvT, ndim, PK2_FIniInvT, nstoch, C_Stoch);
        multTensorTensor3Stochastic(FIniInv, ndim, ndim, PK2_FIniInvT, ndim, FIniInv_PK2_FIniInvT, nstoch, C_Stoch);
        multiSTensorRandomVariableStochastic(FIniInv_PK2_FIniInvT, ndim, ndim, det_FIni, nstoch, C_Stoch, PK2);
    }
        



    computeCauchyStoch(FCurr, ndim, nstoch, C_Stoch, PK1, Cauchy);
    //computeStochasticVMS(Cauchy, ndim, nstoch, VMS, C_Stoch);

    if (flagTanMod) {
        //To be computed... Only works with perturbation atm...
    } 

};

double classStochasticHyperElastic3DCompressibleNeoHookean::defoEnergy(classGPs *GaussPoint) const
{
    // W = 0.5*lambda*log(J)**2 - mu*log(J)+ 0.5*mu*trace(FT*F)
    //To be computed
    return 0.0;
}

