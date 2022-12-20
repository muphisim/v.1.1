//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file 3DprintingThermoelasticityStVenantKirchhoff.cpp
  \brief This constitutive model based on the multiplicative decomposition of F=F^{theta}F^{e} where the elastic stress is defined through St Venant-Kirchhoff model.
   A lineal dependency of the Young's modulus with the temperature is assumed.
*/
#include "3DprintingThermoelasticityStVenantKirchhoff.h"
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
  @param[in] alpha Thermal expansion coefficients
*/

class3DprintingThermoelasticityStVenantKirchhoffBase::class3DprintingThermoelasticityStVenantKirchhoffBase(int ndim,
                                                                                                   string name,
                                                                                                   string name_cons,
                                                                                                   double rho,
                                                                                                   int fieldIndex):
         constitutiveModels(ndim, name, name_cons, rho), _tempratureFieldIndex(fieldIndex)
 {
     
 }           
                                                                                                   
 
/*! \brief This function initialises this law from options
*/
void class3DprintingThermoelasticityStVenantKirchhoffBase::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalExtraDofFullCouplingTerm(std::vector<int>(1,_tempratureFieldIndex));
};


/*! \brief This function actives the elements if the activation time is lower than timeRun
  @param[inout] activation state
*/
void class3DprintingThermoelasticityStVenantKirchhoffBase::checkActivation(classGPs *GaussPoint, double timeRun) const
{
    vector<double>& intVarActivation = GaussPoint->getCurrentState().getInternalVariables();
    vector<double>& intVarActivationPrev = GaussPoint->getPreviousState().getInternalVariables();
    if (intVarActivation[0] <= timeRun) 
    {
        GaussPoint->setActivate(1);
        intVarActivation[_intVarsSize+1] = 1;
        intVarActivationPrev[_intVarsSize+1] = 1;
    } else {
        GaussPoint->setActivate(0);
        intVarActivation[_intVarsSize+1] = 0;
        intVarActivationPrev[_intVarsSize+1] = 0;
    };
};                               

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void class3DprintingThermoelasticityStVenantKirchhoffBase::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    int ndim = getDim();
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    double CurrTemperature = currState.getExtraFields()[_tempratureFieldIndex];
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    //
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();

    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);

    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    
    double& mecaSource = currState.getExtraFieldSources()[_tempratureFieldIndex];    
    //
    double phaseIndex = intVarsCurr[_intVarsSize + 2];
    classTensor2 Ftensor(ndim, FCurr);
    classTensor2 Fini(ndim, GaussPoint->getInitialState().getDeformationGradient());
    // compute rateF
    classTensor2 rateF(ndim,FPrev);
    rateF.scale(-1./dt);
    rateF.axpy(Ftensor,1./dt);
    if (flagTanMod)
    {
        classTensor4 *tanModuli = currState.getTangent();
        
        static classTensor2 PK1tensor(ndim), DPK1tensorDT(ndim);
        static classTensor4 ddPdTdF(ndim);
        static classTensor2 ddPdTdT(ndim);
        compute(GaussPoint, phaseIndex, Fini, Ftensor, CurrTemperature,  PK1tensor, DPK1tensorDT,
                intVarsCurr, intVarsPrev,
                true, tanModuli, &ddPdTdF, &ddPdTdT);
        
        PK1 = PK1tensor.getData();
        mecaSource -= CurrTemperature*TensorOperations::doubleDot(DPK1tensorDT,rateF);
        
        // derivatives
        classTensor2& dPdtheta =  (*currState.getDfirstPiolaKirchhoffDextraFields())[_tempratureFieldIndex];
        dPdtheta = DPK1tensorDT;
        
        classTensor2& DsourceDF = (*currState.getDextraFieldSourcesDdeformationGradient())[_tempratureFieldIndex];
        double& DsourceDT = (*currState.getDextraFieldSourcesDextraFields())[_tempratureFieldIndex][_tempratureFieldIndex];
        
        DsourceDT -= TensorOperations::doubleDot(DPK1tensorDT,rateF)+CurrTemperature*TensorOperations::doubleDot(ddPdTdT,rateF);;
        for (int i = 0; i < ndim; i++) 
        {
            for (int J = 0; J < ndim; J++) 
            {
                DsourceDF(i,J) -= CurrTemperature*DPK1tensorDT(i,J)/dt;
                for (int s = 0; s < ndim; s++) 
                {
                    for (int T = 0; T < ndim; T++) 
                    {
                        DsourceDF(i,J) -= CurrTemperature*ddPdTdF(s,T,i,J)*rateF(s,T);
                    }
                }
            }
        }
    }
    else
    {
        classTensor2 PK1tensor(ndim), DPK1tensorDT(ndim);
        compute(GaussPoint, phaseIndex, Fini, Ftensor, CurrTemperature,  PK1tensor, DPK1tensorDT, intVarsCurr,intVarsPrev);
        PK1 = PK1tensor.getData();
        mecaSource -= CurrTemperature*TensorOperations::doubleDot(DPK1tensorDT,rateF);
        
    };    
};         
                                                            

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young's modulus
  @param[in] nu Poisson ratio
  @param[in] stressState
  @param[in] alpha Thermal expansion coefficients
*/

class3DprintingThermoelasticityStVenantKirchhoff::class3DprintingThermoelasticityStVenantKirchhoff(int ndim,
                                                                                                   string name,
                                                                                                   string name_cons,
                                                                                                   double rho,
                                                                                                   double rhoLiquid,
                                                                                                   double rho_PL,
                                                                                                   double youngRefPowder,
												       double young1_PL,
                                                                                                   double youngRefSolid,
                                                                                                   double young1_SL,
                                                                                                   double young2_SL,
                                                                                                   double nuRefSolid,
                                                                                                   double nu_SL,
                                                                                                   double alphaRefSolid,
												       double logGrowthRate,
                                                                                                   double theta_sigMidpoint,
                                                                                                   double theta_ref,
                                                                                                   double theta_SL,
                                                                                                   double theta_liq,
                                                                                                   int fieldIndex)
        : class3DprintingThermoelasticityStVenantKirchhoffBase(ndim, name, name_cons, rho, fieldIndex) {
    this->rhoNewConf= rhoLiquid;
    this->rho_transNewConf = rho_PL;
    this->youngRefPowder = youngRefPowder;
    this->young1_PL = young1_PL;
    this->youngRefSolid = youngRefSolid;
    this->young1_SL = young1_SL;
    this->young2_SL = young2_SL;
    this->nuRefSolid = nuRefSolid;
    this->nu_SL = nu_SL;
    this->alphaRefSolid = alphaRefSolid;
    this->logGrowthRate = logGrowthRate;
    this->theta_sigMidpoint = theta_sigMidpoint;
    this->theta_ref = theta_ref;
    this->theta_SL = theta_SL;
    this->theta_liq = theta_liq;
    
};


class3DprintingThermoelasticityStVenantKirchhoff::~class3DprintingThermoelasticityStVenantKirchhoff() 
{

};

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double class3DprintingThermoelasticityStVenantKirchhoff::soundSpeed() const
{
    double sound = sqrt(youngRefPowder / _rho);
    return sound;
}

/*! \brief This function predict the Internal Parameters
  @param[out] GaussPoint Gauss Point
  @param[in] dt Time step
  @param[in] timeRun current 
*/
void class3DprintingThermoelasticityStVenantKirchhoff::predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    double T = currState.getExtraFields()[_tempratureFieldIndex];
    
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    
    if ((fabs(intVarsPrev[_intVarsSize + 2]) < 0.5) && (T > theta_liq)) 
    {
        intVarsCurr[_intVarsSize + 2] = 1; //IPhase=0 powder -> IPhase=1 liquid -> IPhase=2 solid
    } 
    else if (fabs(intVarsPrev[_intVarsSize + 2]) > 0.5) 
    {
        if (T > theta_liq) {
            intVarsCurr[_intVarsSize + 2] = 1;
        } else {
            intVarsCurr[_intVarsSize + 2] = 2;
        }
    };
    
    
    // update initial state
    classStates& initState = GaussPoint->getInitialState();
    vector<double>& FIni = initState.getDeformationGradient();
    const vector<double>& FCurr = currState.getDeformationGradient();
    const vector<double>& FPrev = prevState.getDeformationGradient();

    int flagFIni = GaussPoint->getflagFIni();
    if (fabs(intVarsPrev[_intVarsSize + 2]) > 0.5) 
    {
      if (flagFIni == 0 && (fabs(intVarsPrev[_intVarsSize + 2]) > 0.5)) 
        {
            FIni = FPrev;
            GaussPoint->setflagFIni();
        }
    } 
    else 
    {
        static classTensor2 I(getDim(), 1.);
        FIni = I.getData();
    };
};


void class3DprintingThermoelasticityStVenantKirchhoff::computeThermalStretchRatio(classGPs *GaussPoint, double phaseIndex, double theta, double& thermalStretchRatio, 
                                                                                   double& DthermalStretchRatioDT,
                                                                                   double& DDthermalStretchRatioDTDT) const
{

   double thermalStrain;
   double _alpha = alphaRefSolid;
   double theta_Ini;
   
   double initialPhase = GaussPoint->getInitialState().getInternalVariables()[2];
   
   if (fabs(initialPhase) > 1.5)
   {
     
     theta_Ini=theta_ref;
	   
   }else{
	   
     theta_Ini=theta_liq;
     
   }
   
   if (fabs(phaseIndex) < 0.5){
	   
    thermalStretchRatio = 1;
    thermalStrain = 0;
    DthermalStretchRatioDT = 0;
    DDthermalStretchRatioDTDT = 0;
     
   }else{
	   
    thermalStretchRatio = exp(_alpha * (theta - theta_Ini) - _alpha / logGrowthRate * log(1 + exp(logGrowthRate * (theta - theta_sigMidpoint))) + _alpha / logGrowthRate * log(1 + exp(logGrowthRate * (theta_Ini - theta_sigMidpoint))));
    thermalStrain = 1 - thermalStretchRatio;
    DthermalStretchRatioDT = _alpha * exp(_alpha * theta) * pow(1 + exp(logGrowthRate * (theta - theta_sigMidpoint)),-(logGrowthRate + _alpha) / _alpha);
    DDthermalStretchRatioDTDT = _alpha * exp(_alpha * theta) * pow(1 + exp(logGrowthRate * (theta - theta_sigMidpoint)),-(_alpha / logGrowthRate) - 2) * (_alpha - logGrowthRate * exp(logGrowthRate * (theta - theta_sigMidpoint)));
   }
   
};

void class3DprintingThermoelasticityStVenantKirchhoff::computeElasticContants(double phaseIndex, double theta, double& Young, double& nu, double& lambda, double& DlambdaDtheta, double& DDlambdaDthetaDtheta,
                                       double& mu, double& DmuDtheta, double& DDmuDthetaDtheta) const
{
    double DYoungDT, D2YoungDT2;
    double DnuDT, D2nuDT2;
	
    nu = nuRefSolid + nu_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));
    DnuDT = nu_SL * logGrowthRate * (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                 (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
    D2nuDT2 = 2 * nu_SL * logGrowthRate * logGrowthRate *
                  (exp(-2 * logGrowthRate * (theta - theta_sigMidpoint))) /
                  ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                   (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                   (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))))
                  - nu_SL * logGrowthRate * logGrowthRate *
                    (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                    ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                     (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
	
    if (fabs(phaseIndex) < 0.5) 
    {
	Young = youngRefPowder+
                    young1_PL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));
        DYoungDT = young1_PL * logGrowthRate * (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                    ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                     (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
        D2YoungDT2 = 2 * young1_PL * logGrowthRate * logGrowthRate *
                     (exp(-2 * logGrowthRate * (theta - theta_sigMidpoint))) /
                    ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                     (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                     (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))))
                     - young1_PL * logGrowthRate * logGrowthRate *
                     (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                    ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                     (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
		
    } 
    else 
    {
      
        if (theta > theta_SL) 
        {
            Young = youngRefSolid +
                    young1_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));
            DYoungDT = young1_SL * logGrowthRate * (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                       ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                        (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
            D2YoungDT2 = 2 * young1_SL * logGrowthRate * logGrowthRate *
                         (exp(-2 * logGrowthRate * (theta - theta_sigMidpoint))) /
                         ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                          (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                          (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))))
                         - young1_SL * logGrowthRate * logGrowthRate *
                           (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                           ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                            (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
        } else 
        {
            Young = young2_SL * (theta - theta_SL) + youngRefSolid +
                    young1_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));
            DYoungDT = young2_SL +
                       young1_SL * logGrowthRate * (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                       ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                        (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
            D2YoungDT2 = 2 * young1_SL * logGrowthRate * logGrowthRate *
                         (exp(-2 * logGrowthRate * (theta - theta_sigMidpoint))) /
                         ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                          (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                          (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))))
                         - young1_SL * logGrowthRate * logGrowthRate *
                           (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                           ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                            (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
        }
    }
    
    mu = 0.5 * Young / (1. + nu);
    lambda = (Young * nu) / ((1. + nu) * (1. - 2. * nu));
    
    DlambdaDtheta =
            (Young * nu * (DnuDT + 4. * nu * DnuDT)) / ((2. * nu * nu + nu - 1.) * (2. * nu * nu + nu - 1.)) -
            (nu * DYoungDT) / (nu + 2. * nu * nu - 1.) - (Young * DnuDT) / (nu + 2. * nu * nu - 1.);
    DmuDtheta = DYoungDT / (2. * nu + 2.) - (2. * Young * DnuDT) / ((2. * nu + 2.) * (2. * nu + 2.));

    DDlambdaDthetaDtheta = (Young * nu * (4. * pow(DnuDT, 2) + 4. * nu * D2nuDT2 + D2nuDT2)) /
                         ((2. * pow(nu, 2) + nu - 1.) * (2. * pow(nu, 2) + nu - 1.)) -
                         (nu * D2YoungDT2) / (nu + 2. * pow(nu, 2) - 1.) -
                         (2. * DYoungDT * DnuDT) / (nu + 2. * pow(nu, 2) - 1.) -
                         (Young * D2nuDT2) / (nu + 2. * pow(nu, 2) - 1.) -
                         (2. * Young * nu * (DnuDT + 4. * nu * DnuDT) * (DnuDT + 4. * nu * DnuDT)) /
                         ((2. * pow(nu, 2) + nu - 1.) * (2. * pow(nu, 2) + nu - 1.) * (2. * pow(nu, 2) + nu - 1.)) +
                         (2. * Young * DnuDT * (DnuDT + 4. * nu * DnuDT)) /
                         ((2. * pow(nu, 2) + nu - 1.) * (2. * pow(nu, 2) + nu - 1.)) +
                         (2. * nu * DYoungDT * (DnuDT + 4. * nu * DnuDT)) /
                         ((2. * pow(nu, 2) + nu - 1.) * (2. * pow(nu, 2) + nu - 1.)); 


    DDmuDthetaDtheta = D2YoungDT2 / (2. * nu + 2.) +
                     (8. * Young * DnuDT * DnuDT) / ((2. * nu + 2.) * (2. * nu + 2.) * (2. * nu + 2.)) -
                     (2. * Young * D2nuDT2) / ((2. * nu + 2.) * (2. * nu + 2.)) -
                     (4. * DYoungDT * DnuDT) / ((2. * nu + 2.) * (2. * nu + 2.));
};

void class3DprintingThermoelasticityStVenantKirchhoff::computeDensity(double phaseIndex, double theta, double& density, double& DdensityDtheta) const
{
    if (fabs(phaseIndex) > 0.5) 
    {
        // New reference configuration when material turns from powder to liquid
        density = rhoNewConf;
    } 
    else 
    {
        //Transition from powder to liquid
        double expTerm = exp(-logGrowthRate * (theta - theta_sigMidpoint));
        density = _rho + rho_transNewConf * 1 / (1 + expTerm);
        DdensityDtheta = rho_transNewConf * logGrowthRate*expTerm / (1 + expTerm) / (1 + expTerm);
    }
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void class3DprintingThermoelasticityStVenantKirchhoff::initIntVars(vector<double> &intVars) 
{
    _intVarsSize = intVars.size();
    intVars.push_back(0.0);  //activation time
    intVars.push_back(0.0);  //activation state (1 active, 0 inactive)
    //we store the evolution in the mechanical properties
    intVars.push_back(0.0);  //material phase: 0 powder-1Liquid/Solid
    intVars.push_back(0.0);  //young
    intVars.push_back(0.0);  //nu
    intVars.push_back(0.0);  //thermalStretchRatio
};


void class3DprintingThermoelasticityStVenantKirchhoff::compute(classGPs *GaussPoint, double phaseIndex, const classTensor2& Fini, const classTensor2& F, double theta,  
                classTensor2& P, classTensor2& DPDT, vector<double>& intVarsCurr, const vector<double>& intVarsPrev, 
                bool stiff, classTensor4* dPdF,  classTensor4* ddPdTdF, classTensor2* ddPdTdT) const
{
    int ndim = getDim();
    // compute elastic constant
    double lambda, DlambdaDT, DDlambdaDTDT;
    double mu, DmuDT, DDmuDTDT;
    double Young, nu;
    computeElasticContants(phaseIndex, theta, Young, nu, lambda, DlambdaDT, DDlambdaDTDT, mu, DmuDT, DDmuDTDT);
    
    intVarsCurr[_intVarsSize+3] = Young;
    intVarsCurr[_intVarsSize+4] = nu;
    
    // compute thermal stretch
    double thermalStretchRatio, DthermalStretchRatioDT, DDthermalStretchRatioDTDT;
    computeThermalStretchRatio(GaussPoint, phaseIndex, theta, thermalStretchRatio, DthermalStretchRatioDT, DDthermalStretchRatioDTDT);
    intVarsCurr[_intVarsSize+5] = thermalStretchRatio;
    //
    static classTensor2 C(ndim);
    TensorOperations::rightCauchy(F,C); // C= FT*F
    
    classTensor2 invFini(ndim), invFiniT(ndim);
    TensorOperations::inverseTensor2(Fini, invFini);
    TensorOperations::transpose(invFini, invFiniT);
    classTensor2 invFinitTCinvFini(ndim);
    TensorOperations::dotThreeTensors(invFiniT, C, invFini, invFinitTCinvFini);
    // compute Ee = 0.5(invFinitTCinvFini/v**2- I)
    static classTensor2 Ee(ndim);
    Ee.setAll(0);
    double fact = 0.5/thermalStretchRatio/thermalStretchRatio;
    Ee.axpy(invFinitTCinvFini,fact);
    Ee.addDiagonal(-0.5);
        
    // compute Se as a function of Ee 
    static classTensor2 Se(ndim);  
    double traceEe = Ee.trace();
    static classTensor2 I(ndim, 1.); // unity
    Se.setAll(0);
    Se.axpy(I, traceEe*lambda);
    Se.axpy(Ee, 2.*mu);
    
    double Jini = TensorOperations::determinant(Fini);
    // compute second PK
    classTensor2 S(ndim);
    TensorOperations::dotThreeTensors(invFini, Se, invFiniT, S, thermalStretchRatio * Jini);
    TensorOperations::dot(F, S, P);
        
    // compute DEeDT DSeDT, 
    static classTensor2 DEeDT(ndim);
    double DfactDT = -DthermalStretchRatioDT/thermalStretchRatio/thermalStretchRatio/thermalStretchRatio;
    DEeDT.setAll(0);
    DEeDT.axpy(invFinitTCinvFini,DfactDT);
    
    double DtraceEeDT = DEeDT.trace();
    static classTensor2 DSeDT(ndim);
    DSeDT.setAll(0);
    DSeDT.axpy(I, DtraceEeDT * lambda);
    DSeDT.axpy(I, traceEe * DlambdaDT);
    DSeDT.axpy(Ee, 2.* DmuDT);
    DSeDT.axpy(DEeDT, 2.* mu);
    // compute DSDT
    classTensor2 DSDT(ndim);
    TensorOperations::dotThreeTensors(invFini, DSeDT, invFiniT, DSDT, thermalStretchRatio*Jini);
    TensorOperations::dotThreeTensorsAdd(invFini, Se, invFiniT, DSDT, DthermalStretchRatioDT*Jini);
    TensorOperations::dot(F, DSDT, DPDT);
    
    if (stiff)
    {
        // compute by perturbation
        double tol = 1e-6;
        classTensor2 Pplus(ndim), DPDTplus(ndim);
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< ndim; j++)
            {
                classTensor2 Fplus(F);
                Fplus(i,j) += tol;
                compute(GaussPoint, phaseIndex, Fini, Fplus, theta, Pplus, DPDTplus, intVarsCurr, intVarsPrev, false);
                
                for (int k=0; k < ndim; k++)
                {
                    for (int l=0; l< ndim; l++)
                    {
                        (*dPdF)(k,l,i,j) = (Pplus(k,l)- P(k,l))/tol;
                        (*ddPdTdF)(k,l,i,j) = (DPDTplus(k,l)- DPDT(k,l))/tol;
                    }
                }
                
            }
        }
        
        double tolTemp = tol*(theta+1.);
        compute(GaussPoint, phaseIndex, Fini, F, theta+tolTemp, Pplus, DPDTplus, intVarsCurr, intVarsPrev, false);
        for (int k=0; k < ndim; k++)
        {
            for (int l=0; l< ndim; l++)
            {
                (*ddPdTdT)(k,l) = (Pplus(k,l)- P(k,l))/tolTemp;
            }
        }            
    };
};
