//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\temperature.cpp
  \brief This file contains all functions related to the thermal analysis, uncoupled  and fully-coupled Thermo-Mechanical Analysis.
*/

#include "temperature.h"
#include "maths.h"
#include "classGPs.h"
#include "TensorOperations.h"
#include "3DprintingThermoelasticityStVenantKirchhoff.h"


/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] name_cons Name of the constitutive model
  @param[in] rho density
  @param[in] mechCouling True if a full coupling is considered
*/
ThermalLawBase::ThermalLawBase(int ndim, string name, string name_cons, double rho, bool mechCoupling) : 
                constitutiveModels(ndim, name, name_cons, rho), _mechCoupling(mechCoupling)
{

}

ThermalLawBase::~ThermalLawBase()
{
   
}

void ThermalLawBase::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    if (_mechCoupling)
    {
        _term = new ExtraDofMechanicalFullCouplingTerm(_extraFieldIndexes[0]);
    }
    else
    {
        _term = new ExtraDofTerm(_extraFieldIndexes[0]);
    }
}

double ThermalLawBase::soundSpeed() const {
    return 0;
}

void ThermalLawBase::initIntVars(vector<double> &intVars) 
{
    _intVarsSize = intVars.size(); //Index of current internal variable
}


/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void ThermalLawBase::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    double T = currState.getExtraFields()[_extraFieldIndexes[0]];
    double Tprev = prevState.getExtraFields()[_extraFieldIndexes[0]];
    const classTensor1& gradT = currState.getExtraFieldGradients()[_extraFieldIndexes[0]];
    
    // temperature -dependent of conductivity
    double K, DKDT;
    computeK(GaussPoint,K,DKDT);
    
    int ndim = getDim();
    static classTensor2 Ktensor(ndim);
    Ktensor.setAll(0);
    static classTensor4 dKtensordF(ndim);
    static classTensor2 dKtensordT(ndim);
    dKtensordT.setAll(0);
    // K = det(F)*K*RCinv
    if (_mechCoupling)
    {
        const vector<double>&F = currState.getDeformationGradient();
        classTensor2 Ftensor(ndim, F); // construct tensor from vector
        if (flagTanMod)
        {
            static classTensor2 RC(ndim), RCinv(ndim), DJDF(ndim);
            static classTensor4 DRCDF(ndim);
            TensorOperations::rightCauchy(Ftensor,RC,&DRCDF); // C = FT.F        
            double J = TensorOperations::determinant(Ftensor,&DJDF);
            static classTensor4 DRCinvDRC(ndim);
            TensorOperations::inverseTensor2(RC,RCinv,&DRCinvDRC,true);
           
            Ktensor.axpy(RCinv,K*J);
            //
            dKtensordT.axpy(RCinv,DKDT*J);
            //
            TensorOperations::doubleDot(DRCinvDRC,DRCDF,dKtensordF,K*J);
            TensorOperations::dyadicProductAdd(RCinv,DJDF,dKtensordF,K);
        }
        else
        {
            static classTensor2 RC(ndim), RCinv(ndim);
            TensorOperations::rightCauchy(Ftensor,RC); // C = FT.F        
            double J = TensorOperations::determinant(Ftensor);
            TensorOperations::inverseTensor2(RC,RCinv);
            Ktensor.axpy(RCinv,K*J);
        }
    }
    else
    {
       Ktensor.addDiagonal(K);
       dKtensordT.addDiagonal(DKDT);
    }
    
    classTensor1& flux = currState.getExtraFieldFluxes()[_extraFieldIndexes[0]];
    double& source = currState.getExtraFieldSources()[_extraFieldIndexes[0]];
    double& fieldDensity = currState.getExtraFieldDensities()[_extraFieldIndexes[0]];
    
    double Rho, DRhoDT;
    computeRho(GaussPoint, Rho,DRhoDT);
    double C, DCDT;
    computeC(GaussPoint,C,DCDT);
    
    fieldDensity = Rho*C;
    TensorOperations::dot(Ktensor,gradT, flux, -1.); //q= -Ktensor * gradT
    
    // source must be added up to existing one
    source += Rho*C*(T-Tprev)/dt;
    
    if (flagTanMod)
    {
        double& DsourceDT = (*currState.getDextraFieldSourcesDextraFields())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        classTensor1& DfluxDT = (*currState.getDextraFieldFluxesDextraFields())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        classTensor2& DfluxDgradT = (*currState.getDextraFieldFluxesDextraFieldGrads())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        
        // add-up
        DsourceDT += Rho*C/dt + (DRhoDT*C+Rho*DCDT)*(T-Tprev)/dt;
        TensorOperations::dot(dKtensordT,gradT, DfluxDT, -1.); 
        DfluxDgradT = Ktensor;
        DfluxDgradT.scale(-1.);
        
        if (_mechCoupling)
        {
            classTensor3& DfluxDF = (*currState.getDextraFieldFluxesDdeformationGradient())[_extraFieldIndexes[0]];
            for (int i=0; i< ndim; i++)
            {
                for (int j=0; j< ndim; j++)
                {
                    for (int p=0; p< ndim; p++)
                    {
                        for (int q=0; q < ndim; q++)
                        {
                            DfluxDF(i,p,q) -= dKtensordF(i,j,p,q)*gradT(j);
                        }
                    }
                }
            }
            
        }
    }
};


/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] C_0 Specific heat at reference temperature
  @param[in] C_1 Temperature dependence of specific heat
  @param[in] k_0 Conductivity at reference temperature (current configuration)
  @param[in] k_1 Temperature dependence of specific heat of conductivity (current configuration)
  @param[in] r Volumetric heat supplied externally into the body
  @param[in] theta_in Initial temperature at which the material is defined as undeformed in the absence of mechanical loads
  @param[in] theta_ref1 reference temperature of C_0
  @param[in] theta_ref2 reference temperature of k_0
  @param[in] mechCouling True if a full coupling is considered
*/
classtemperature::classtemperature(int ndim, string name, string name_cons, double rho, double C_0,
                                                   double C_1, double k_0, double k_1, double theta_ini,
                                                   double theta_ref1, double theta_ref2, 
                                                   bool mechCoupling) : 
                                    ThermalLawBase(ndim, name, name_cons, rho, mechCoupling) {
    _C_0 = C_0;
    _k_0 = k_0;
    _C_1 = C_1;
    _k_1 = k_1;
    _theta_ini = theta_ini;
    _theta_ref1 = theta_ref1;
    _theta_ref2 = theta_ref2;
}

classtemperature::~classtemperature()
{
   
}

void classtemperature::computeRho(classGPs *GaussPoint, double& Rho, double& dRhodtheta) const
{
    Rho = getRho();
    dRhodtheta = 0.;
}
void classtemperature::computeK(classGPs *GaussPoint, double& K, double& dKdtheta) const
{
    double T = GaussPoint->getCurrentState().getExtraFields()[_extraFieldIndexes[0]];
    K = _k_0 + _k_1 * (T - _theta_ref2);
    dKdtheta = _k_1;
}
void classtemperature::computeC(classGPs *GaussPoint, double& C, double& dCdtheta) const
{
    double T = GaussPoint->getCurrentState().getExtraFields()[_extraFieldIndexes[0]];
    C = _C_0 + _C_1 * (T - _theta_ref1);
    dCdtheta = _C_1;
};


/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] C_0 Specific heat at reference temperature
  @param[in] C_1 Temperature dependence of specific heat
  @param[in] k_0 Conductivity at reference temperature (current configuration)
  @param[in] k_1 Temperature dependence of specific heat of conductivity (current configuration)
  @param[in] r Volumetric heat supplied externally into the body
  @param[in] theta_in Initial temperature at which the material is defined as undeformed in the absence of mechanical loads
  @param[in] theta_ref1 reference temperature of C_0
  @param[in] theta_ref2 reference temperature of k_0
*/
classtemperature3Dprinting::classtemperature3Dprinting(int ndim, string name, string name_cons, double rho,
                                                       double CRefSolid,  double C1_SL,
                                                       double C2_SL, double latentHeat1_PL, double latentHeat2_PL,
                                                       double latentHeat3_PL, double latentHeat1_SL,
                                                       double latentHeat2_SL, double latentHeat3_SL,
                                                       double kRefPowder, double k1_PL,
                                                       double kRefSolid, double k1_SL, double k2_SL,
                                                       double k3_SL, double logGrowthRate,
                                                       double theta_sigMidpoint, double theta_ref,
                                                       double theta_powder, double theta_SL, double theta_liq,
                                                       bool mechCoupling)
        : ThermalLawBase(ndim, name, name_cons, rho, mechCoupling) {

    this->CRefSolid = CRefSolid;
    this->C1_SL = C1_SL;
    this->C2_SL = C2_SL;
    this->latentHeat1_PL = latentHeat1_PL;
    this->latentHeat2_PL = latentHeat2_PL;
    this->latentHeat3_PL = latentHeat3_PL;
    this->latentHeat1_SL = latentHeat1_SL;
    this->latentHeat2_SL = latentHeat2_SL;
    this->latentHeat3_SL = latentHeat3_SL;
    this->kRefPowder = kRefPowder;
    this->k1_PL = k1_PL;
    this->kRefSolid = kRefSolid;
    this->k1_SL = k1_SL;
    this->k2_SL = k2_SL;
    this->k3_SL = k3_SL;
    this->logGrowthRate = logGrowthRate;
    this->theta_sigMidpoint = theta_sigMidpoint;
    this->theta_ref = theta_ref;
    this->theta_powder = theta_powder;
    this->theta_SL= theta_SL;
    this->theta_liq = theta_liq;
    
};

void classtemperature3Dprinting::computeRho(classGPs *GaussPoint, double& Rho, double& dRhodtheta) const
{
    const constitutiveModels* firstLaw = GaussPoint->getConstitutiveManager().getConsModel()[0];
    const class3DprintingThermoelasticityStVenantKirchhoffBase* thermoLaw = dynamic_cast<const class3DprintingThermoelasticityStVenantKirchhoffBase*>(firstLaw);
    if (firstLaw == NULL)
    {
        ERROR("class3DprintingThermoelasticityStVenantKirchhoffBase must be used with classtemperature3Dprinting");
        exit(-1);
    }
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    double phaseIndex = intVarsCurr[2];    
    double theta = currState.getExtraFields()[_extraFieldIndexes[0]];
    thermoLaw->computeDensity(phaseIndex, theta, Rho, dRhodtheta);
    intVarsCurr[_intVarsSize+2] = Rho;
}

void classtemperature3Dprinting::computeK(classGPs *GaussPoint, double& k, double& dkdtheta) const
{
    const constitutiveModels* firstLaw = GaussPoint->getConstitutiveManager().getConsModel()[0];
    const class3DprintingThermoelasticityStVenantKirchhoffBase* thermoLaw = dynamic_cast<const class3DprintingThermoelasticityStVenantKirchhoffBase*>(firstLaw);
    if (firstLaw == NULL)
    {
        ERROR("class3DprintingThermoelasticityStVenantKirchhoffBase must be used with classtemperature3Dprinting");
        exit(-1);
    }
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    double phaseIndex = intVarsCurr[2];    
    double theta = currState.getExtraFields()[_extraFieldIndexes[0]];
    
    if (fabs(phaseIndex) < 0.5) {

        k = kRefPowder + k1_PL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));
        dkdtheta = k1_PL  * logGrowthRate * (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                   ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                    (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
     
      
    } 
 else if (fabs(phaseIndex ) > 0.5) {
        
    if (theta < theta_ref){
        k = kRefSolid + k1_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));
        dkdtheta = k1_SL* logGrowthRate * (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                   ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                    (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
    } else if (theta <= theta_SL){ 
	    k = kRefSolid + k1_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) + k2_SL * (theta - theta_ref);
        dkdtheta = k1_SL* logGrowthRate * (exp(-logGrowthRate * (theta - theta_sigMidpoint))) /
                   ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) *
                    (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))))+ k2_SL;
    } else { 
	    k = k3_SL + k1_SL *1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));
        dkdtheta = k1_SL * logGrowthRate * exp(theta - theta_sigMidpoint) / ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)))*(1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
    }
     
      }
    intVarsCurr[_intVarsSize+0] = k;
};
void classtemperature3Dprinting::computeC(classGPs *GaussPoint, double& C, double& dCdtheta) const
{
    const constitutiveModels* firstLaw = GaussPoint->getConstitutiveManager().getConsModel()[0];
    const class3DprintingThermoelasticityStVenantKirchhoffBase* thermoLaw = dynamic_cast<const class3DprintingThermoelasticityStVenantKirchhoffBase*>(firstLaw);
    if (firstLaw == NULL)
    {
        ERROR("class3DprintingThermoelasticityStVenantKirchhoffBase must be used with classtemperature3Dprinting");
        exit(-1);
    }
    
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    double phaseIndex = intVarsCurr[2];   // TODO change this to couple this law with any mechancial law
    double theta = currState.getExtraFields()[_extraFieldIndexes[0]];
    
    if (fabs(phaseIndex) < 0.5) {
        if (theta >= theta_ref &&  theta <= theta_powder) {
            C = CRefSolid + C1_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) + C2_SL * (theta - theta_ref);
            dCdtheta = C1_SL * logGrowthRate * exp(theta - theta_sigMidpoint) / ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) * (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)))) + C2_SL;

        } else if (theta > theta_powder && theta < theta_liq) {
            C = latentHeat1_PL * theta * theta + latentHeat2_PL * theta + latentHeat3_PL;
            dCdtheta = latentHeat1_PL * 2 * theta + latentHeat2_PL;

        } else {
            C = CRefSolid + C1_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));  
            dCdtheta = C1_SL * logGrowthRate * exp(theta - theta_sigMidpoint) / ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) * (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
        }

    } else if (fabs(phaseIndex) > 0.5) {
        if (theta >= theta_ref &&  theta <= theta_SL) {
            C = CRefSolid + C1_SL * 1 / (1 + exp(-logGrowthRate * (theta-theta_sigMidpoint))) + C2_SL * (theta - theta_ref);
            dCdtheta = C1_SL * logGrowthRate * exp(theta - theta_sigMidpoint) / ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) * (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)))) + C2_SL;

        } else if (theta > theta_SL && theta < theta_liq) {
            C = latentHeat1_SL * theta * theta + latentHeat2_SL * theta + latentHeat3_SL;
            dCdtheta = latentHeat1_SL * 2 * theta + latentHeat2_SL;

        } else {
            C = CRefSolid + C1_SL * 1 / (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint)));  
            dCdtheta = C1_SL * logGrowthRate * exp(theta - theta_sigMidpoint) / ((1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))) * (1 + exp(-logGrowthRate * (theta - theta_sigMidpoint))));
        }
    }
    intVarsCurr[_intVarsSize+1] = C;
};


void classtemperature3Dprinting::initIntVars(vector<double> &intVars) 
{
    _intVarsSize = intVars.size(); //Index of current internal variable
    intVars.push_back(0.0);  //Conductivity  value
    intVars.push_back(0.0);  //Specific heat value
    intVars.push_back(0.0);  //Density
}
