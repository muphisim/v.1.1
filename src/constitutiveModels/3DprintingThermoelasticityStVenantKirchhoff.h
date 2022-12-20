//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _hyperElasticStVenantKirchhoffthermomecha_H_
#define _hyperElasticStVenantKirchhoffthermomecha_H_

#include "constitutiveModels.h"


class class3DprintingThermoelasticityStVenantKirchhoffBase : public constitutiveModels 
{
    protected:
        int _intVarsSize; //Size of intVars before adding the internal variables of the extra model
        int _tempratureFieldIndex;  /*! index of temperature field in the list of extraDofs*/
        
    public:
        class3DprintingThermoelasticityStVenantKirchhoffBase(int ndim, string name, string name_cons, double rho, int fieldIndex); 
    
        virtual ~class3DprintingThermoelasticityStVenantKirchhoffBase(){};
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const = 0;
        virtual void initIntVars(vector<double> &intVars) = 0;        
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const;
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const = 0;
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        
    
    public:
        // constitutive law with temperature dependence will be specif
        virtual void computeDensity(double phaseIndex, double theta, double& density, double& DdensityDtheta) const = 0;
        virtual void compute(classGPs *GaussPoint, double phaseIndex, const classTensor2& Fini, const classTensor2& F, double theta,  classTensor2& P, classTensor2& DPDT,
                    vector<double>& intVarsCurr, const vector<double>& intVarsPrev, 
                    bool stiff=false, classTensor4* dPdF=NULL,  classTensor4* ddPdTdF=NULL, classTensor2* ddPdTdT=NULL) const = 0;
};

/*! \brief Hyper Elastic St-Venant - Kirchhoff constitutive model.
*/
class class3DprintingThermoelasticityStVenantKirchhoff : public class3DprintingThermoelasticityStVenantKirchhoffBase {
    protected:
        double rhoNewConf;
        double rho_transNewConf;
        double youngRefPowder;
        double young1_PL;
        double youngRefSolid;
        double young1_SL;
        double young2_SL;
        double nuRefSolid;
        double nu_SL;
        double alphaRefSolid;
        double logGrowthRate;
        double theta_sigMidpoint;
        double theta_ref;
        double theta_SL;
        double theta_liq;
    
    public:
        class3DprintingThermoelasticityStVenantKirchhoff(int ndim, string name, string name_cons, double rho,
                                                         double rhoLiquid, double rho_PL,
                                                         double youngRefPowder, double young1_PL,
                                                         double youngRefSolid, double young1_SL,
                                                         double young2_SL,double nuRefSolid, double nu_SL, 
                                                         double alphaRefSolid, double logGrowthRate,
                                                         double theta_sigMidpoint, double theta_ref,
                                                         double theta_SL, double theta_liq,
                                                         int fieldIndex=0);
        
        virtual ~class3DprintingThermoelasticityStVenantKirchhoff();
        virtual double soundSpeed() const;
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const;
        virtual void initIntVars(vector<double> &intVars);  
        virtual void computeDensity(double phaseIndex, double theta, double& density, double& DdensityDtheta) const;
   
        virtual void compute(classGPs *GaussPoint, double phaseIndex, const classTensor2& Fini, const classTensor2& F, double theta,  classTensor2& P, classTensor2& DPDT,
                    vector<double>& intVarsCurr, const vector<double>& intVarsPrev, 
                    bool stiff=false, classTensor4* dPdF=NULL,  classTensor4* ddPdTdF=NULL, classTensor2* ddPdTdT=NULL) const;
     
     private:           
        void computeElasticContants(double phaseIndex, double T, double& Young, double& nu, double& lambda, double& DlambdaDT, double& DDlambdaDTDT,
                                           double& mu, double& DmuDT, double& DDmuDTDT) const;
        void computeThermalStretchRatio(classGPs *GaussPoint, double phaseIndex, double theta, double& thermalStretchRatio, 
                                           double& DthermalStretchRatioDT, double& DDthermalStretchRatioDTDT) const;
        
};

#endif

