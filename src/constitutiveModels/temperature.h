//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _temperature_H_
#define _temperature_H_

#include "constitutiveModels.h"

/*! \brief  This file contains all functions related to the thermal analysis, uncoupled  and fully-coupled Thermo-Mechanical Analysis.
*/

class ThermalLawBase : public constitutiveModels 
{
    protected:
        int _intVarsSize; //Size of intVars before adding the internal variables of the extra model
        int _mechCoupling;
        
    public:
        ThermalLawBase(int ndim, string name, string name_cons, double rho, bool mechCoupling);
        virtual ~ThermalLawBase();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return 1;}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        
    protected:
        virtual void computeRho(classGPs *GaussPoint, double& Rho, double& dRhodtheta) const = 0;
        virtual void computeK(classGPs *GaussPoint, double& K, double& dKdtheta) const = 0;
        virtual void computeC(classGPs *GaussPoint, double& C, double& dCdtheta) const = 0;
 
};

// properties with linear dependency
class classtemperature : public ThermalLawBase 
{
    protected:
        double _k_0; 
        double _k_1;
        double _C_0;
        double _C_1; 
        double _theta_ini;
        double _theta_ref1;
        double _theta_ref2;

    public:
        classtemperature(int ndim, string name, string name_cons, double rho, double C_0, double C_1,
                                 double k_0, double k_1, double theta_ini, double theta_ref1, double theta_ref2,
                                 bool mechCouling);
        virtual ~classtemperature();
     
     protected:
        virtual void computeRho(classGPs *GaussPoint, double& Rho, double& dRhodtheta) const;
        virtual void computeK(classGPs *GaussPoint, double& K, double& dKdtheta) const;
        virtual void computeC(classGPs *GaussPoint, double& C, double& dCdtheta) const;
};

/*! \brief  This file contains all functions related to the thermal analysis, uncoupled  and fully-coupled Thermo-Mechanical Analysis.
*/
class classtemperature3Dprinting : public ThermalLawBase 
{
    protected:
        double CRefSolid;
        double C1_SL;
        double C2_SL;
        double latentHeat1_PL;
        double latentHeat2_PL;
        double latentHeat3_PL;
        double latentHeat1_SL;
        double latentHeat2_SL;
        double latentHeat3_SL;
        double kRefPowder;
        double k1_PL;
        double kRefSolid;
        double k1_SL;
        double k2_SL;
        double k3_SL;
        double logGrowthRate;
        double theta_sigMidpoint;
        double theta_ref;
        double theta_powder;
        double theta_SL;
        double theta_liq;

    public:
        classtemperature3Dprinting(int ndim, string name, string name_cons, double rho, double CRefSolid,
                                   double C1_SL, double C2_SL, double latentHeat1_PL, 
                                   double latentHeat2_PL, double latentHeat3_PL,
                                   double latentHeat1_SL, double latentHeat2_SL, double latentHeat3_SL,
                                   double kRefPowder, double k1_PL, double kRefSolid,
                                   double k1_SL, double k2_SL, double k3_SL, double logGrowthRate,
                                   double theta_sigMidpoint, double theta_ref, double theta_powder,
                                   double theta_SL, double theta_liq, bool mechCoupling);
        virtual ~classtemperature3Dprinting(){}
        
        virtual void initIntVars(vector<double> &intVars);
        
    protected:
        virtual void computeRho(classGPs *GaussPoint, double& Rho, double& dRhodtheta) const;
        virtual void computeK(classGPs *GaussPoint, double& K, double& dKdtheta) const;
        virtual void computeC(classGPs *GaussPoint, double& C, double& dCdtheta) const;

};
#endif
