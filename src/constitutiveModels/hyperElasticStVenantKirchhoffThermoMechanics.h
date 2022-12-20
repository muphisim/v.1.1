//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _thermoMechanics_H_
#define _thermoMechanics_H_

#include "hyperElasticStVenantKirchhoff.h"
/*! \brief This function is a wrapper to couple Hyper Elastic St-Venant-Kirchhoff model
 * with a thermomechanical couling term
 * */
 
class classHyperElasticStVenantKirchhoffThermoMechanics : public classHyperElasticStVenantKirchhoff {
    protected:
        classTensor2 _alphatensor;
        double _theta_ref; // reference temperature
        int _stressState; //0 for plain stress, 1 for plain strain
        int _intVarsSize; //Size of intVars before adding the internal variables of the extra model
        int _tempratureFieldIndex;  /*! index of temperature field in the list of extraDofs*/

    public:
        classHyperElasticStVenantKirchhoffThermoMechanics(int ndim, string name, string name_cons, double rho, double young,
                                                          double nu, int stressState, double alpha, double theta_ref, 
                                                          int fieldIndex=0);
        virtual  ~classHyperElasticStVenantKirchhoffThermoMechanics();
        virtual void initLawFromOptions();
        virtual void initIntVars(vector<double> &intVars);
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
};

#endif
