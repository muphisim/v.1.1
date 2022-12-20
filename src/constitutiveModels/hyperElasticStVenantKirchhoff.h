//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _hyperElasticStVenantKirchhoff_H_
#define _hyperElasticStVenantKirchhoff_H_

#include "constitutiveModels.h"

/*! \brief Hyper Elastic St-Venant - Kirchhoff constitutive model.
*/
class classHyperElasticStVenantKirchhoff : public constitutiveModels 
{
    protected:
        classTensor4 _elasticStiffTensor;
        double _Young;
        double _nu;
        int _stressState; //0 for plain stress, 1 for plain strain

    public:
        classHyperElasticStVenantKirchhoff(int ndim, string name, string name_cons, double rho, double young, double nu,
                                           int stressState);
        virtual ~classHyperElasticStVenantKirchhoff();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const;
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        virtual double defoEnergy(classGPs *GaussPoint) const;
    
    protected:
        void fillElasticStiffTensor(int ndim, classTensor4& elasticStiffTensor) const;
        void stress(vector<double> &Sig, const vector<double> &E) const;
};

#endif
