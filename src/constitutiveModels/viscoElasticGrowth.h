//
//
// File authors:  see Authors.txt
// Description: nonlinear viscoelasticity (Simo J.C., Computational Inelasticity, Springer) + Growth constitutive models
// 
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _viscoElasticGrowth_H_
#define _viscoElasticGrowth_H_

#include "constitutiveModels.h"

/*! \brief This function is a wrapper to couple viscoelastic and growth models. See the references in inside those classes */
class classViscoElasticGrowth : public constitutiveModels {
    protected:
        constitutiveModels *viscoElastic, *growth;
        int _finalSize;
        int _N;
        double _Young;
        double _nu;
        
    public:
        classViscoElasticGrowth(int ndim, string name, string name_cons, double rho, double K, int N, vector<double> mus,
                                vector<double> etas, double Gc);
        virtual ~classViscoElasticGrowth();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
};

#endif
