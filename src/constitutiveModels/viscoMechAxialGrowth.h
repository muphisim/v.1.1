//
//
// File authors:  see Authors.txt
// Description: nonlinear viscoelasticity (Simo J.C., Computational Inelasticity, Springer) + Growth constitutive models
// 
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _viscoMechAxialGrowth_H_
#define _viscoMechAxialGrowth_H_

#include "constitutiveModels.h"

/*! \brief This function is a wrapper to couple viscoelastic and growth models. See the references in inside those classes */
class classViscoMechAxialGrowth : public constitutiveModels {
    protected:
        constitutiveModels *viscoElastic; // The calculation of the growth part will be inside the update of this constitutive model. I am not calling the growth const model

        int _finalSize = 0;
        int _N = 0;
        double Young = 0;
        double _K = 0, _sum_mu = 0;
        double nu = 0;
        double lambda = 0; /*! Lame's first parameter*/
        double mu = 0; /*! Shear modulus*/
        double Gc0 = 0; /*! Growth multiplier*/ // Same notation than Ellen Kuhl
        vector<double> n0; /*! Normal to the plane of growth */
        double Kp = 0, Kd = 0; // Polymerisation and depolymerisation rate constants
        double rhoMicro = 0; // Proportion of the density that is microtubules
    public:
        classViscoMechAxialGrowth(int ndim, string name, string name_cons, double rho, double K, int N, vector<double> mus,
                                  vector<double> etas, double Gc, vector<double> n0, double Kp, double Kd, double rhoMicro);
        virtual ~classViscoMechAxialGrowth();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
};

#endif
