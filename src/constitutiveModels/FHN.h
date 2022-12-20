//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _FHN_H_
#define _FHN_H_

#include "constitutiveModels.h"
#include "maths.h"


/*! \brief Ellen Kuhl formulation for the area growth, Journal of Mechanical Behavior of Biomedical Materials, 29 (2014), 529-543.
Decomposition of the deformation gradient into the elastic (Hyperelastic) and growth contribution: F=Fe*Fg
*/
class classFHN : public constitutiveModels {
    protected:
        double _vdotequi;
        double _aelec;
        double _b;
        double _epsi;
        double _vequi;
        classTensor1 _n0;
        double _D_iso;
        double _D_ani;
        double _stimulusCurrent;

        int _intVarsSize; //Size of intVars before adding the internal variables of the extra model
        bool _mechCouling; // true if performing mechanical coupling

    public:
        classFHN(int ndim, string name, string name_cons, double rho, double aelec, double b, double epsi, double vequi,
                 const vector<double> &n0, double D_iso, double D_ani, double stimulusCurrent, 
                 bool mechCouling =false);
        virtual ~classFHN();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return 1;}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const {}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
};
#endif