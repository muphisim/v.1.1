//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _volMechDrivenGrowth_H_
#define _volMechDrivenGrowth_H_

#include "constitutiveModels.h"

/*! \brief Ellen Kuhl formulation for the area growth, Journal of Mechanical Behavior of Biomedical Materials, 29 (2014), 529-543.
Decomposition of the deformation gradient into the elastic (Hyperelastic) and growth contribution: F=Fe*Fg
*/
class classVolMechDrivenGrowth : public constitutiveModels {
    protected:
        double _lambda; /*! Lame's first parameter*/
        double _mu; /*! Shear modulus*/
        double _Gc; /*! Growth multiplier*/ // Same notation than Ellen Kuhl
        vector<double> _n0; /*! Normal to the plane of growth */
        double _Young; /*! Young modulus*/
        double _nu; /*! Poisson ratio*/

    public:
        classVolMechDrivenGrowth(int ndim, string name, string name_cons, double rho, double Young, double nu, double Gc);
        virtual ~classVolMechDrivenGrowth();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
};

#endif
