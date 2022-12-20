//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _isoMorphogeneticGrowth_H_
#define _isoMorphogeneticGrowth_H_

#include "constitutiveModels.h"

/*! \brief Ellen Kuhl formulation for the area growth, Journal of Mechanical Behavior of Biomedical Materials, 29 (2014), 529-543.
Decomposition of the deformation gradient into the elastic (Hyperelastic) and growth contribution: F=Fe*Fg
*/
class classIsoMorphogeneticGrowth : public constitutiveModels {
    protected:
        double _Young;
        double _nu;
        double _lambda;
        double _mu;
        double _Gc; // Same notation than Ellen Kuhl

    public:
        classIsoMorphogeneticGrowth(int ndim, string name, string name_cons, double rho, double Young, double nu,
                                    double Gc);
        virtual ~classIsoMorphogeneticGrowth();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}; 
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);  

        // For external calls
        void getFg(const vector<double> &FCurr, const vector<double> &FPrev, const vector<double> &intVarsPrev,
                    vector<double> &intVarsCurr, double dt, double timeRun, bool flagTanMod,
                    const vector<double> &PK1, classTensor4 *tanModuli, vector<double> &Fg) const ;
        void getE(int ndim, vector<double> &E, const vector<double> &Fcurr, const vector<double> &Fprev, const vector<double> &intVars,
                  double dt, double timeRun, int nStep) const;
};

#endif
