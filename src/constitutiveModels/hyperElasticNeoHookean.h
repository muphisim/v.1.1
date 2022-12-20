//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _hyperElasticNeoHookean_H_
#define _hyperElasticNeoHookean_H_

#include "constitutiveModels.h"

/*! \brief Ted_Belytschko,_Wing_Kam_Liu,_Brian_Moran 2000. Nonlinear finite elements for continua and structures.
*/
class classHyperElasticNeoHookean : public constitutiveModels {
    protected:
        double _Young;
        double _nu;
        double _lambda;
        double _mu;
        double _K1;// Bulk modulus
    public:
        classHyperElasticNeoHookean(int ndim, string name, string name_cons, double rho, double Young, double nu);
        virtual ~classHyperElasticNeoHookean();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const {}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
};


/*! \brief Arroyo, Ortiz. Local maximum-entropy approximation schemes: a seamless
bridge between ﬁnite elements and meshfree methods
*/
class classHyperElastic3DCompressibleNeoHookean : public constitutiveModels {
    protected:
        double _Young;
        double _nu;
        double _lambda;
        double _mu;
        double _K1;// Bulk modulus
    public:
        classHyperElastic3DCompressibleNeoHookean(int ndim, string name, string name_cons, double rho, double Young, double nu);
        virtual ~classHyperElastic3DCompressibleNeoHookean();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const {}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        virtual double defoEnergy(classGPs *GaussPoint) const;
};

#endif
