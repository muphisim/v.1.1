//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _stochasticHyperElasticNeoHookean_H_
#define _stochasticHyperElasticNeoHookean_H_

#include "constitutiveModels.h"
#include "maths.h"

/*! \brief Hyper Elastic Neo Hookean constitutive model.
*/
class classStochasticHyperElasticNeoHookean : public constitutiveModels {
    protected:
        classTensor6 *elasticStiffTensor;
        vector<double> Young;
        vector<double> C;
        vector<double> nu;
        vector<double> mu;
        vector<double> lambda;
        vector<double> K1;
        int nstoch;
        vector<double> C_Stoch;
        bool _stochasticMapping;
        
    public:
        classStochasticHyperElasticNeoHookean(int ndim, string name, string name_cons, double rho, double young,
                                                     double nu, vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution);
        virtual ~classStochasticHyperElasticNeoHookean();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return nstoch*getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        
        const vector<double>& getC_Stoch() const { return C_Stoch; }
        
    private:
        void stress(vector<double> &PK2, vector<double> &PK1, const vector<double> &FCurr);
       
};

/*! \brief HyperElastic3DCompressibleNeoHookean stochastic
*/
class classStochasticHyperElastic3DCompressibleNeoHookean : public constitutiveModels {
    protected:
        classTensor6 *elasticStiffTensor;
        vector<double> Young;
        vector<double> C;
        vector<double> nu;
        vector<double> mu;
        vector<double> lambda;
        vector<double> K1;
        int nstoch;
        vector<double> C_Stoch;
        bool _stochasticMapping;
    public:
        classStochasticHyperElastic3DCompressibleNeoHookean(int ndim, string name, string name_cons, double rho, double young,
                                                     double nu, vector<double> Stochastic_numbers,
                                                     vector<double> Stochastic_parameters,
                                                     vector<double> Stochastic_function, string approximation, int order, int resolution, bool stochasticMapping, string distribution);
        virtual ~classStochasticHyperElastic3DCompressibleNeoHookean();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return nstoch*getDim();}
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const {}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        virtual double defoEnergy(classGPs *GaussPoint) const;
};

#endif
