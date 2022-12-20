//
//
// File authors:  see Authors.txt
// Description: nonlinear viscoelasticity (Simo J.C., Computational Inelasticity, Springer)
// 
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _viscoElastic_H_
#define _viscoElastic_H_

#include "constitutiveModels.h"

/*! \brief Nonlinear viscoelasticity (Simo J.C., Computational Inelasticity, Springer)
 */
class classViscoElastic : public constitutiveModels {

    protected:
        double _rho;                 // density
        double _E;                   // Young's modulus
        double _nu;                  // Poisson's ratio
        double _K;                   // bulk modulus
        int _N;                      // number of internal variables
        double eps = 1.e-8;
        double tol = 1.e-6;
        double _timeStep;
        int _finalSize;

        vector<double> _eta;    // Viscosity vector, it starts with eta1 for the first classViscoelastic branch (spring+damper)
        vector<double> _mu;     // shear modulus vector, it starts with mu_inf at t = inf
        vector<double> _tau;    // relaxation time
        
    public:

        classViscoElastic(int ndim, string name, string name_cons, double rho, double K, int N, 
                          const vector<double>& mus,
                          const vector<double>& etas);
        virtual ~classViscoElastic();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}
        virtual double soundSpeed() const; 
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{};
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{};
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);
        
        
        void updateConstitutive(const vector<double> &Fn,  const vector<double> &F0, const vector<double> &intVarsPrev, vector<double> &intVarsCurr,
                       double dT, double timeRun, bool stiff, vector<double> &P, const vector<double> &PK1Prev, classTensor4 *tanModuli);
                       
    private:
        void update(const vector<double> &F0,const vector<double> &F, vector<double> &P, vector<double> &IPvis,
                const vector<double> &IPvisprev, classTensor4 *&tanModuli, bool stiff) const;

        void updateIP(const vector<double> &F1, const vector<double> &kirchhoff1_, vector<double> &IPvis,
                                 const vector<double> &IPvisprev, double &g, vector<double> &hn, vector<double> &r,
                                 classTensor4 *&tanPartD, classTensor4 *&ICinv, classTensor4 *&CinvTensorCinv,
                                 const vector<double> &C, const vector<double> &Cinv) const;

        void stress(vector<double> &Sig_, const vector<double> &Fn, vector<double> &IPvis, const vector<double> &IPvisprev,
                    classTensor4 *&tanModulvectori) const;

        void initialKirchhoffStressTensor(const vector<double> &Fn, vector<double> &Kirchhoff_) const;
    
};

#endif
