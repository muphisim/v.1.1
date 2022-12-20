//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _J2plasticityFiniteStrain_H_
#define _J2plasticityFiniteStrain_H_

#include "constitutiveModels.h"

/*! \brief J2 plasticity model based on the bi-logarithmic strain energy and isotropic hardening 
*/

class IsotropicHardening
{
    public:
        static IsotropicHardening* allocate(const char what[], const std::vector<double>& data);
    public:
        IsotropicHardening(){}
        IsotropicHardening(const IsotropicHardening& src){} 
        virtual ~IsotropicHardening(){}
        
        // get initial yield stress
        virtual double getR0() const = 0;
        // from plastic deformation, get yield stress, hardening modulus
        virtual void hardening(double p, double& R, double& dR) const = 0;
        // to identify the law
        virtual std::string getStrName() const = 0;
        // clone
        virtual IsotropicHardening* clone() const = 0;
        
};

class linearIsotropicHardening: public IsotropicHardening
{
    protected:
        double _sy0, _H;
    
    public:
        linearIsotropicHardening(double sy0, double H): IsotropicHardening(), _sy0(sy0), _H(H){}
        linearIsotropicHardening(const linearIsotropicHardening& src): IsotropicHardening(src),
                   _sy0(src._sy0), _H(src._H) {} 
        virtual ~linearIsotropicHardening(){};
        
         // get initial yield stress
        virtual double getR0() const {return _sy0;};
        // from plastic deformation, get yield stress, hardening modulus
        virtual void hardening(double p, double& R, double& dR) const 
        {
            R = _sy0 + _H*p;
            dR = _H;
        };
        // to identify the law
        virtual std::string getStrName() const {return "LinISO";};
        
        virtual IsotropicHardening* clone() const {return new linearIsotropicHardening(*this);}
};

class SwiftIsotropicHardening: public IsotropicHardening
{
    protected:
        double _sy0, _eps0, _n;
    
    public:
        SwiftIsotropicHardening(double sy0, double p0, double n): IsotropicHardening(), _sy0(sy0), _eps0(p0), _n(n){}
        SwiftIsotropicHardening(const SwiftIsotropicHardening& src): IsotropicHardening(src),
                   _sy0(src._sy0), _eps0(src._eps0), _n(src._n) {} 
        virtual ~SwiftIsotropicHardening(){};
        
         // get initial yield stress
        virtual double getR0() const {return _sy0;};
        // from plastic deformation, get yield stress, hardening modulus
        virtual void hardening(double p, double& R, double& dR) const 
        {
            R = _sy0*pow(1+p/_eps0, _n);
            dR = _sy0*_n*pow(1+p/_eps0, _n-1)/_eps0;
        };
        // to identify the law
        virtual std::string getStrName() const {return "SwiftISO";};
        
        virtual IsotropicHardening* clone() const {return new SwiftIsotropicHardening(*this);}
};

class VoceIsotropicHardening: public IsotropicHardening
{
    protected:
        double _sy0, _K, _C;
    
    public:
        VoceIsotropicHardening(double sy0, double K, double C): IsotropicHardening(), _sy0(sy0), _K(K), _C(C){}
        VoceIsotropicHardening(const VoceIsotropicHardening& src): IsotropicHardening(src),
                   _sy0(src._sy0), _K(src._K), _C(src._C) {} 
        virtual ~VoceIsotropicHardening(){};
        
         // get initial yield stress
        virtual double getR0() const {return _sy0;};
        // from plastic deformation, get yield stress, hardening modulus
        virtual void hardening(double p, double& R, double& dR) const 
        {
            R = _sy0 + _K*(1.-exp(-_C*p));
            dR = _K*_C*exp(-_C*p);
        };
        // to identify the law
        virtual std::string getStrName() const {return "VoceISO";};
        
        virtual IsotropicHardening* clone() const {return new VoceIsotropicHardening(*this);}
};



class J2plasticityFiniteStrain : public constitutiveModels
{
    protected:
        double _E, _nu, _K, _mu; // young modulus and Poisson ratio
        IsotropicHardening* _isoHardening;
        int _intVarsSize; // first location in the internal variable
        int _order; // strain order
        classTensor4 _Cel, _I4dev;
        
    public:
        J2plasticityFiniteStrain(int ndim, string name, string name_cons, double rho, 
                                 double E, double nu, const IsotropicHardening& hardenLaw,
                                 int order);
        
        virtual ~J2plasticityFiniteStrain();
        virtual void initLawFromOptions();
        virtual int getNbrDofConsMod() const {return getDim();}; 
        virtual double soundSpeed() const;
        virtual void initIntVars(vector<double> &intVars);
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const{}
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const{}
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod);  
};
#endif //_J2plasticityFiniteStrain_H_