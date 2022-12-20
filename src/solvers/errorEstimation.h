//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file errorEstimation.h
  \brief This file contains all functions related the error analysis
*/

#ifndef _errorEstimation_H_
#define _errorEstimation_H_


#include <vector>
#include "configuration.h"

class trueDisplacementField
{
    public:
        virtual ~trueDisplacementField(){}
        virtual std::string getName() const = 0;
        
        virtual double getL2Error(const vector<double> &uK, 
                         const vector<classGPs *> &GPs, 
                         const vector<classNodes *> &nodes) const = 0;
        virtual double getH1Error(const vector<double> &uK, 
                         const vector<classGPs *> &GPs, 
                         const vector<classNodes *> &nodes) const = 0;
                         
        static trueDisplacementField* allocate(const char what[], const std::vector<double>& data);
};


/***
 * 
 * analytic solution of Cantilever beam subjected to a parabolic traction
 * 
 * */
class  CantileverBeamUnderParabolicTraction : public trueDisplacementField
{
    protected:
        double _L; // length
        double _c; // halfwidth
        double _E; // young modulus
        double _nu; // Poisson ratio
        double _P; // applied load
        
    public:
        CantileverBeamUnderParabolicTraction(double L, double c, double E, double nu, double P):
                            _L(L),_c(c),_E(E),_nu(nu),_P(P)
        {
            INFO("analytic solution of a cantilever beam subjected to a parabolic traction\n L=%g c=%g E=%g nu=%g, P=%g",_L,_c,_E,_nu,_P);
        }
        virtual ~CantileverBeamUnderParabolicTraction(){}
        virtual std::string getName() const {return "CantileverBeamUnderParabolicTraction";}
        
        virtual double getL2Error(const vector<double> &uK, 
                         const vector<classGPs *> &GPs, 
                         const vector<classNodes *> &nodes) const;
        virtual double getH1Error(const vector<double> &uK, 
                         const vector<classGPs *> &GPs, 
                         const vector<classNodes *> &nodes) const;
                         
    private:
        void getField(double x, double y, double& u, double& v) const;
        void getGradField(double x, double y, double& dudx, double& dudy, double& dvdx, double& dvdy) const;
};

#endif //_errorEstimation_H_