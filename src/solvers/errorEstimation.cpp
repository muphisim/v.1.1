//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file errorEstimation.cpp
  \brief This file contains all functions related the error analysis
*/

#include "errorEstimation.h"
#include "classGPs.h"

trueDisplacementField* trueDisplacementField::allocate(const char what[], const std::vector<double>& data)
{
    INFO("error estimation for %s",what);
    printVector(data,"input data");
    if (strcmp(what,"CantileverBeamUnderParabolicTraction") == 0 )
    {
        if (data.size() <2 )
        {
            ERROR("The number of parameters must be equal to 2 with CantileverBeamUnderParabolicTraction");
            exit(-1);
        }
        return new CantileverBeamUnderParabolicTraction(data[0],data[1],data[2],data[3],data[4]);
    }
    else
    {
        ERROR("trueDisplacementField cannot be allocated with %s",what);
        exit(-1);
    }
};


void CantileverBeamUnderParabolicTraction::getField(double x, double y, double& u, double& v) const
{
    double I = 2.*_c*_c*_c/3.;
    u = (_P*y/(6*_E*I))*((6*_L-3*x)*x+(2+_nu)*(y*y-_c*_c));
    v = -(_P/(6*_E*I))*(3*_nu*y*y*(_L-x)+(4+5*_nu)*_c*_c*x+(3*_L-x)*x*x);
};

void CantileverBeamUnderParabolicTraction::getGradField(double x, double y, double& dudx, double& dudy, double& dvdx, double& dvdy) const
{
    double I = 2.*_c*_c*_c/3.;
    dudx = (_P*y/(6*_E*I))*((-3)*x+(6*_L-3*x));;
    dudy = (_P/(6*_E*I))*((6*_L-3*x)*x+(2+_nu)*(y*y-_c*_c)) + (_P*y/(6*_E*I))*((2+_nu)*(2*y));
    dvdx = -(_P/(6*_E*I))*(3*_nu*y*y*(-1)+(4+5*_nu)*_c*_c+(-1)*x*x+(3*_L-x)*2*x);
    dvdy = -(_P/(6*_E*I))*(3*_nu*2*y*(_L-x));
};

double CantileverBeamUnderParabolicTraction::getL2Error(const vector<double> &uK, 
                         const vector<classGPs *> &GPs, 
                         const vector<classNodes *> &nodes) const
{
    double error = 0;
    double Unorm = 0;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        const classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        const vector<double>& xyz = GaussPoint->getXYZ();
        const vector<double>& phi = GaussPoint->getPhi();
        const vector<int>& neighbours = GaussPoint->getNeighbours();
        double weightJ = GaussPoint->getWeightJ();
        int numNodes = nodes.size();
        double ucur=0.; 
        double vcur=0.;
        int numNeighbours = neighbours.size();
        for (int a = 0; a < numNeighbours; a++) 
        {
            // FE approximation
            ucur += uK[neighbours[a]+ 0 * numNodes] * phi[a];
            vcur += uK[neighbours[a]+ 1 * numNodes] * phi[a];
        }
        double uAnalytic=0;
        double vAnalytic=0;
        getField(xyz[0],xyz[1],uAnalytic,vAnalytic);
        //INFO("x = %g y = %g ucur = %g vcur = %g ua = %g va = %g",xyz[0],xyz[1],ucur,vcur,uAnalytic,vAnalytic);
        error += weightJ*((ucur - uAnalytic)*(ucur - uAnalytic) + (vcur-vAnalytic)*(vcur-vAnalytic));
        Unorm += weightJ*((uAnalytic )*(uAnalytic) + (vAnalytic)*(vAnalytic));
    }
    return sqrt(error/Unorm);
};

double CantileverBeamUnderParabolicTraction::getH1Error(const vector<double> &uK, 
                 const vector<classGPs *> &GPs, 
                 const vector<classNodes *> &nodes) const
{
 
    double error = 0;
    double Enorm = 0;
    for (vector<classGPs *>::const_iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        const classGPs *GaussPoint = *ite;
        if ((GaussPoint->getActivate()) == 0)
            continue;
        const vector<double>& xyz = GaussPoint->getXYZ();
        const vector<double>& phi = GaussPoint->getPhi();
        const vector<double>& Dphi = GaussPoint->getDphi();
        const vector<int>& neighbours = GaussPoint->getNeighbours();
        double weightJ = GaussPoint->getWeightJ();
        int numNodes = nodes.size();
        int ndim = GaussPoint->getDimension();
        double ucur(0.), DuDx(0.), DuDy(0.); 
        double vcur(0.), DvDx(0.), DvDy(0.);
        int numNeighbours = neighbours.size();
        for (int a = 0; a < numNeighbours; a++) 
        {
            // FE approximation
            ucur += uK[neighbours[a]+ 0 * numNodes] * phi[a];
            
            DuDx +=  uK[neighbours[a]+ 0 * numNodes]*Dphi[a * ndim + 0];
            DuDy +=  uK[neighbours[a]+ 0 * numNodes]*Dphi[a * ndim + 1];
            
             vcur += uK[neighbours[a]+ 1 * numNodes] * phi[a];    
            DvDx +=  uK[neighbours[a]+ 1 * numNodes]*Dphi[a * ndim + 0];
            DvDy +=  uK[neighbours[a]+ 1 * numNodes]*Dphi[a * ndim + 1];
            
        }
        double uAnalytic=0;
        double vAnalytic=0;
        getField(xyz[0],xyz[1],uAnalytic,vAnalytic);
        
        double DuAnalyticDx, DuAnalyticDy, DvAnalyticDx, DvAnalyticDy;
        getGradField(xyz[0],xyz[1],DuAnalyticDx, DuAnalyticDy, DvAnalyticDx, DvAnalyticDy);
        
        error += weightJ*((ucur - uAnalytic)*(ucur - uAnalytic) + 
                          (vcur-vAnalytic)*(vcur-vAnalytic)+
                          _L*_L*(DuDx-DuAnalyticDx)*(DuDx-DuAnalyticDx)+
                          _L*_L*(DuDy-DuAnalyticDy)*(DuDy-DuAnalyticDy)+
                          _L*_L*(DvDx-DvAnalyticDx)*(DvDx-DvAnalyticDx)+
                          _L*_L*(DvDy-DvAnalyticDy)*(DvDy-DvAnalyticDy));
        Enorm += weightJ*((uAnalytic)*(uAnalytic) + 
                          (vAnalytic)*(vAnalytic)+
                          _L*_L*(DuAnalyticDx)*(DuAnalyticDx)+
                          _L*_L*(DuAnalyticDy)*(DuAnalyticDy)+
                          _L*_L*(DvAnalyticDx)*(DvAnalyticDx)+
                          _L*_L*(DvAnalyticDy)*(DvAnalyticDy));
        //INFO("error = %g",error);
    }
    return sqrt(error/Enorm);
};