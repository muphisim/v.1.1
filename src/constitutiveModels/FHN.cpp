//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\FHN.cpp
  \brief This file contains all functions related to the FHN based formulation to describe the evolution of the electrophysiological behaviour at the tissue scale.
*/

#include "FHN.h"
#include "classGPs.h"
#include "TensorOperations.h"


/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] a Model parameter
  @param[in] b Model parameter
  @param[in] epsi Model parameter
  @param[in] vequi Model parameter
  @param[in] N Number of internal variables (r)

*/
classFHN::classFHN(int ndim, string name, string name_cons, double rho, double aelec, double b, double epsi,
                   double vequi, const vector<double> &n0, double D_iso, double D_ani, double stimulusCurrent, 
                   bool mechCouling)
        : constitutiveModels(ndim, name, name_cons, rho),
    _vdotequi(vequi * (vequi - aelec) * (1. - vequi)),
    _aelec(aelec),
    _b(b),
    _epsi(epsi),
    _vequi(vequi),
    _n0(ndim,n0),
    _D_iso(D_iso),
    _D_ani(D_ani),
    _stimulusCurrent(stimulusCurrent),
    _intVarsSize(0),
    _mechCouling(mechCouling)
{
};

classFHN::~classFHN()
{
}

/*! \brief This function initialises this law from options
*/
void classFHN::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    if (_mechCouling)
    {
        _term = new ExtraDofMechanicalFullCouplingTerm(_extraFieldIndexes[0]);
    }
    else
    {
        _term = new ExtraDofTerm(_extraFieldIndexes[0]);
    }
}

double classFHN::soundSpeed() const
{
    double rate = _D_iso + _D_ani;
    return rate;
}

void classFHN::initIntVars(vector<double> &intVars) 
{
    _intVarsSize = intVars.size(); //Index of current FHN internal variable
    intVars.push_back(0); //Recovery current
};

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classFHN::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double>& F = currState.getDeformationGradient();
    double Phi = currState.getExtraFields()[_extraFieldIndexes[0]];
    double Phiprev = prevState.getExtraFields()[_extraFieldIndexes[0]];
    const classTensor1& gradPhi = currState.getExtraFieldGradients()[_extraFieldIndexes[0]];
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    vector<double> &intVarsCurr = currState.getInternalVariables();
    const vector<double> &intVarsprev = prevState.getInternalVariables();
    // update recovery variable
    intVarsCurr[_intVarsSize] = (_epsi * _b * Phi * dt + intVarsprev[_intVarsSize]) / (1. + _epsi * dt);
    int ndim = getDim();
    static classTensor2 D(ndim);
    D.setAll(0);
    static classTensor4 dDdF(ndim);
    static classTensor2 a0a0(ndim);
    TensorOperations::dyadicProduct(_n0,_n0,a0a0);

    computeGL(F, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);


    if (_mechCouling)
    {
        classTensor2 Ftensor(ndim, F); // construct tensor from vector
        if (flagTanMod)
        {
            static classTensor2 C(ndim), Cinv(ndim), DJDF(ndim);
            static classTensor4 DCDF(ndim);
            TensorOperations::rightCauchy(Ftensor,C,&DCDF); // C = FT.F        
            double J = TensorOperations::determinant(Ftensor,&DJDF);
            static classTensor4 DCinvDC(ndim);
            TensorOperations::inverseTensor2(C,Cinv,&DCinvDC,true);
            double invariant4 = TensorOperations::doubleDot(a0a0,C);
            D.axpy(Cinv, _D_iso*J);
            D.axpy(a0a0,J*_D_ani/invariant4);
            static classTensor2 Dinvariant4DF(ndim);
            TensorOperations::doubleDot(a0a0,DCDF,Dinvariant4DF);
            
            // dDdF need to be computed
            TensorOperations::doubleDot(DCinvDC,DCDF,dDdF,_D_iso*J);
            TensorOperations::dyadicProductAdd(Cinv,DJDF,dDdF,_D_iso);
            TensorOperations::dyadicProductAdd(a0a0,DJDF,dDdF,_D_ani/invariant4);
            TensorOperations::dyadicProductAdd(a0a0,Dinvariant4DF,dDdF,-J*_D_ani/invariant4/invariant4);
        }
        else
        {
            static classTensor2 C(ndim), Cinv(ndim);
            TensorOperations::rightCauchy(Ftensor,C); // C = FT.F        
            double J = TensorOperations::determinant(Ftensor);
            TensorOperations::inverseTensor2(C,Cinv);
            double invariant4 = TensorOperations::doubleDot(a0a0,C);
            D.axpy(Cinv, _D_iso*J);
            D.axpy(a0a0,J * _D_ani/invariant4);
        }
    }
    else
    {
       D.addDiagonal(_D_iso);
       D.axpy(a0a0,_D_ani);
    }
        
    classTensor1& flux = currState.getExtraFieldFluxes()[_extraFieldIndexes[0]];
    double& source = currState.getExtraFieldSources()[_extraFieldIndexes[0]];
    double& fieldDensity = currState.getExtraFieldDensities()[_extraFieldIndexes[0]];
    
    fieldDensity = 1;
    flux.setAll(0);
    TensorOperations::dot(D,gradPhi,flux,-1.);
    
    //ionic current
    double Fv = (Phi + _vequi) * (Phi + _vequi - _aelec) * (1 - Phi - _vequi) - _vdotequi - intVarsCurr[_intVarsSize] + _stimulusCurrent;

    // source must be added up to existing one
    source += (Phi-Phiprev)/dt - Fv;
    
    if (flagTanMod)
    {
        double dFvdv = -3. * Phi * Phi + 2. * Phi * (1. - 3. * _vequi + _aelec) - 3 * _vequi * _vequi + 2. * _vequi +
                        2. * _vequi * _aelec - _aelec - _epsi * _b * dt / (1. + _epsi * dt);
    
        double& DsourceDT = (*currState.getDextraFieldSourcesDextraFields())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        classTensor1& DfluxDT = (*currState.getDextraFieldFluxesDextraFields())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        classTensor2& DfluxDgradT = (*currState.getDextraFieldFluxesDextraFieldGrads())[_extraFieldIndexes[0]][_extraFieldIndexes[0]];
        // add-up
        DsourceDT += 1./dt - dFvdv;
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< ndim; j++)
            {
                DfluxDgradT(i,j) -= D(i,j);
            }
        }
        
        if (_mechCouling)
        {
            classTensor3& DfluxDF = (*currState.getDextraFieldFluxesDdeformationGradient())[_extraFieldIndexes[0]];
            for (int i=0; i< ndim; i++)
            {
                for (int j=0; j< ndim; j++)
                {
                    for (int p=0; p< ndim; p++)
                    {
                        for (int q=0; q < ndim; q++)
                        {
                            DfluxDF(i,p,q) -= dDdF(i,j,p,q)*gradPhi(j);
                        }
                    }
                }
            }
            
        }
    }
};
