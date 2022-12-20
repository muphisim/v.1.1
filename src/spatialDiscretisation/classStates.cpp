//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
#include "classStates.h"

/*! Constructor for classState
  @param[in] ndim dimension of the problem
  @param[in] numExtraDof number of extraDof
  @param[in] nstoch number of stochastic dofs except the mean ones
*/
classStates::classStates(int ndim):
    _ndim(ndim), _numMechDof(_ndim), _numExtraDof(0),
    _deformationGradient(ndim * ndim, 0.),
    _firstPiolaKirchhoff(ndim * ndim, 0.),
    _cauchy(ndim*ndim,0),
    _gl(ndim * ndim, 0),
    _VMS(1,0),
    _equivalentStrain(1,0),
    _volumetricStrain(1,0),
    _intVars(0),
    _extraFields(0),
    _extraFieldSources(0),
    _extraFieldGrads(0),
    _extraFieldFluxes(0),
    _extraFieldDensities(0),
    _tangentAllocated(false),
    _tangent(NULL),
    _tangentStochastic(NULL),
    _DfirstPiolaKirchhoffDextraFields(NULL),
    _DextraFieldFluxesDextraFields(NULL),
    _DextraFieldFluxesDextraFieldGrads(NULL),
    _DextraFieldFluxesDdeformationGradient(NULL),
    _DextraFieldSourcesDextraFields(NULL),
    _DextraFieldSourcesDdeformationGradient(NULL)
{
    for (int i = 0; i < ndim; i++) 
    {
        _deformationGradient[i * ndim + i] = 1.0; // diagonal one
    }
}
classStates::classStates(const classStates& src):
    _ndim(src._ndim), _numMechDof(src._numMechDof), _numExtraDof(src._numExtraDof),
    _deformationGradient(src._deformationGradient),
    _firstPiolaKirchhoff(src._firstPiolaKirchhoff),
    _intVars(src._intVars),
    _extraFields(src._extraFields),
    _extraFieldSources(src._extraFieldSources),
    _extraFieldGrads(src._extraFieldGrads),
    _extraFieldFluxes(src._extraFieldFluxes), 
    _extraFieldDensities(src._extraFieldDensities),
    _tangentAllocated(src._tangentAllocated),
    _tangent(NULL),
    _tangentStochastic(NULL),
    _DfirstPiolaKirchhoffDextraFields(NULL),
    _DextraFieldFluxesDextraFields(NULL),
    _DextraFieldFluxesDextraFieldGrads(NULL),
    _DextraFieldFluxesDdeformationGradient(NULL),
    _DextraFieldSourcesDextraFields(NULL),
    _DextraFieldSourcesDdeformationGradient(NULL)
{
    if (src._tangent) 
        _tangent = new classTensor4(*src._tangent);;
    if (src._tangentStochastic) 
        _tangentStochastic = new classTensor6(*src._tangentStochastic);
    if (src._DfirstPiolaKirchhoffDextraFields) 
        _DfirstPiolaKirchhoffDextraFields = new vector<classTensor2 >(*src._DfirstPiolaKirchhoffDextraFields);
    if (src._DextraFieldFluxesDextraFields)
        _DextraFieldFluxesDextraFields = new vector<vector<classTensor1> >(*src._DextraFieldFluxesDextraFields);
    if (src._DextraFieldFluxesDextraFieldGrads)
        _DextraFieldFluxesDextraFieldGrads = new vector<vector<classTensor2> >(*src._DextraFieldFluxesDextraFieldGrads);
    if (src._DextraFieldFluxesDdeformationGradient)
        _DextraFieldFluxesDdeformationGradient = new vector<classTensor3>(*src._DextraFieldFluxesDdeformationGradient);
    if (src._DextraFieldSourcesDextraFields) 
        _DextraFieldSourcesDextraFields = new vector<vector<double> >(*src._DextraFieldSourcesDextraFields);
    if (src._DextraFieldSourcesDdeformationGradient) 
        _DextraFieldSourcesDdeformationGradient = new vector<classTensor2>(*src._DextraFieldSourcesDdeformationGradient);
}
classStates& classStates::operator = (const classStates& src)
{
    _ndim = (src._ndim);
    _numMechDof = (src._numMechDof);
    _numExtraDof = (src._numExtraDof);
    _deformationGradient = (src._deformationGradient);
    _firstPiolaKirchhoff = (src._firstPiolaKirchhoff);
    _intVars = (src._intVars);
    _extraFields = (src._extraFields);
    _extraFieldSources = (src._extraFieldSources);
    _extraFieldGrads = (src._extraFieldGrads);
    _extraFieldFluxes = (src._extraFieldFluxes);
    _extraFieldDensities = (src._extraFieldDensities);
    _tangentAllocated = src._tangentAllocated;
    if (src._tangent)
    { 
        if (_tangent)
            (*_tangent) = (*src._tangent);
        else
            _tangent = new classTensor4(*src._tangent);
    }
    else
    {
        if (_tangent) delete _tangent;
        _tangent = NULL;
    }
    if (src._tangentStochastic)
    {
        if (_tangentStochastic)
            (*_tangentStochastic) = (*src._tangentStochastic);
        else
            _tangentStochastic = new classTensor6(*src._tangentStochastic);
    } 
    else
    {
        if (_tangentStochastic) delete _tangentStochastic;
        _tangentStochastic = NULL;
    }
        
    if (src._DfirstPiolaKirchhoffDextraFields) 
    {
        if (_DfirstPiolaKirchhoffDextraFields)
            (*_DfirstPiolaKirchhoffDextraFields) = (*src._DfirstPiolaKirchhoffDextraFields);
        else
            _DfirstPiolaKirchhoffDextraFields = new vector<classTensor2 >(*src._DfirstPiolaKirchhoffDextraFields);
    }
    else
    {
        if (_DfirstPiolaKirchhoffDextraFields) delete _DfirstPiolaKirchhoffDextraFields;
        _DfirstPiolaKirchhoffDextraFields = NULL;
    }
    
    if (src._DextraFieldFluxesDextraFields)
    {
        if (_DextraFieldFluxesDextraFields)
            (*_DextraFieldFluxesDextraFields) = (*src._DextraFieldFluxesDextraFields);
        else
            _DextraFieldFluxesDextraFields = new vector<vector<classTensor1> >(*src._DextraFieldFluxesDextraFields);
    }
    else
    {
        if (_DextraFieldFluxesDextraFields) delete _DextraFieldFluxesDextraFields;
        _DextraFieldFluxesDextraFields = NULL;
    }
    
    if (src._DextraFieldFluxesDextraFieldGrads)
    {
        if (_DextraFieldFluxesDextraFieldGrads)
            (*_DextraFieldFluxesDextraFieldGrads) = (*src._DextraFieldFluxesDextraFieldGrads);
        else
            _DextraFieldFluxesDextraFieldGrads = new vector<vector<classTensor2> >(*src._DextraFieldFluxesDextraFieldGrads);
    }
    else
    {
        if (_DextraFieldFluxesDextraFieldGrads) delete _DextraFieldFluxesDextraFieldGrads;
        _DextraFieldFluxesDextraFieldGrads = NULL;
    }
    
    if (src._DextraFieldFluxesDdeformationGradient)
    {
        if (_DextraFieldFluxesDdeformationGradient)
            (*_DextraFieldFluxesDdeformationGradient) = (*src._DextraFieldFluxesDdeformationGradient);
        else
            _DextraFieldFluxesDdeformationGradient = new vector<classTensor3>(*src._DextraFieldFluxesDdeformationGradient);
    }
    else
    {
        if (_DextraFieldFluxesDdeformationGradient) delete _DextraFieldFluxesDdeformationGradient;
        _DextraFieldFluxesDdeformationGradient = NULL;
    }
    
    if (src._DextraFieldSourcesDextraFields) 
    {
        if (_DextraFieldSourcesDextraFields)
            (*_DextraFieldSourcesDextraFields) = (*src._DextraFieldSourcesDextraFields);
        else
            _DextraFieldSourcesDextraFields = new vector<vector<double> >(*src._DextraFieldSourcesDextraFields);
    }
    else
    {
        if (_DextraFieldSourcesDextraFields) delete _DextraFieldSourcesDextraFields;
        _DextraFieldSourcesDextraFields = NULL;
    }
    
    if (src._DextraFieldSourcesDdeformationGradient) 
    {
        if (_DextraFieldSourcesDdeformationGradient)
            (*_DextraFieldSourcesDdeformationGradient) = (*src._DextraFieldSourcesDdeformationGradient);
        else
            _DextraFieldSourcesDdeformationGradient = new vector<classTensor2>(*src._DextraFieldSourcesDdeformationGradient);
    }   
    else
    {
        if (_DextraFieldSourcesDdeformationGradient) delete _DextraFieldSourcesDdeformationGradient;
        _DextraFieldSourcesDdeformationGradient = NULL;
    }
    return *this;
};
classStates::~classStates()
{
    if (_tangent) delete _tangent;
    if (_tangentStochastic) delete _tangentStochastic;
    if (_DfirstPiolaKirchhoffDextraFields) delete _DfirstPiolaKirchhoffDextraFields;
    if (_DextraFieldFluxesDextraFields) delete _DextraFieldFluxesDextraFields;
    if (_DextraFieldFluxesDextraFieldGrads) delete _DextraFieldFluxesDextraFieldGrads;
    if (_DextraFieldFluxesDdeformationGradient) delete _DextraFieldFluxesDdeformationGradient;
    if (_DextraFieldSourcesDextraFields) delete _DextraFieldSourcesDextraFields;
    if (_DextraFieldSourcesDdeformationGradient) delete _DextraFieldSourcesDdeformationGradient;
    
}

void classStates::startFromOtherState(const classStates& src)
{
    _intVars = (src._intVars);
    setAll(_firstPiolaKirchhoff,0);
    for (int i=0; i < _numExtraDof; i++)
    {
        _extraFieldFluxes[i].setAll(0);
    }
    setAll(_extraFieldSources,0.);
    setAll(_extraFieldDensities,0.);
    
    if (_tangentAllocated)
    {
        if (_tangent) _tangent->setAll(0);
        if (_tangentStochastic) _tangentStochastic->setAll(0.);
        
        for (int j=0; j< _numExtraDof; j++)
        {
            (*_DfirstPiolaKirchhoffDextraFields)[j].setAll(0.);
            (*_DextraFieldFluxesDdeformationGradient)[j].setAll(0.);
            (*_DextraFieldSourcesDdeformationGradient)[j].setAll(0.);
            
            for (int k=0; k < _numExtraDof; k++)
            {
                (*_DextraFieldFluxesDextraFields)[j][k].setAll(0.);
                (*_DextraFieldFluxesDextraFieldGrads)[j][k].setAll(0.);
                (*_DextraFieldSourcesDextraFields)[j][k]= 0.;
                
            }
        }
        
    }
}

void classStates::allocateMechanicalField(int numMechDof)
{
    _numMechDof = numMechDof;
    int nstochField = (numMechDof)/_ndim;
    _deformationGradient.clear();
    _deformationGradient.resize((nstochField) * _ndim * _ndim, 0.);
    for (int i = 0; i < _ndim; i++) 
    {
        _deformationGradient[i * _ndim + i] = 1.0; // diagonal one
    }
    _firstPiolaKirchhoff.clear();
    _firstPiolaKirchhoff.resize((nstochField) * _ndim * _ndim, 0.);
    _gl.clear();
    _gl.resize((nstochField) * _ndim * _ndim, 0.);
    _VMS.clear();
    _VMS.resize(nstochField, 0.);
    _cauchy.clear();
    _cauchy.resize((nstochField) * _ndim * _ndim, 0.);
    _equivalentStrain.clear();
    _equivalentStrain.resize(nstochField, 0.);
    _volumetricStrain.clear();
    _volumetricStrain.resize(nstochField, 0.);
}
void classStates::allocateExtraField(int numExtraDof)
{
    _numExtraDof = numExtraDof;
    _extraFields.clear();
    _extraFieldSources.clear();
    _extraFieldGrads.clear();
    _extraFieldFluxes.clear();
    
    _extraFields.resize(_numExtraDof,0);
    _extraFieldSources.resize(_numExtraDof,0);
    _extraFieldGrads.resize(_numExtraDof,classTensor1(_ndim));
    _extraFieldFluxes.resize(_numExtraDof,classTensor1(_ndim));
    _extraFieldDensities.resize(_numExtraDof,0.);
};

void classStates::allocateTangent()
{
    if (_tangentAllocated) return;
    _tangentAllocated = true;
    
    _tangent = new classTensor4(_ndim);
    int nStochastic = _numMechDof/_ndim;
    if (nStochastic >1)
    {
        _tangentStochastic = new classTensor6(_ndim,nStochastic);
    }
    else
    {
        _tangentStochastic = NULL;
    }
    _DfirstPiolaKirchhoffDextraFields = new vector<classTensor2 >(_numExtraDof,classTensor2(_ndim));
    _DextraFieldFluxesDextraFields = new vector<vector<classTensor1> >(_numExtraDof, vector<classTensor1>(_numExtraDof,classTensor1(_ndim)));
    _DextraFieldFluxesDextraFieldGrads = new vector<vector<classTensor2> >(_numExtraDof, vector<classTensor2>(_numExtraDof,classTensor2(_ndim)));
    _DextraFieldFluxesDdeformationGradient = new vector<classTensor3>(_numExtraDof,classTensor3(_ndim));
    _DextraFieldSourcesDextraFields = new vector<vector<double> >(_numExtraDof,std::vector<double>(_numExtraDof,0.));
    _DextraFieldSourcesDdeformationGradient = new vector<classTensor2>(_numExtraDof,classTensor2(_ndim));
};


void classStates::getCauchy(vector<double>& Cauchy) const 
{
    setAll(Cauchy,0);
    Cauchy.resize(_ndim *_ndim, 0);
    double detJinv = 1./determinantTensor(_deformationGradient, _ndim);
    multSTensor3SecondTranspose(_firstPiolaKirchhoff, _ndim, _ndim, _deformationGradient, _ndim, Cauchy);
    multiSTensorScalar(Cauchy,_ndim,_ndim,detJinv);
};
