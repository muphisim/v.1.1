//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#ifndef _classStates_H_
#define _classStates_H_

#include "maths.h"

class classStates
{
    protected:
        int _ndim; /*!< \brief dimension of the problem*/
        int _numMechDof; /*!< \brief number of mechanical dof which is larger than _ndim with stochastic*/
        int _numExtraDof; /*!< \brief number of extraDof*/
        vector<double> _deformationGradient;  /*!< \brief deformation gradient*/
        vector<double> _firstPiolaKirchhoff; /*!< \brief first PK stress*/
        vector<double> _intVars; /*!< \brief internal variables*/
        vector<double> _extraFields; /*!< \brief values of extra-fields*/
        vector<double> _extraFieldSources; /*!< \brief field source*/ 
        vector<classTensor1> _extraFieldGrads; /*!< \brief gradients of extra-fields*/ 
        vector<classTensor1>_extraFieldFluxes; /*!< \brief fluxes of extra-fields*/    
        vector<double> _gl; /*!<\brief GreenLagrangeTensor*/
        vector<double> _VMS; /*!<\brief VMS stress, scalar if deterministic, vector otherwise*/
        vector<double> _cauchy; /*!<\brief Cauchy Tensor */
        vector<double> _equivalentStrain; /*!<\brief Cauchy Tensor */
        vector<double> _volumetricStrain; /*!<\brief Cauchy Tensor */
        vector<double> _extraFieldDensities; /*!< \brief quantity a in the extra-dof equation a*dotT + grad(q) - r= 0*/
        
        // tangents
        bool _tangentAllocated;
        classTensor4 *_tangent;
        classTensor6 *_tangentStochastic;
        // Stress depends on defo and fields
        vector<classTensor2 >* _DfirstPiolaKirchhoffDextraFields;
        // extrafield fluxes depend on fields, field grads, and 
        vector<vector<classTensor1> >*  _DextraFieldFluxesDextraFields;
        vector<vector<classTensor2> >* _DextraFieldFluxesDextraFieldGrads;
        vector<classTensor3>* _DextraFieldFluxesDdeformationGradient;
        // extraField source depends on others fields and deformation
        vector<vector<double> >* _DextraFieldSourcesDextraFields;
        vector<classTensor2>* _DextraFieldSourcesDdeformationGradient;
        
    public:
        classStates(int ndim); /*!< \brief constructor*/
        classStates(const classStates& src);  /*!< \brief copy constructor*/
        classStates& operator = (const classStates& src); /*!< \brief equal operator*/
        virtual ~classStates();
        
        void startFromOtherState(const classStates& src);
        
        int getDimension() const {return _ndim;};
        int getNumMechanicalDofsPerNode() const {return _numMechDof;};
        int getNumExtraDofsPerNode() const {return _numExtraDof;};
        int getTotalDofsPerNode() const {return _numMechDof+_numExtraDof;};
        
        void allocateMechanicalField(int numMechDof);
        void allocateExtraField(int numExtraDof);
        
        bool tangentAllocated() const {return _tangentAllocated;};
        void allocateTangent();
        
        vector<double> & getDeformationGradient() {return _deformationGradient;};
        const vector<double> & getDeformationGradient() const {return _deformationGradient;};
        
        vector<double> & getFirstPiolaKirchhoffStress() {return _firstPiolaKirchhoff;};
        const vector<double> & getFirstPiolaKirchhoffStress() const {return _firstPiolaKirchhoff;};
        
        vector<double> & getInternalVariables() {return _intVars;};
        const vector<double> & getInternalVariables() const {return _intVars;};
        
        vector<double>& getExtraFields() {return _extraFields;};
        const vector<double>& getExtraFields() const {return _extraFields;};
        
        vector<double>& getExtraFieldSources() {return _extraFieldSources;};
        const vector<double>& getExtraFieldSources() const {return _extraFieldSources;};
        
        vector<classTensor1>& getExtraFieldGradients() {return _extraFieldGrads;};
        const vector<classTensor1>& getExtraFieldGradients() const {return _extraFieldGrads;};
        
        vector<classTensor1>& getExtraFieldFluxes() {return _extraFieldFluxes;};
        const vector<classTensor1>& getExtraFieldFluxes() const {return _extraFieldFluxes;};


        vector<double>& getGL(){return _gl;};
        vector<double>& getVMS(){return _VMS;};
        vector<double>& getCauchy(){return _cauchy;};
        vector<double>& getEquivalentStrain(){return _equivalentStrain;};
        vector<double>& getVolumetricStrain(){return _volumetricStrain;};



        vector<double>& getExtraFieldDensities() {return _extraFieldDensities;};
        const vector<double>& getExtraFieldDensities() const {return _extraFieldDensities;};
                
        void getCauchy(vector<double>& Cauchy) const;
        
        classTensor4 *getTangent() {return _tangent;};
        const classTensor4 *getTangent() const {return _tangent;};
        
        classTensor6 *getTangentStochastic() {return _tangentStochastic;};
        const classTensor6 *getTangentStochastic() const {return _tangentStochastic;};
        
        vector<classTensor2 >* getDfirstPiolaKirchhoffDextraFields() {return _DfirstPiolaKirchhoffDextraFields;};
        const vector<classTensor2 >* getDfirstPiolaKirchhoffDextraFields() const {return _DfirstPiolaKirchhoffDextraFields;};

        vector<vector<classTensor1> >* getDextraFieldFluxesDextraFields() {return _DextraFieldFluxesDextraFields;}
        const vector<vector<classTensor1> >* getDextraFieldFluxesDextraFields() const {return _DextraFieldFluxesDextraFields;}
        
        vector<vector<classTensor2> >* getDextraFieldFluxesDextraFieldGrads() {return _DextraFieldFluxesDextraFieldGrads;}
        const vector<vector<classTensor2> >* getDextraFieldFluxesDextraFieldGrads() const {return _DextraFieldFluxesDextraFieldGrads;}
        
        vector<classTensor3>* getDextraFieldFluxesDdeformationGradient() {return _DextraFieldFluxesDdeformationGradient;};
        const vector<classTensor3>* getDextraFieldFluxesDdeformationGradient() const {return _DextraFieldFluxesDdeformationGradient;};
        
        // extraField source depends on others fields and deformation
        vector<vector<double> >* getDextraFieldSourcesDextraFields() {return _DextraFieldSourcesDextraFields;};
        const vector<vector<double> >* getDextraFieldSourcesDextraFields() const {return _DextraFieldSourcesDextraFields;};
        
        vector<classTensor2>* getDextraFieldSourcesDdeformationGradient() {return _DextraFieldSourcesDdeformationGradient;};
        const vector<classTensor2>* getDextraFieldSourcesDdeformationGradient() const {return _DextraFieldSourcesDdeformationGradient;};
        
};
#endif //_classStates_H_
