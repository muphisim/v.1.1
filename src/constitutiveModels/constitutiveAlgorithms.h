//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _CONSTITUTIVEALGORITHMS_H_
#define _CONSTITUTIVEALGORITHMS_H_

#include "configuration.h"
#include "maths.h"

/*! \brief  This is the wrapper to get forces and stiffness, mass matrices at each Gausspoint. 
 * This base class implements general and virtual functions. 
 * ALL VIRTUAL FUNCTIONS MUST BE IMPLEMENTE BY THE USER */
class Term
{
    public:
        virtual ~Term(){};
        /*! \brief This function is to get the index of Dofs at the Gauss point in the local vector
          @param[in] GaussPoint Gauss Point
          @param[out] dofPosition Array of all indexes of the related dofs define by terms
        */
        virtual void getDofsLocalPosition(const classGPs *GaussPoint, vector<int>& dofPosition) const = 0;
        /*! \brief This function to get the internal force
          @param[in] GaussPoint Gauss Point
          @param[out] vForce force vector
        */
        virtual void getForceVector(const classGPs *GaussPoint, mVector& vForce) const = 0;
        /*! \brief This function to get the stiffness matrix
          @param[in] GaussPoint Gauss Point
          @param[out] stiff stiffness matrix
        */
        virtual void getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const = 0;
        /*! \brief This function to get the mass matrix
          @param[in] GaussPoint Gauss Point
          @param[out] mass mass matrix
        */
        virtual void getMassMatrix(const classGPs *GaussPoint, mMatrix& mass) const =0;
        
        /*! \brief This function to get the lumped mass vector
          @param[in] GaussPoint Gauss Point
          @param[out] mass mass vector
        */
        virtual void getLumpedMassVector(const classGPs *GaussPoint, mVector& mass) const = 0;
        /*! \brief This function to to know if coupling ocurs
          If the coupling occurs, the stiffness matrix has the number of cols equal to number of existing dofs in GP
          otherwise the stiffness matrix has a square form of the related field
        */
        virtual bool withCoupling() const = 0;
        
        /*! \brief This function to get fNint in explicit scheme with extraDof
          @param[in] GaussPoint Gauss Point
          @param[out] vForce force vector
          @param[in] timeRun current time
          @param[in] dt time step
        */
        virtual void getForceVectorExplicitScheme(const classGPs *GaussPoint, mVector& vForce, double timeRun, double dt) const = 0;
};

/*! \brief  This is the wrapper for pure mechanical response. 
*/
class MechanicalTerm : public Term
{
    public:
        MechanicalTerm(){}
        virtual ~MechanicalTerm(){}
        virtual void getDofsLocalPosition(const classGPs *GaussPoint, vector<int>& dofPosition) const;
        virtual void getForceVector(const classGPs *GaussPoint, mVector& vForce) const;
        virtual void getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const;
        virtual void getMassMatrix(const classGPs *GaussPoint, mMatrix& mass) const;
        virtual void getLumpedMassVector(const classGPs *GaussPoint, mVector& mass) const;
        virtual void getForceVectorExplicitScheme(const classGPs *GaussPoint, mVector& vForce, double timeRun, double dt) const;
        virtual bool withCoupling() const {return false;}
};

/*! \brief  This is the wrapper for pure stochastic mechanical response. 
*/
class MechanicalTermStochastic : public Term
{
    public:
        MechanicalTermStochastic(){}
        virtual ~MechanicalTermStochastic(){}
        virtual void getDofsLocalPosition(const classGPs *GaussPoint, vector<int>& dofPosition) const;
        virtual void getForceVector(const classGPs *GaussPoint, mVector& vForce) const;
        virtual void getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const;
        virtual void getMassMatrix(const classGPs *GaussPoint, mMatrix& stiff) const;
        virtual void getLumpedMassVector(const classGPs *GaussPoint, mVector& mass) const;
        virtual void getForceVectorExplicitScheme(const classGPs *GaussPoint, mVector& vForce, double timeRun, double dt) const;
        virtual bool withCoupling() const {return false;}
};

/*! \brief  This is the wrapper for pure extraDof response. 
*/
class ExtraDofTerm : public Term
{
    protected:
        int _fieldIndex;    
            
    public:
        ExtraDofTerm(int fieldIndex);
        virtual ~ExtraDofTerm(){}
        virtual void getDofsLocalPosition(const classGPs *GaussPoint, vector<int>& dofPosition) const;
        virtual void getForceVector(const classGPs *GaussPoint, mVector& vForce) const;
        virtual void getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const;
        virtual void getMassMatrix(const classGPs *GaussPoint, mMatrix& stiff) const;
        virtual void getLumpedMassVector(const classGPs *GaussPoint, mVector& mass) const;
        virtual void getForceVectorExplicitScheme(const classGPs *GaussPoint, mVector& vForce, double timeRun, double dt) const;
        virtual bool withCoupling() const {return false;}
};


/*! \brief  This is the wrapper for coupling mechanical response coupled extraDofs response. 
*/
class  MechanicalExtraDofFullCouplingTerm: public MechanicalTerm
{
    protected:
        vector<int> _fieldIndexes;    
            
    public:
        MechanicalExtraDofFullCouplingTerm(const vector<int>& fieldIndexes);
        virtual ~MechanicalExtraDofFullCouplingTerm(){}
        virtual void getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const;
        virtual bool withCoupling() const {return true;}
};


/*! \brief  This is the wrapper for coupling  extraDofs - mechanical response. 
*/
class  ExtraDofMechanicalFullCouplingTerm: public ExtraDofTerm
{
    public:
        ExtraDofMechanicalFullCouplingTerm(int fieldIndex);
        virtual ~ExtraDofMechanicalFullCouplingTerm(){}
        virtual void getStiffnessMatrix(const classGPs *GaussPoint, mMatrix& stiff) const;
        virtual bool withCoupling() const {return true;}
};

#endif //_CONSTITUTIVEALGORITHMS_H_