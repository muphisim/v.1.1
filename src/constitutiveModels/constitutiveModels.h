//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _constitutiveModels_H_
#define _constitutiveModels_H_

#include "configuration.h"
#include "constitutiveAlgorithms.h"

/*! \brief  This is the wrapper for all constitutive models. This base class implements general and virtual functions. ALL VIRTUAL FUNCTIONS MUST BE IMPLEMENTE BY THE USER */
class constitutiveModels 
{
    protected:
        string _name; // name given in the input filed
        string _name_cons; // name of the constitutive lew
        int _ndim; // dimension of the problem
        double _rho; // density
        vector<int> _extraFieldIndexes; // number of extraFields in this constitutive law, zero vector for mechancis, an extra-filed can contain multiple dofs
        bool _isInitialised; // true if some thing needs to be set from options
        Term* _term;
    
    public:
        constitutiveModels(int ndim, string name, string name_cons, double rho);
        virtual ~constitutiveModels();
        inline string getName() const{ return _name; };
        inline string getNameCons() const { return _name_cons; };
        inline double getRho() const { return _rho; };
        inline int getDim() const {return _ndim;};
        inline void addExtraField(int field) {_extraFieldIndexes.push_back(field);};
        inline void clearExtraFieldIndexes() {_extraFieldIndexes.clear();};
        inline const vector<int>& getExtraFieldIndexes() const {return _extraFieldIndexes;}; 
        inline bool isInitialised() const {return _isInitialised;};
        const Term* getTerm() const {return _term;}
        
        /*! \brief This function initialises the law form options, 
         * eg - allocate term 
        */
        virtual void initLawFromOptions() = 0;
        
        /*! \brief This function returns the number of unknows at each nodes. 
         * This function is mandatory and it must be implemented by the user.
          @return number of dofs at each node
        */
        virtual int getNbrDofConsMod() const = 0; 

        /*! \brief This function provides the sounds velocity in the material to calculate the critical time step. 
         * This function is mandatory and it must be implemented by the user.
          @return Sound speed in the material
        */
        virtual double soundSpeed() const = 0;

        /*! \brief This function initialise the internal variables this constitutive model. 
         * This function is mandatory and it must be implemented by the user.
          @param[inout] intVars Array of the internal variables
        */
        virtual void initIntVars(vector<double> &intVars) = 0;
        
        /*! \brief This function checks activation of a GP associated associated with this constitutive model. 
         * This function is mandatory and it must be implemented by the user.
         @param[inout] GaussPoint Gauss point
         @param[in] timeRun Current time
         */
        virtual void checkActivation(classGPs *GaussPoint, double timeRun) const =0;
        
        /*! \brief This function predict the current state a GP before calling updateConstitutive. 
         * Eg, check phase change and update initial deformation in 3D printing model 
         * This function is mandatory and it must be implemented by the user.
          @param[inout] GaussPoint Gauss point
          @param[in] dt Time step
          @param[in] timeRun Current time
         */
        virtual void predictIntVars(classGPs *GaussPoint, double dt, double timeRun) const =0;
        
        /*! \brief This function updates the constitutive response. 
         * This function is mandatory and it must be implemented by the user.
          @param[inout] GaussPoint Gauss point
          @param[in] dt Time step
          @param[in] timeRun Current time
          @param[in] flagTanMod Flag to calculate the tangent modulus
        */
        virtual void updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) =0;
        
        /*! \brief This function return the elastic potential enery, 0 by defaut. 
          @param[inout] GaussPoint Gauss point
        */
        virtual double defoEnergy(classGPs *GaussPoint) const 
        {
          return 0.;
        }

};

#endif
