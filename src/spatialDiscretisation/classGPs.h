//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classGPs_H_
#define _classGPs_H_

#include "configuration.h"
#include "classSpatialPoint.h"
#include "constitutiveModels.h"
#include "maths.h"
#include "classStates.h"
#include "constitutiveManager.h"

/*! \brief  This class defines each GP in the domain*/
class classGPs : public classSpatialPoint {

protected:
    string _eleType; /*!< \brief String with the type of element in which this GP belongs to*/
    int cell;
    int bulkElement;
    double w; /*!< \brief Weight in the quadrature*/
    double J; /*!< \brief Jacobian of the quadrature*/
    
    double _wJ;  /*!< \brief Weight*Jacobian at the quadrature*/
    constitutiveManager _allaws;  /*!< \brief all material laws at GPs*/
    classStates _currentState; /*!< \brief current state*/
    classStates _previousState; /*!< \brief previous state*/
    classStates _initialState; /*!< \brief initial state*/
    
    classStates* _tmpState; // for checkpoint purpose;
    
    int activate; /*! \brief =1 means is active, =0 means is inactive*/
    int flagFIni;

    
public:

    classGPs(int me, const vector<double>& xyz, double weight, double Jacobian, int ndim);
    
    inline double getW() const { return w; }
    inline double getJ() const { return J; }
    
    inline int getActivate() const { return activate; };
    inline void setActivate(int _activate) { this->activate = _activate; };
    
    double getWeightJ() const { return _wJ; };// Common
    constitutiveManager& getConstitutiveManager() {return _allaws;};
    const constitutiveManager& getConstitutiveManager() const {return _allaws;};
    
    classStates& getCurrentState() {return _currentState;};
    const classStates& getCurrentState() const {return _currentState;};
    
    classStates& getPreviousState() {return _previousState;};
    const classStates& getPreviousState() const {return _previousState;};
    
    classStates& getInitialState() {return _initialState;};
    const classStates& getInitialState() const {return _initialState;};
    
    void getLocalDofs(int numNodes, std::vector<int>& vdofs) const;
    
    void nextStep() 
    { 
        _previousState= _currentState;
    };
    void resetToPreviousState()
    {
        _currentState = _previousState;
    };
    void resetToInitialState()
    {
        _previousState = _initialState;
        _currentState = _initialState;
    };

    inline vector<double> getCauchy() const {
        vector<double> Cauchy;
        _currentState.getCauchy(Cauchy);
        return Cauchy;
    };


    inline void setflagFIni() {
        this->flagFIni = 1;
    };

    inline int getflagFIni() const{
        return this->flagFIni;
    };

    inline void setBulkElement(int _bulkElement) {
        this->bulkElement = _bulkElement;
    };

    inline void setCell(int cell, string eleType) {
        this->cell = cell;
        this->_eleType = eleType;
    };
    inline string getElementType() const {return _eleType;};

    void setConsModel(vector<constitutiveModels *>& constitutiveModel);
    
    inline int getCell() const { return cell; };// Common
    inline int getBulkElement() const { return bulkElement; };
    
    void makeACheckPoint();
    void restoreACheckPoint();

    /*! Destructor
     */
    virtual ~classGPs() {
        if (_tmpState) delete _tmpState;
    };
};

#endif

