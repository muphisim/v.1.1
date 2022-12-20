//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classNodes_H
#define _classNodes_H

#include "configuration.h"
#include "classSpatialPoint.h"
#include "constitutiveModels.h"

/*! \brief  This class defines each node in the domain*/
class classNodes : public classSpatialPoint {

protected:
    double tol = 1E-10; /*!< \brief This tol is quite important. It is used to avoid numerical errors and thus problems in the calculation of the shape functions at the nodes (because of the preciosion, they can be out of the convex hull of the discretization */
    double U[3];/*!< \brief Displacements*/
    vector<constitutiveModels *> constitutiveModel;/*!< \brief It is only for the nodal addition */
public:

    classNodes();

    classNodes(int ndim, int me, const vector<double>& xyz);

    inline void setConsModel(vector<constitutiveModels *>& constitutiveModel) {
        this->constitutiveModel = constitutiveModel;
    }; // Common

    inline void setConsModel(constitutiveModels *constitutiveModel) {
        this->constitutiveModel.push_back(constitutiveModel);
    }; // Common

    inline vector<constitutiveModels *>& getConsModel() {
        return constitutiveModel;
    }; // Common

    /*! Destructor
     */
    ~classNodes() {};
};

/*! \brief  This class defines each Dirichlet node in the domain*/
class classDirichlet : public classNodes {

private:
    double prescribedU[3];
    double UU;
    int myNode;
    int direction;
    int meOld;
    int flagDiri; //This value defines ramp (-1) or instantaneous (1). It will be 0 after applying instantaneous Dirichlet BC

public:
    classDirichlet();

    classDirichlet(int ndim, int me, const vector<double>& xyz, double Ux, double Uy, double Uz);

    classDirichlet(int ndim, int me, const vector<double>& xyz, double U, int direction, int numNodes);

    classDirichlet(int ndim, int me, const vector<double>& xyz, double U, int direction, int numNodes, double duration);

    inline double *getPrescribedU() { return prescribedU; }; // Common
    inline int getMyNode() const { return myNode; };// Common
    inline int getDir() const { return direction; };// Common
    inline void setMe(int numNodes) {
        this->meOld = this->me;
        this->me = 0;
        this->me = myNode + direction * numNodes;
    };// Common
    inline int getMeOld() const{ return meOld; }; // Common
    inline double getUU() const{
        return UU;
    }; // Common
    inline int getflagDiri() const{ return flagDiri; };

    inline void resetflagDiri() { this->flagDiri = 0; };// Once the BC has been applied, we resets it
    /*! Destructor
     */
    ~classDirichlet() {};
};

/*! \brief  This class defines each Neumann node in the domain*/
class classNeumann : public classNodes {

private:
    vector<double> _tractionsf;
    string _evolution;
    double _t0, _tf;
    vector<double> _m;


public:
    classNeumann();

    classNeumann(int ndim, int me, const vector<double>& xyz, const vector<double>& Ff, string type, double t0, double tf);

    inline vector<double> getTracAmount(double timeRun, double dt) const{
        vector<double> F(ndim, 0);
        if (timeRun >= _t0 && timeRun <= _tf) {
            for (int k = 0; k < ndim; k++) {
                F[k] = _m[k] * timeRun;//so that the force bc is calculated in the same way as pressure bc
            }
        }
        return F;
    }; // Common
    /*! Destructor
     */
    ~classNeumann() {};
};

/*! \brief  This structure sort is to be used with std::sort function. Property is the direction to sort 0=x, 1=y and 2=z */
struct EntityComp {
    int property;

    EntityComp(int property) { this->property = property; }

    bool operator()(const classNeumann *s1, const classNeumann *s2) const {
        vector<double> xyz1 = s1->getXYZ();
        vector<double> xyz2 = s2->getXYZ();
        if (property == 0) {
            return xyz1[0] < xyz2[0];
        } else if (property == 1) {
            return xyz1[1] < xyz2[1];
        } else if (property == 2) {
            return xyz1[2] < xyz2[2];
        } else {
            cout << "ERROR, you are sorting strange things" << endl;
            exit(-1);
        }
    }
};

#endif
