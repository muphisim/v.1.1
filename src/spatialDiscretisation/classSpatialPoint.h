//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classSpatialPoint_H_
#define _classSpatialPoint_H_

#include "configuration.h"
#include "constitutiveModels.h"
#include "maths.h"

/*! \brief  This class defines the common parts of all spatial points */
class classSpatialPoint {

protected:
    int ndim;
    vector<double> xyz; /*!< \brief Position of the GP*/
    vector<double> xyzCurr;/*!< \brief Position of the node in the current configuration*/
    int me; /*!< \brief Who are you?. Extra information */
    vector<int> neighbours; /*!< \brief This is the array that contains the nodes in which this GP relies*/
    // Carefull with the loops to access to these arrays
    vector<double> phi; /*!< \brief These are the shape functions of all neighboring nodes at this GP. This array has a size of numNeighbours */
    vector<double> dphi; /*!< \brief These are the derivatives of shape functions of all neighboring nodes at this GP. ndim x numNeighbours */
    vector<double> ddphi; /*!< \brief These are the second derivatives shape functions of all neighboring nodes at this GP. ndim x ndim x numNeighbours  */
    double rmax = 0; /*!< \brief This is the radius of the support domain of the node*/
    // Carefull with the loops to access to these arrays
    vector<double> phiMLS; /*!< \brief These are the shape functions of all neighboring nodes at this GP. This array has a size of numNeighbours */
    vector<int> neighboursMLS; /*!< \brief This is the array that contains the nodes in which this GP relies*/
public:
    classSpatialPoint();

    classSpatialPoint(const vector<double>& xyz, int ndim, int me);
    int getDimension() const {return ndim;};
    inline const vector<double>& getXYZ() const { return xyz; };// Common
    void printBasisFunctions(int ndim) const;

    void getF(const vector<double> &uK, int numNodes, int ndim, vector<double> &Fout) const;

    double getNodeToGP(const vector<double> &var, int ndim) const;

    double getAverageNode(const vector<double> &var) const;
    //////////////////////////////////////////////////////////////

    inline const vector<int>& getNeighbours() const { return neighbours; };// Common
    inline const vector<double>& getPhiMLS() const { return phiMLS; };// Common
    inline const vector<double>& getPhi() const { return phi; };// Common
    inline const vector<double>& getDphi() const { return dphi; };// Common
    inline const vector<double>& getDdphi() const { return ddphi; };// Common
    inline void setNeighbours(const vector<int> &tneighbours, int numNeighbours) {
        neighbours.clear();
        neighbours.assign(tneighbours.begin(), tneighbours.end());
    };

    inline void setShapeFunctions(double tphi[], double tdphi[], double tddphi[], int ndim, int numNeighbours) {
        phi.clear();
        dphi.clear();

        phi.assign(tphi, tphi + numNeighbours);
        for (int i = 0; i < numNeighbours; i++) {
            for (int j = 0; j < ndim; j++) {
                dphi.push_back(tdphi[ndim * i + j]); // From fortran notation to c++
            }
        }
    };

    inline void clearPhiMLS() {
        phiMLS.clear();
        neighboursMLS.clear();
    }

    inline void setPhiMLS(double tphi[], int ndim, int numNeighbours) {
        phiMLS.clear();
        phiMLS.assign(tphi, tphi + numNeighbours);
    };

    inline void setNeighboursMLS(const vector<int> &tneighbours, int numNeighbours) {
        neighboursMLS.clear();
        neighboursMLS.assign(tneighbours.begin(), tneighbours.end());
    };

    inline const vector<int>& getNeighboursMLS() const { return neighboursMLS; };// Common


    inline void setXYZCurr(const vector<double>& xyz) {
        this->xyzCurr = xyz;
    };

    inline const vector<double>& getXYZCurr() const { return xyzCurr; }; // Common
    inline double getRmax() const { return rmax; };// Common
    inline void setRmax(double rmax) {
        this->rmax = rmax;
    };// Common
    inline int getMe() const { return me; };// Common

};

#endif
