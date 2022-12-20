//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classSpatialPoint.cpp
  \brief This file contains all common functions for spatial points
*/

#include "classSpatialPoint.h"

/*! Constructor for classSpatialPoint. It should not be called.

*/
classSpatialPoint::classSpatialPoint() {
    cout << "This constructor should not be called" << endl;
    exit(-1);
};

/*! Constructor for classSpatialPoint
  @param[in] xyz Position of the node
  @param[in] ndim Dimension of the domain
  @param[in] me Label of the point 
*/
classSpatialPoint::classSpatialPoint(const vector<double>& xyz, int ndim, int me) {
    this->me = me;
    this->xyz = xyz;
    this->ndim = ndim;
    this->xyzCurr = xyz;
};

/*! This function evaluates the deformation gradient at the spatial point
  @param[in] uK Displacements of all nodes
  @param[in] numNodes Number of nodes in the domain
  @param[in] ndim Dimension of the domain 
  @param[out] Fout Deformation gradient 
*/
void classSpatialPoint::getF(const vector<double> &uK, int numNodes, int ndim, vector<double> &Fout) const{

    int numNeighbours = neighbours.size();
    int deltak[ndim * ndim];
    memset(deltak, 0, sizeof deltak);
    for (int i = 0; i < ndim; i++) {
        deltak[i * ndim + i] = 1; // defining the kronecker delta
    }
    double dphis[ndim];
    double disp[ndim];
    for (int a = 0; a < numNeighbours; a++) {
        for (int l = 0; l < ndim; l++) {
            dphis[l] = dphi[a * ndim + l];
        }
        int ind = neighbours[a];
        for (int i = 0; i < ndim; i++) {
            disp[i] = uK[ind + i * numNodes];
        }
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
                Fout[i * ndim + j] += dphis[j] * disp[i];
            }
        }
    }
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fout[i * ndim + j] += deltak[i * ndim + j];
        }
    }
};

/*! This function prints the shape functions at this spatial point
  @param[in] ndim Dimension of the domain 
*/
void classSpatialPoint::printBasisFunctions(int ndim) const {

    INFO("Begining with the printing of the information of %i", this->me);
    int tempNumNeighbours = neighbours.size();

    INFO("Number of Neighbours = %i", tempNumNeighbours);
    for (int i = 0; i < tempNumNeighbours; i++) {
        cout << neighbours[i] << endl;
    }
    INFO("Basis functions");
    for (int i = 0; i < tempNumNeighbours; i++) {
        cout << "i=" << i << " " << phi[i] << endl;
        if (!std::isfinite(phi[i])) {
            cout << "There is a nan in the shape functions. Check convex hull and try new rmax" << endl;
            exit(0);
        }
    }
    INFO("Derivatives of the basis functions");
    // dphi
    for (int i = 0; i < tempNumNeighbours; i++) {
        for (int l = 0; l < ndim; l++) {
            cout << dphi[i * ndim + l] << " ";
        }
        cout << endl;
    }
}

/*! This function evaluates a given variable defined at the nodes and provides its value at a Gauss point
  @param[in] var variable of all nodes
  @param[in] numNodes Number of nodes in the domain
  @param[in] ndim Dimension of the domain 
  @return variable at GP 
*/
double classSpatialPoint::getNodeToGP(const vector<double> &var, int ndim) const{
    double varout = 0.0;
    int numNeighbours = neighbours.size();
    for (int a = 0; a < numNeighbours; a++) {
        int ind = neighbours[a];
        varout += var[ind] * phi[a];
    }
    return varout;
};

/*! This function provides the average of a given variable defined at the nodes
  @param[in] var variable of all nodes
  @param[in] numNodes Number of nodes in the domain
  @param[in] ndim Dimension of the domain
  @return variable at GP
*/
double classSpatialPoint::getAverageNode(const vector<double> &var) const{
    double varout = 0.0;
    int numNeighbours = neighbours.size();
    for (int a = 0; a < numNeighbours; a++) {
        int ind = neighbours[a];
        varout += var[ind] * 1. / numNeighbours;
    }
    return varout;
};
