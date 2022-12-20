//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classNodes.cpp
  \brief This file contains the definition of the functions used of the class Node
*/

#include "classNodes.h"

/*! Constructor for classNodes. It should not be called.

*/
classNodes::classNodes() {
    cout << "This constructor should not be called" << endl;
    exit(-1);
};

/*! Constructor for classNodes
  @param[in] ndim Dimension of the domain 
  @param[in] me Label for the node
  @param[in] xyz Position of the node
*/
classNodes::classNodes(int ndim, int me, const vector<double>& xyz) : classSpatialPoint(xyz, ndim, me) {
    this->U[0] = 0;
    this->U[1] = 0;
    this->U[2] = 0;
};

/*! Constructor for classDirichlet
  @param[in] ndim Dimension of the domain 
  @param[in] me Label for the node
  @param[in] xyz Position of the node
  @param[in] Ux Displacement in the x direction
  @param[in] Uy Displacement in the y direction
  @param[in] Uz Displacement in the z direction
*/
classDirichlet::classDirichlet(int ndim, int me, const vector<double>& xyz, double Ux, double Uy, double Uz) : classNodes(ndim,
                                                                                                                   me,
                                                                                                                   xyz) {
    this->prescribedU[0] = Ux;
    this->prescribedU[1] = Uy;
    this->prescribedU[2] = Uz;
};

/*! Constructor for classDirichlet
  @param[in] ndim Dimension of the domain 
  @param[in] me Label for the node
  @param[in] xyz Position of the node
  @param[in] U Displacement
  @param[in] direction Direction for the given displacement
  @param[in] numNodes Number of nodes in the domain
*/
classDirichlet::classDirichlet(int ndim, int me, const vector<double>& xyz, double U, int direction, int numNodes)
        : classNodes(ndim, me, xyz) {
    this->myNode = me;
    this->direction = direction;
    this->me = me + direction * numNodes;
    this->UU = U;
    this->flagDiri = 1;
};

classDirichlet::classDirichlet(int ndim, int me, const vector<double>& xyz, double U, int direction, int numNodes,
                               double duration) : classNodes(ndim, me, xyz) {
    this->myNode = me;
    this->direction = direction;
    this->me = me + direction * numNodes;
    this->UU = U / duration;
    this->flagDiri = -1;
};

/*! Constructor for classDirichlet
  @param[in] ndim Dimension of the domain 
  @param[in] me Label for the node
  @param[in] xyz Position of the node
  @param[in] Ff Force
  @param[in] type Type of Neumann BCs
  @param[in] t0 Init time to apply the pressure
  @param[in] tf Final time to apply the pressure
*/
classNeumann::classNeumann(int ndim, int me, const vector<double>& xyz, const vector<double>& Ff, string type, double t0, double tf)
        : classNodes(ndim, me, xyz) {

    _tractionsf = Ff;
    _evolution = type;
    _t0 = t0;
    _tf = tf;
    _m.resize(ndim, 0);
    for (int k = 0; k < ndim; k++) {
        _m[k] = (_tractionsf[k]) / (_tf - _t0);
    }
};
