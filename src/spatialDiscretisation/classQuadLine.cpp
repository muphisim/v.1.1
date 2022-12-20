//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classQuadLine.cpp
  \brief This file contains the definition of the functions for classQuadLines
*/

#include "classQuadLine.h"

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classQuadLine::getn0(vector<classNodes *>& nod, int ndim) { //should work assuming midnode has index 2
    vector<double> n0(ndim, 0);
    const vector<double>& xyz1 = nod[myNodes[0]]->getXYZ();
    const vector<double>& xyz2 = nod[myNodes[1]]->getXYZ();

    n0[0] = (xyz2[1] - xyz1[1]);
    n0[1] = (xyz2[0] - xyz1[0]) * -1;

    double no = norm2(n0, ndim);
    for (int i = 0; i < ndim; i++) {
        n0[i] = n0[i] / no;
    }
    return n0;
}

/*! Constructor for classElements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] nod Array with all nodes in the domain
  @param[in] bulkElement index of mother element
*/
classQuadLine::classQuadLine(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
                             int bulkElement) : classElements(me, type, myNodes, ndim, bulkElement) {
    int numNodes = this->myNodes.size();
    vector<vector<double>> nodesVec;
    for (int j = 0; j < numNodes; j++) {
        const vector<double>& currentNode = nod[myNodes[j]]->getXYZ();
        nodesVec.push_back(currentNode);
    }
    this->charLength = this->characteristicLength(ndim, nodesVec);
};

/*! This function calculates the characteristic lenght of the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ndim Dimension of the domain
  @param[in] nodesVec array containing postion vectors of nodes
*/
double classQuadLine::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    const vector<double>& xyz1 = nodesVec[0];
    const vector<double>& xyz2 = nodesVec[1];
    const vector<double>& xyz3 = nodesVec[2];
    double dist1 = calcDist(xyz1, xyz2, ndim);
    double dist2 = calcDist(xyz1, xyz3, ndim);
    double dist3 = calcDist(xyz2, xyz3, ndim);
    double dist = std::min(dist1, dist2);
    dist = std::min(dist, dist3);
    return dist;
};

/*! This function defines just the positions, weights and Jacobians of the GPs inside the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] orderd This value defines the number of GPs inside the element (dramatically important for MM simulations)
  @param[inout] mmGPs Array with all GPs in the MM domain
*/
void
classQuadLine::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int orderd,
                                   vector<classGPs *> &mmGPs) {
    //GPs of subelements are defined using defineGPs function (as in FEM)
    ERROR("How did you even get here?");
    ERROR("MM does not work with quadratic elements. Please use linear");
    exit(-1);

};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
void classQuadLine::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim) {
    int localNdim = this->localNdim;
    double xw[2][2];
    int order = 2;
    int nGPs = 0; // Number of GPs per elements. It is not the total number of GPs
    GPsLine1DValues(order, nGPs, xw); // Integration order, not elements order
    vector<int> neighbours = this->myNodes;
    int numNodes = neighbours.size();
    vector<vector<double>> xyz;
    vector<double> xyzCoorLocal(localNdim * numNodes, 0);
    for (int i = 0; i < numNodes; i++) {
        vector<double> xyzTemp = nodes[myNodes[i]]->getXYZ();
        xyz.push_back(xyzTemp);
    }
    xyzCoorLocal[0] = 0;
    xyzCoorLocal[1] = calcDist(xyz[0], xyz[1], 2);
    xyzCoorLocal[2] = calcDist(xyz[0], xyz[2], 2);//midnode


    for (int j = 0; j < nGPs; j++) {
        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * localNdim, 0);
        vector<double> Jacobian(localNdim * localNdim, 0);
        vector<double> invJacobian(localNdim * localNdim, 0);
        getInterpolationParameters(phi, natDev, xw, j);
        multGeneralMatrices(xyzCoorLocal, localNdim, numNodes, natDev, numNodes, localNdim, Jacobian);
        double J = determinantTensor(Jacobian, localNdim);
        vector<double> xyzIn(ndim, 0);
        for (int i = 0; i < ndim; i++) {
            for (int k = 0; k < numNodes; k++) {
                xyzIn[i] += xyz[k][i] * phi[k];
            }
        }

        double weight = 0.5 * xw[j][1];
        classGPs *tempGP = new classGPs(GPs.size(), xyzIn, weight, J, ndim);
        int numNeighbours = neighbours.size();
        tempGP->setNeighbours(neighbours, numNeighbours);
        tempGP->setCell(ind, this->getType());
        tempGP->setConsModel(constitutiveModel);
        double dphi[localNdim * numNeighbours];   // first derivatives of basis functions
        double ddphi[localNdim * localNdim * numNeighbours];// second derivatives of basis functions

        vector<double> result(numNeighbours * localNdim, 0);
        multGeneralMatrices(natDev, numNeighbours, localNdim, invJacobian, localNdim, localNdim, result);
        for (int i = 0; i < numNeighbours; i++) {
            for (int j = 0; j < localNdim; j++) {
                dphi[i * localNdim + j] = result[i * localNdim + j];
            }
        }

        tempGP->setShapeFunctions(phi, dphi, ddphi, ndim, numNeighbours);
        GPs.push_back(tempGP);
    }
};

void classQuadLine::getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][2], int j) {
    //shape functions c++ style:
    phi[0] = 0.5 * xw[j][0] * xw[j][0] - 0.5 * xw[j][0];
    phi[1] = 0.5 * xw[j][0] * xw[j][0] + 0.5 * xw[j][0];
    phi[2] = -1.0 * xw[j][0] * xw[j][0] + 1.0;


    //natural derivatives c++ style (i*ndim+j convention):
    natDev[0] = 1.0 * xw[j][0] - 0.5;
    natDev[1] = 1.0 * xw[j][0] + 0.5;
    natDev[2] = -2.0 * xw[j][0] + 0.0;


};


void
classQuadLine::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int bulkElement) {
    int currentIndex = GPs.size();
    defineGPs(ind, GPs, nodes, ndim);
    for (int j = currentIndex; j < GPs.size(); j++) {
        GPs[j]->setBulkElement(bulkElement);
    }
};
