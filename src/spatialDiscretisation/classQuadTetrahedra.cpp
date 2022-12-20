//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classQuadTetrahedra.cpp
  \brief This file contains the definition of the functions used of the classQuadTetrahedra
*/

#include "classQuadTetrahedra.h"

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classQuadTetrahedra::getn0(vector<classNodes *>& nod, int ndim) {
    ERROR("Tetrahedra cannot be a subelement for the application of the pressure BCs");
    exit(-1);
}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classQuadTetrahedra::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
                                        int bulkElement, int surface_order) {
    int numSs = myNodes.size();
    if (S > numSs) {
        ERROR("You are trying to access to a subElement that does not exist in your original element");
        exit(-1);
    }
    vector<int> mySubNodes;

    switch (S) {
        case 1: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[5]);
            mySubNodes.push_back(myNodes[6]);
            break;
        }
        case 2: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[8]);
            mySubNodes.push_back(myNodes[4]);
            break;
        }
        case 3: {
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[8]);
            mySubNodes.push_back(myNodes[9]);
            mySubNodes.push_back(myNodes[5]);
            break;
        }
        case 4: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[9]);
            mySubNodes.push_back(myNodes[7]);
            break;
        }
    }
    //for a given S_a:
    //quad tri cornernodes (counterclockwise) a0, a1, a2
    //quad tri midnodes (counterclockwise, first node right adjecent to a0) a3, a4, a5
    //so subelements S_a1 = {a0, a3, a5}, S_a2={a3, a1, a4}, S_a3={a5, a4, a2}, S_a4={a5, a3, a4}
    if (surface_order == 2) {
        classElements *subElement2;
        subElement2 = new classQuadTriangles(me, "2dQuadTriangle", mySubNodes, 2, nod, bulkElement);
        subElements.push_back(subElement2);
    } else if (surface_order == 1) {
        vector<int> S1 = {mySubNodes[0], mySubNodes[3], mySubNodes[5]};
        vector<int> S2 = {mySubNodes[3], mySubNodes[1], mySubNodes[4]};
        vector<int> S3 = {mySubNodes[5], mySubNodes[4], mySubNodes[2]};
        vector<int> S4 = {mySubNodes[5], mySubNodes[3], mySubNodes[4]};
        classElements *sub1 = new classLinearTriangles(me, "2dLinearTriangle", S1, 2, nod, bulkElement);
        classElements *sub2 = new classLinearTriangles(me + 1, "2dLinearTriangle", S2, 2, nod, bulkElement);
        classElements *sub3 = new classLinearTriangles(me + 2, "2dLinearTriangle", S3, 2, nod, bulkElement);
        classElements *sub4 = new classLinearTriangles(me + 3, "2dLinearTriangle", S4, 2, nod, bulkElement);
        subElements.push_back(sub1);
        subElements.push_back(sub2);
        subElements.push_back(sub3);
        subElements.push_back(sub4);
    }
}

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] tempGP Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] Point Point to be evaluated
*/
void classQuadTetrahedra::defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim,
                                              const vector<double>& Point) {
    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(ndim * numNodes, 0);
    double x[10] = {}, y[10] = {}, z[10] = {};
    for (int j = 0; j < numNodes; j++) {
        vector<double> xyz = nodes[myNodes.at(j)]->getXYZ();
        x[j] = xyz[0];
        y[j] = xyz[1];
        z[j] = xyz[2];
    }
    for (int k = 0; k < numNodes; k++) {
        xyzCoor[k] = x[k];
        xyzCoor[numNodes + k] = y[k];
        xyzCoor[2 * numNodes + k] = z[k];
    }

    vector<double> xyzIn(ndim, 0);
    xyzIn[0] = Point[0];
    xyzIn[1] = Point[1];
    xyzIn[2] = Point[2];
    // I have to know chi and eta
    vector<double> A(ndim * ndim, 0);
    vector<double> Ainv(ndim * ndim, 0);
    vector<double> X(ndim, 0);
    vector<double> Xres(ndim, 0);
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    x1 = nodes[myNodes[0]]->getXYZ()[0];
    x2 = nodes[myNodes[1]]->getXYZ()[0];
    x3 = nodes[myNodes[2]]->getXYZ()[0];
    x4 = nodes[myNodes[3]]->getXYZ()[0];
    y1 = nodes[myNodes[0]]->getXYZ()[1];
    y2 = nodes[myNodes[1]]->getXYZ()[1];
    y3 = nodes[myNodes[2]]->getXYZ()[1];
    y4 = nodes[myNodes[3]]->getXYZ()[1];
    z1 = nodes[myNodes[0]]->getXYZ()[2];
    z2 = nodes[myNodes[1]]->getXYZ()[2];
    z3 = nodes[myNodes[2]]->getXYZ()[2];
    z4 = nodes[myNodes[3]]->getXYZ()[2];
    //form the transformation matrix to map isoparametric to reference frame
    A[0] = x2 - x1;
    A[1] = x3 - x1;
    A[2] = x4 - x1;
    A[3] = y2 - y1;
    A[4] = y3 - y1;
    A[5] = y4 - y1;
    A[6] = z2 - z1;
    A[7] = z3 - z1;
    A[8] = z4 - z1;

    X[0] = Point[0] - x1;
    X[1] = Point[1] - y1;
    X[2] = Point[2] - z1;
    inverse(A, ndim, Ainv);
    multTensorVector(Ainv, ndim, ndim, X, ndim, Xres);
    double xw[1][4];
    xw[0][0] = Xres[0];
    xw[0][1] = Xres[1];
    xw[0][2] = Xres[2];
    xw[0][2] = 0;//not used

    double phi[numNodes];     // basis functions
    vector<double> natDev(numNodes * ndim, 0);
    vector<double> Jacobian(ndim * ndim, 0);
    vector<double> invJacobian(ndim * ndim, 0);
    getInterpolationParameters(phi, natDev, xw, 0);
    multGeneralMatrices(xyzCoor, ndim, numNodes, natDev, numNodes, ndim, Jacobian);
    double J = determinantTensor(Jacobian, ndim);
    inverse(Jacobian, ndim, invJacobian);

    double weight = 1. / 24.;
    tempGP = new classGPs(ind, xyzIn, weight, J, ndim);
    int numNeighbours = neighbours.size();
    tempGP->setNeighbours(neighbours, numNeighbours);
    tempGP->setCell(ind, this->getType());
    tempGP->setConsModel(constitutiveModel);
    double dphi[ndim * numNeighbours];   // first derivatives of basis functions
    double ddphi[ndim * ndim * numNeighbours];// second derivatives of basis functions
    inverse(Jacobian, ndim, invJacobian);
    vector<double> result(numNeighbours * ndim, 0);
    multGeneralMatrices(natDev, numNeighbours, ndim, invJacobian, ndim, ndim, result);
    for (int i = 0; i < numNeighbours; i++) {
        for (int j = 0; j < ndim; j++) {
            dphi[i * ndim + j] = result[i * ndim + j];
        }
    }
    tempGP->setShapeFunctions(phi, dphi, ddphi, ndim, numNeighbours);
};

/*! Constructor for classElements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] nod Array with all nodes in the domain
  @param[in] bulkElement index of mother element
*/
classQuadTetrahedra::classQuadTetrahedra(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
                                         int bulkElement) : classElements(me, type, myNodes, ndim, bulkElement) {
    int numNodes = this->myNodes.size();
    vector<vector<double> > nodesVec;
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
double classQuadTetrahedra::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    double minDist = 1E25;
    //get all nodes in an array

    //divide tetrahedron into subelements, calculate radius of the inscribed sphere for each subelement and pick the smallest
    vector<vector<int> > subelements = {{6, 5, 2, 9},
                                        {4, 8, 1, 5},
                                        {0, 4, 6, 7},
                                        {3, 5, 6, 9},
                                        {8, 3, 4, 5},
                                        {7, 3, 4, 6},
                                        {0, 4, 6, 7}};
    double radii[7];
    for (int j = 0; j < 7; j++)//tetrahedron subdivides into 7 subtetrahedras with cornernodes and midnodes
    {
        //each subelement can again be subdivided into 4 smaller subsubelements (all subelements are tetrahedras)
        //the subdivision of each tetrahedral element is such that each of the subsubelements has a base corresponding to one of
        //the triangular surfaces of the subelement, and the last node is center of the inscribed sphere. The sum of volumes of
        //subsubelements is equal to the volume of the subelement, and since v=1/3*a*h, where h of a subsubelement is the radius of the
        //sphere, we know that total volume of the subelement = 1/3*(s1+s2+s3+s4)*r
        const vector<double>& p0 = nodesVec[subelements[j][0]];
        const vector<double>& p1 = nodesVec[subelements[j][1]];
        const vector<double>& p2 = nodesVec[subelements[j][2]];
        const vector<double>& p3 = nodesVec[subelements[j][3]];//vertecies of current subelement in the loop, for convenience
        vector<double> tetraNodes;
        tetraNodes.resize(ndim * 4);
        for (int i = 0; i < ndim; i++) {
            tetraNodes[i * 4] = p0[i];
            tetraNodes[i * 4 + 1] = p1[i];
            tetraNodes[i * 4 + 2] = p2[i];
            tetraNodes[i * 4 + 3] = p3[i];

        }
        radii[j] = inscribedSphere(tetraNodes, ndim);
    }
    double rmin = radii[0];
    for (int i = 1; i < 7; i++) {
        if (rmin > radii[i])
            rmin = radii[i];
    }
    if (rmin < minDist) {
        minDist = 2 * rmin;//diameter!!! radius was way slow
    }
    return minDist;
};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
void classQuadTetrahedra::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim) {

    int nGPs = 0;
    int order = 2;
    if (order == 1) {
        nGPs = 1;
    } else if (order == 2) {
        nGPs = 4;
    } else {
        ERROR("The integration order %i for tetrahedra is not available. Orders 1 and 2 are the ones available in MuPhiSim",
              order);
        exit(-1);
    }
    double xw[4][4];//1gp, 3dim+w
    GPsTetrahedron3DValues(order, nGPs, xw);

    vector<int> neighbours; // Temporal list of nodes in which a given GPs lie
    int numNodes = this->myNodes.size();
    neighbours = myNodes;
    vector<double> xyzCoor(ndim * numNodes, 0);


    vector<vector<double>> xyz;
    for (int i = 0; i < numNodes; i++) {
        vector<double> xyzTemp = nodes[myNodes[i]]->getXYZ();
        xyz.push_back(xyzTemp);
        xyzCoor[i] = xyzTemp[0];
        xyzCoor[numNodes + i] = xyzTemp[1];
        xyzCoor[2 * numNodes + i] = xyzTemp[2];
    }

    for (int j = 0; j < nGPs; j++) {
        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * ndim, 0);
        vector<double> Jacobian(ndim * ndim, 0);
        vector<double> invJacobian(ndim * ndim, 0);
        getInterpolationParameters(phi, natDev, xw, j);
        multGeneralMatrices(xyzCoor, ndim, numNodes, natDev, numNodes, ndim, Jacobian);
        double J = determinantTensor(Jacobian, ndim);
        inverse(Jacobian, ndim, invJacobian);


        vector<double> xyzIn(ndim, 0);
        for (int i = 0; i < ndim; i++) {
            for (int k = 0; k < numNodes; k++) {
                xyzIn[i] += xyz[k][i] * phi[k];
            }
        }

        double weight = xw[j][3];// It is already divided by 1/6
        classGPs *tempGP = new classGPs(GPs.size(), xyzIn, weight, J, ndim);
        int numNeighbours = neighbours.size();
        tempGP->setNeighbours(neighbours, numNeighbours);

        double dphi[ndim * numNeighbours];   // first derivatives of basis functions
        double ddphi[ndim * ndim * numNeighbours];// second derivatives of basis functions
        vector<double> result(numNeighbours * ndim, 0);
        multGeneralMatrices(natDev, numNeighbours, ndim, invJacobian, ndim, ndim, result);
        for (int i = 0; i < numNeighbours; i++) {
            for (int j = 0; j < ndim; j++) {
                dphi[i * ndim + j] = result[i * ndim + j];
            }
        }
        tempGP->setShapeFunctions(phi, dphi, ddphi, ndim, numNeighbours);
        tempGP->setCell(ind, this->getType());
        tempGP->setConsModel(constitutiveModel);
        GPs.push_back(tempGP);
    }
};

void classQuadTetrahedra::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                    int bulkElement) {
    INFO("should not be used");
};

void classQuadTetrahedra::getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][4], int j) {
    //shape functions c++ style:
    phi[0] = 2.0 * xw[j][2] * xw[j][2] - 1.0 * xw[j][2];
    phi[1] = 2.0 * xw[j][1] * xw[j][1] - 1.0 * xw[j][1];
    phi[2] = 2.0 * xw[j][0] * xw[j][0] - 1.0 * xw[j][0];
    phi[3] = 2.0 * xw[j][0] * xw[j][0] + 4.0 * xw[j][0] * xw[j][1] + 4.0 * xw[j][0] * xw[j][2] - 3.0 * xw[j][0] +
             2.0 * xw[j][1] * xw[j][1] + 4.0 * xw[j][1] * xw[j][2] - 3.0 * xw[j][1] + 2.0 * xw[j][2] * xw[j][2] -
             3.0 * xw[j][2] + 1.0;
    phi[4] = 4.0 * xw[j][1] * xw[j][2];
    phi[5] = 4.0 * xw[j][0] * xw[j][1];
    phi[6] = 4.0 * xw[j][0] * xw[j][2];
    phi[7] = -4.0 * xw[j][0] * xw[j][2] - 4.0 * xw[j][1] * xw[j][2] - 4.0 * xw[j][2] * xw[j][2] + 4.0 * xw[j][2];
    phi[8] = -4.0 * xw[j][0] * xw[j][1] - 4.0 * xw[j][1] * xw[j][1] - 4.0 * xw[j][1] * xw[j][2] + 4.0 * xw[j][1];
    phi[9] = -4.0 * xw[j][0] * xw[j][0] - 4.0 * xw[j][0] * xw[j][1] - 4.0 * xw[j][0] * xw[j][2] + 4.0 * xw[j][0];


    //natural derivatives c++ style (i*ndim+j convention):
    natDev[0] = 0.0 - 0.0;
    natDev[3] = 0.0 - 0.0;
    natDev[6] = 4.0 * xw[j][0] - 1.0;
    natDev[9] = 4.0 * xw[j][0] + 4.0 * xw[j][1] + 4.0 * xw[j][2] - 3.0 + 0.0 + 0.0 - 0.0 + 0.0 - 0.0 + 0.0;
    natDev[12] = 0.0;
    natDev[15] = 4.0 * xw[j][1];
    natDev[18] = 4.0 * xw[j][2];
    natDev[21] = -4.0 * xw[j][2] - 0.0 - 0.0 + 0.0;
    natDev[24] = -4.0 * xw[j][1] - 0.0 - 0.0 + 0.0;
    natDev[27] = -8.0 * xw[j][0] - 4.0 * xw[j][1] - 4.0 * xw[j][2] + 4.0;
    natDev[1] = 0.0 - 0.0;
    natDev[4] = 4.0 * xw[j][1] - 1.0;
    natDev[7] = 0.0 - 0.0;
    natDev[10] = 0.0 + 4.0 * xw[j][0] + 0.0 - 0.0 + 4.0 * xw[j][1] + 4.0 * xw[j][2] - 3.0 + 0.0 - 0.0 + 0.0;
    natDev[13] = 4.0 * xw[j][2];
    natDev[16] = 4.0 * xw[j][0];
    natDev[19] = 0.0;
    natDev[22] = -0.0 - 4.0 * xw[j][2] - 0.0 + 0.0;
    natDev[25] = -4.0 * xw[j][0] - 8.0 * xw[j][1] - 4.0 * xw[j][2] + 4.0;
    natDev[28] = -0.0 - 4.0 * xw[j][0] - 0.0 + 0.0;
    natDev[2] = 4.0 * xw[j][2] - 1.0;
    natDev[5] = 0.0 - 0.0;
    natDev[8] = 0.0 - 0.0;
    natDev[11] = 0.0 + 0.0 + 4.0 * xw[j][0] - 0.0 + 0.0 + 4.0 * xw[j][1] - 0.0 + 4.0 * xw[j][2] - 3.0 + 0.0;
    natDev[14] = 4.0 * xw[j][1];
    natDev[17] = 0.0;
    natDev[20] = 4.0 * xw[j][0];
    natDev[23] = -4.0 * xw[j][0] - 4.0 * xw[j][1] - 8.0 * xw[j][2] + 4.0;
    natDev[26] = -0.0 - 0.0 - 4.0 * xw[j][1] + 0.0;
    natDev[29] = -0.0 - 0.0 - 4.0 * xw[j][0] + 0.0;


};

/*! This function defines just the positions, weights and Jacobians of the GPs inside the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] orderd This value defines the number of GPs inside the element (dramatically important for MM simulations)
  @param[inout] mmGPs Array with all GPs in the MM domain
*/
void classQuadTetrahedra::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                              int orderd, vector<classGPs *> &mmGPs) {
    ERROR("MM should not be used with quadratic elements. Please use linear.");
    exit(-1);

}
