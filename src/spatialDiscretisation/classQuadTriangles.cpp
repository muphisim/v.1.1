//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classQuadTriangles.cpp
  \brief This file contains the definition of the functions for classQuadTriangles
*/

#include "classQuadTriangles.h"

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classQuadTriangles::getn0(vector<classNodes *>& nod, int ndim) {
    const vector<double>& xyz1 = nod[myNodes[0]]->getXYZ();
    const vector<double>& xyz2 = nod[myNodes[1]]->getXYZ();
    const vector<double>& xyz3 = nod[myNodes[2]]->getXYZ();

    vector<double> n0(ndim, 0);

    for (int i = 0; i < ndim; i++) {
        int coord1, coord2;
        if (i == 0) //nx
        {
            coord1 = 1;
            coord2 = 2;
        } else if (i == 1) //ny
        {
            coord1 = 2;
            coord2 = 0;
        } else {
            coord1 = 0;
            coord2 = 1;
        }
        n0[i] = (xyz2[coord1] - xyz1[coord1]) * (xyz3[coord2] - xyz1[coord2]) -
                (xyz3[coord1] - xyz1[coord1]) * (xyz2[coord2] - xyz1[coord2]);
    }
    double no = norm2(n0, ndim);
    for (int i = 0; i < ndim; i++) {
        n0[i] = n0[i] / no;
    }
    return n0;

}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classQuadTriangles::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
                                       int bulkElement, int surface_order) {
    // I know I am a triangle
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
            mySubNodes.push_back(myNodes[3]);
            break;
        }
        case 2: {
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[4]);
            break;
        }
        case 3: {
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[5]);
            break;
        }
    }

    if (surface_order == 2) {
        classElements *subElement2;
        subElement2 = new classQuadLine(me, "XXX1dQuad", mySubNodes, 1, nod,
                                        bulkElement);//XXX call 1D element, changed from Linear to Quad
        subElements.push_back(subElement2);
    } else if (surface_order == 1) {
        vector<int> S1 = {mySubNodes[0], mySubNodes[2]};
        vector<int> S2 = {mySubNodes[2], mySubNodes[1]};
        classElements *sub1 = new classLinearLine(me, "XXX1dLinear", S1, 1, nod, bulkElement);
        classElements *sub2 = new classLinearLine(me + 1, "XXX1dLinear", S2, 1, nod, bulkElement);
        subElements.push_back(sub1);
        subElements.push_back(sub2);
    }

}


/*! Constructor for classElements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] nod Array with all nodes in the domain
  @param[in] bulkElement index of mother element
*/
classQuadTriangles::classQuadTriangles(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
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
double classQuadTriangles::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    double minDist = 1E25;
    // get radii of circles inscribed by triangles formed by nodes (0,3,5),(3,1,4),(4,2,5),(3,4,5)
    //get coordinates of all nodes
    double radii[4];
    vector<vector<int>> subelements = {{0, 3, 5},
                                       {3, 1, 4},
                                       {4, 2, 5},
                                       {3, 4, 5}};
    for (int j = 0; j < 4; j++)//subdivides into 4 subtriangles
    {
        vector<double> triangleNodes;
        triangleNodes.resize(ndim * 3, 0);

        const vector<double>& p0 = nodesVec[subelements[j][0]];
        const vector<double>& p1 = nodesVec[subelements[j][1]];
        const vector<double>& p2 = nodesVec[subelements[j][2]];

        for (int i = 0; i < ndim; i++) {
            triangleNodes[i * 3] = p0[i];
            triangleNodes[i * 3 + 1] = p1[i];
            triangleNodes[i * 3 + 2] = p2[i];
        }
        radii[j] = inscribedCircle(triangleNodes, ndim);

    }
    minDist = radii[0];
    for (int i = 1; i < 4; i++) {
        if (radii[i] < minDist) {
            minDist = radii[i];
        }
    }


    return 2 * minDist;
};

/*! This function defines just the positions, weights and Jacobians of the GPs inside the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] orderd This value defines the number of GPs inside the element (dramatically important for MM simulations)
  @param[inout] mmGPs Array with all GPs in the MM domain
*/
void classQuadTriangles::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                             int orderd, vector<classGPs *> &mmGPs) {

    ERROR("Quadratic elements are not supported in MM. Please use linear");
    exit(-1);
};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of the element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
void classQuadTriangles::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim) {
    int localNdim = this->localNdim;
    double xw[3][3];
    int order = 2;
    int nGPs = 0; // Number of GPs per elements. It is not the total number of GPs
    GPsTriangle2DValues(order, nGPs, xw); // Integration order, not elements order
    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(localNdim * numNodes, 0);

    vector<vector<double>> xyz;
    for (int i = 0; i < numNodes; i++) {
        vector<double> xyzTemp = nodes[myNodes[i]]->getXYZ();
        xyz.push_back(xyzTemp);
    }

    if (ndim == 2) {
        for (int j = 0; j < numNodes; j++) {
            xyzCoor[j] = xyz[j][0];
            xyzCoor[numNodes + j] = xyz[j][1];
        }
    } else if (ndim == 3)//need to form a local coordinate system R2 on triangle in R3
    {
        vector<double> p1 = xyz[0];
        vector<double> p2 = xyz[1];
        vector<double> p3 = xyz[2];

        vector<double> a(ndim, 0);
        vector<double> b(ndim, 0);
        vectorArithmetic(p2, p1, a, ndim, false);
        vectorArithmetic(p3, p1, b, ndim, false);//b=p3-p1
        vector<double> map;
        createOrthogonalBasisOnPlane(a, b, ndim, map);
        //all points to be mapped must also be offset by p1
        vector<double> xyzCoorTemp(ndim * numNodes, 0);
        for (int k = 0; k < numNodes; k++) {
            vector<double> xyzOff(ndim, 0);
            vectorArithmetic(xyz[k], p1, xyzOff, ndim, false);
            xyzCoorTemp[k] = xyzOff[0];
            xyzCoorTemp[numNodes + k] = xyzOff[1];
            xyzCoorTemp[2 * numNodes + k] = xyzOff[2];
        }
        multGeneralMatrices(map, localNdim, ndim, xyzCoorTemp, ndim, numNodes, xyzCoor);
    }

    for (int j = 0; j < nGPs; j++) {

        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * localNdim, 0);
        vector<double> Jacobian(localNdim * localNdim, 0);
        vector<double> invJacobian(localNdim * localNdim, 0);
        getInterpolationParameters(phi, natDev, xw, j);
        multGeneralMatrices(xyzCoor, localNdim, numNodes, natDev, numNodes, localNdim, Jacobian);
        double J = determinantTensor(Jacobian, localNdim);
        inverse(Jacobian, localNdim, invJacobian);

        vector<double> xyzIn(ndim, 0);
        for (int i = 0; i < ndim; i++) {
            for (int k = 0; k < numNodes; k++) {
                xyzIn[i] += xyz[k][i] * phi[k];
            }
        }
        //xyzIn[0] = x[0]*(1-xw[j][0]-xw[j][1])+x[1]*xw[j][0]+x[2]*xw[j][1];
        //xyzIn[1] = y[0]*(1-xw[j][0]-xw[j][1])+y[1]*xw[j][0]+y[2]*xw[j][1];
        double weight = 0.5 * xw[j][2];
        classGPs *tempGP = new classGPs(GPs.size(), xyzIn, weight, J, ndim);
        int numNeighbours = neighbours.size();
        tempGP->setNeighbours(neighbours, numNeighbours);
        tempGP->setCell(ind, this->getType());
        tempGP->setConsModel(constitutiveModel);

        double dphi[localNdim * numNeighbours];   // first derivatives of basis functions
        double ddphi[localNdim * localNdim * numNeighbours];// second derivatives of basis functions;
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

void classQuadTriangles::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                   int bulkElement) {
    int currentIndex = GPs.size();
    defineGPs(ind, GPs, nodes, ndim);
    for (int j = currentIndex; j < GPs.size(); j++) {
        GPs[j]->setBulkElement(bulkElement);
    }
};

void classQuadTriangles::getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][3], int j) {

    //shape functions c++ style:
    phi[0] = 2.0 * xw[j][0] * xw[j][0] + 4.0 * xw[j][0] * xw[j][1] - 3.0 * xw[j][0] + 2.0 * xw[j][1] * xw[j][1] -
             3.0 * xw[j][1] + 1.0;
    phi[1] = 2.0 * xw[j][0] * xw[j][0] - 1.0 * xw[j][0];
    phi[2] = 2.0 * xw[j][1] * xw[j][1] - 1.0 * xw[j][1];
    phi[3] = -4.0 * xw[j][0] * xw[j][0] - 4.0 * xw[j][0] * xw[j][1] + 4.0 * xw[j][0];
    phi[4] = 4.0 * xw[j][0] * xw[j][1];
    phi[5] = -4.0 * xw[j][0] * xw[j][1] - 4.0 * xw[j][1] * xw[j][1] + 4.0 * xw[j][1];


//natural derivatives c++ style (i*ndim+j convention):
    natDev[0] = 4.0 * xw[j][0] + 4.0 * xw[j][1] - 3.0 + 0.0 - 0.0 + 0.0;
    natDev[2] = 4.0 * xw[j][0] - 1.0;
    natDev[4] = 0.0 - 0.0;
    natDev[6] = -8.0 * xw[j][0] - 4.0 * xw[j][1] + 4.0;
    natDev[8] = 4.0 * xw[j][1];
    natDev[10] = -4.0 * xw[j][1] - 0.0 + 0.0;
    natDev[1] = 0.0 + 4.0 * xw[j][0] - 0.0 + 4.0 * xw[j][1] - 3.0 + 0.0;
    natDev[3] = 0.0 - 0.0;
    natDev[5] = 4.0 * xw[j][1] - 1.0;
    natDev[7] = -0.0 - 4.0 * xw[j][0] + 0.0;
    natDev[9] = 4.0 * xw[j][0];
    natDev[11] = -4.0 * xw[j][0] - 8.0 * xw[j][1] + 4.0;

}

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] tempGP Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] Point Point to be evaluated
*/

void classQuadTriangles::defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim,
                                             const vector<double>& Point) {

    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(ndim * numNodes, 0);
    double x[6] = {}, y[6] = {}; // They are all of them triangles at this point
    for (int j = 0; j < numNodes; j++) {
        vector<double> xyz = nodes[myNodes[j]]->getXYZ();
        x[j] = xyz[0];
        y[j] = xyz[1];
    }
    for (int k = 0; k < numNodes; k++) {
        xyzCoor[k] = x[k];
        xyzCoor[numNodes + k] = y[k];
    }

    vector<double> xyzIn(ndim, 0);
    xyzIn[0] = Point[0];
    xyzIn[1] = Point[1];
    // I have to know chi and eta
    vector<double> A(ndim * ndim, 0);
    vector<double> Ainv(ndim * ndim, 0);
    vector<double> X(ndim, 0);
    vector<double> Xres(ndim, 0);
    double x1, x2, x3, y1, y2, y3;
    x1 = nodes[myNodes[0]]->getXYZ()[0];
    x2 = nodes[myNodes[1]]->getXYZ()[0];
    x3 = nodes[myNodes[2]]->getXYZ()[0];
    y1 = nodes[myNodes[0]]->getXYZ()[1];
    y2 = nodes[myNodes[1]]->getXYZ()[1];
    y3 = nodes[myNodes[2]]->getXYZ()[1];
    //form the transformation matrix to map isoparametric to reference frame
    A[0] = x2 - x1;
    A[1] = x3 - x1;
    A[2] = y2 - y1;
    A[3] = y3 - y1;
    X[0] = Point[0] - x1;
    X[1] = Point[1] - y1;
    inverse(A, ndim, Ainv);
    multTensorVector(Ainv, ndim, ndim, X, ndim, Xres);
    double xw[1][3];
    xw[0][0] = Xres[0];
    xw[0][1] = Xres[1];
    xw[0][2] = 0;//not used

    double phi[numNodes];     // basis functions
    vector<double> natDev(numNodes * ndim, 0);
    vector<double> Jacobian(ndim * ndim, 0);
    vector<double> invJacobian(ndim * ndim, 0);
    getInterpolationParameters(phi, natDev, xw, 0);
    multGeneralMatrices(xyzCoor, ndim, numNodes, natDev, numNodes, ndim, Jacobian);
    double J = determinantTensor(Jacobian, ndim);
    inverse(Jacobian, ndim, invJacobian);


    double weight = 0;
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
