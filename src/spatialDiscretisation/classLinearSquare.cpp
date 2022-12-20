//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classLinearSquare.cpp
  \brief This file contains the definition of the functions for classLinearSquare
*/

#include "classLinearSquare.h"

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classLinearSquare::getn0(vector<classNodes *>& nod, int ndim) {
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
void classLinearSquare::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
                                         int bulkElement, int surface_order) {
    // I know I am a Square
    int numSs = myNodes.size();
    if (S > numSs) {
        ERROR("You are trying to access to a subElement that does not exist in your original element");
        exit(-1);
    }
    vector<int> mySubNodes;
    if (S != 4) {
        mySubNodes.push_back(myNodes[S - 1]);
        mySubNodes.push_back(myNodes[S]);
    } else {
        mySubNodes.push_back(myNodes[S - 1]);
        mySubNodes.push_back(myNodes[0]);
    }
    classElements *subElement2;
    subElement2 = new classLinearLine(me, "XXX1dLinear", mySubNodes, 1, nod, bulkElement);//XXX call 1D element
    subElements.push_back(subElement2);
}

/*! Constructor for classElements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] nod Array with all nodes in the domain
  @param[in] bulkElement index of mother element
*/
classLinearSquare::classLinearSquare(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
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
double classLinearSquare::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    vector<double> squareNodes;
    squareNodes.resize(ndim * 4, 0);

    const vector<double>& p0 = nodesVec[0];
    const vector<double>& p1 = nodesVec[1];
    const vector<double>& p2 = nodesVec[2];
    const vector<double>& p3 = nodesVec[3];

    for (int i = 0; i < ndim; i++) {
        squareNodes[i * 4] = p0[i];
        squareNodes[i * 4 + 1] = p1[i];
        squareNodes[i * 4 + 2] = p2[i];
        squareNodes[i * 4 + 3] = p3[i];
    }
    double radius = inscribedCircleQuad(squareNodes, ndim);

    return 2 * radius;
};

/*! This function defines just the positions, weights and Jacobians of the GPs inside the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] orderd This value defines the number of GPs inside the element (dramatically important for MM simulations)
  @param[inout] mmGPs Array with all GPs in the MM domain
*/
void classLinearSquare::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                               int orderd, vector<classGPs *> &mmGPs) {
    //int localNdim=this->localNdim;
    int order = 3;
    if (orderd != 0) {
        order = orderd;
    }
    double xw[12][3];
    int nGPs = 0;
    GPsSquare2DValues(order, nGPs, xw);

    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(localNdim *numNodes,
    0);

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
    } else if (ndim == 3) {
        ERROR("MM calls defineGPs for subelements, dont know how you got here");//not necessary to rotate
        exit(-1);
    }

    double area = areaSquare(xyzCoor, ndim);
    this->setArea(area);
    for (int j = 0; j < nGPs; j++) {
        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * localNdim, 0);
        vector<double> Jacobian(localNdim *localNdim,
        0);
        getInterpolationParameters(phi, natDev, xw, j);
        multGeneralMatrices(xyzCoor, localNdim, numNodes, natDev, numNodes, localNdim, Jacobian);
        double J = determinantTensor(Jacobian, localNdim);
        vector<double> xyzIn(ndim, 0);
        //xyzIn[0] = x[0]*(1-xw[j][0]-xw[j][1])+x[1]*xw[j][0]+x[2]*xw[j][1];
        //xyzIn[1] = y[0]*(1-xw[j][0]-xw[j][1])+y[1]*xw[j][0]+y[2]*xw[j][1];
        for (int i = 0; i < ndim; i++) {
            for (int k = 0; k < numNodes; k++) {
                xyzIn[i] += xyz[k][i] * phi[k];
            }
        }
        double weight = 0.5 * xw[j][2];
        classGPs *tempMMGPs = new classGPs(GPs.size(), xyzIn, weight, J, ndim);
        tempMMGPs->setCell(ind, this->getType());
        tempMMGPs->setConsModel(constitutiveModel);
        this->addMyGPs(GPs.size());
        mmGPs.push_back(tempMMGPs);
        GPs.push_back(tempMMGPs);
    }

};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
void classLinearSquare::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim) {

    // int localNdim=this->localNdim;
    double xw[4][3];//3 gp, 2dim+1weight
    int order = 2;
    int nGPs = 0; // Number of GPs per elements. It is not the total number of GPs
    GPsSquare2DValues(order, nGPs, xw); // Integration order, not elements order
    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(localNdim *numNodes,
    0);

    vector<vector<double>> xyz;
    for (int i = 0; i < numNodes; i++) {
        vector<double> xyzTemp = nodes[myNodes[i]]->getXYZ();
        xyz.push_back(xyzTemp);
    }

    if (ndim == 2) {
        for (int j = 0; j < numNodes; j++) {
            //vector<double> xyzT=nodes[myNodes[j]]->getXYZ();
            xyzCoor[j] = xyz[j][0];
            xyzCoor[numNodes + j] = xyz[j][1];
        }
    } else if (ndim == 3)//need to form a local coordinate system R2 on square in R3
    {
        vector<double> p1 = xyz[0]; // Vertices
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

    // Defining the indices of the GPs only for the initialisation
    bool flagInit = false;
    vector<int> listGPs = this->getMyGPs();
    if (listGPs.size() == 0) { flagInit = true; };

    for (int j = 0; j < nGPs; j++) {
        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * localNdim, 0);
        vector<double> Jacobian(localNdim *localNdim,
        0);
        vector<double> invJacobian(localNdim *localNdim,
        0);
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
        double weight = 1 * xw[j][2];
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
        if (flagInit) {
            this->addMyGPs(GPs.size());
        }
        GPs.push_back(tempGP);
    }
};


void classLinearSquare::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                     int bulkElement) {
    int currentIndex = GPs.size();
    defineGPs(ind, GPs, nodes, ndim);
    for (int j = currentIndex; j < GPs.size(); j++) {
        GPs[j]->setBulkElement(bulkElement);
    }
};


void classLinearSquare::getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][3], int j) {
    //shape functions c++ style:
    phi[0] = 0.25 * (1 - xw[j][0]) * (1 - xw[j][1]);
    phi[1] = 0.25 * (1 + xw[j][0]) * (1 - xw[j][1]);
    phi[2] = 0.25 * (1 + xw[j][0]) * (1 + xw[j][1]);
    phi[3] = 0.25 * (1 - xw[j][0]) * (1 + xw[j][1]);


    //natural derivatives c++ style (i*ndim+j convention):
    natDev[0] = -0.25 * (1 - xw[j][1]);
    natDev[2] = 0.25 * (1 - xw[j][1]);
    natDev[4] = 0.25 * (1 + xw[j][1]);
    natDev[6] = -0.25 * (1 + xw[j][1]);
    natDev[1] = -0.25 * (1 - xw[j][0]) ;
    natDev[3] = -0.25 * (1 + xw[j][0]) ;
    natDev[5] = 0.25 * (1 + xw[j][0]) ;
    natDev[7] = 0.25 * (1 - xw[j][0]) ;

};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] tempGP Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] Point Point to be evaluated
*/

void classLinearSquare::defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim,
                                               const vector<double>& Point) {

    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(ndim * numNodes, 0);
    double x[3] = {}, y[3] = {}; // They are all of them Square at this point
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
    x3 = nodes[myNodes[3]]->getXYZ()[0];
    y1 = nodes[myNodes[0]]->getXYZ()[1];
    y2 = nodes[myNodes[1]]->getXYZ()[1];
    y3 = nodes[myNodes[3]]->getXYZ()[1];
//  A[0]=-x1+x2;
//  A[1]=-x1+x3;
//  A[2]=-y1+y3;
//  A[3]=-y1+y2;
//  X[0]=Point[0]-x1;
//  X[1]=Point[1]-y1;
//  inverse(A,ndim,Ainv);
//  multTensorVector(Ainv,ndim,ndim,X,ndim,Xres);
    A[0] = x2 - x1;//=v1x
    A[1] = x3 - x1;//=v2x
    A[2] = y2 - y1;//=v1y
    A[3] = y3 - y1;//=v2y
    //vectors v1 and v2 are nodes 2 and 3 after node 1 is translated to the origin
    //offset "Point" in the same manner
    X[0] = Point[0] - x1;
    X[1] = Point[1] - y1;
    //matrix A=[v1 v2]=[v1x v2x; v1y v2y] maps the isoparametric coordinate system
    //to the actual configuration, more precisely it maps vector [1 0] to v1 and point [0 1] to v2
    //solving A*[1 0]=v1 and A*[0 1]=v2 sets up our mapping matrix
    inverse(A, ndim, Ainv);
    //A*xi=x, so to get xi, invert A
    multTensorVector(Ainv, ndim, ndim, X, ndim, Xres);
    double xw[4][3];
    xw[0][0] = 0.577350269; //gp1 x
    xw[0][1] = 0.577350269; //gp1 y
    xw[0][2] = 1; //gp1 weight
    xw[1][0] = -0.577350269;
    xw[1][1] = 0.577350269;
    xw[1][2] = 1;
    xw[2][0] = 0.577350269;
    xw[2][1] = -0.577350269;
    xw[2][2] = 1;
    xw[3][0] = -0.577350269;
    xw[3][1] = -0.577350269;
    xw[3][2] = 1;//not used
    //xi=Xres[0];
    //eta=Xres[1];
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
