//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classLinearHexahedra.cpp
  \brief This file contains the definition of the functions used of the classLinearHexahedra
*/

#include "classLinearHexahedra.h"

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classLinearHexahedra::getn0(vector<classNodes *>& nod, int ndim){
    ERROR("Hexahedra cannot be a subelement for the application of the pressure BCs");
    exit(-1);
}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classLinearHexahedra::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
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
            mySubNodes.push_back(myNodes[3]);
            break;
        }
        case 2: {
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[5]);
            break;
        }
        case 3: {
            mySubNodes.push_back(myNodes[0]);
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[5]);
            mySubNodes.push_back(myNodes[1]);
            break;
        }
        case 4: {
            mySubNodes.push_back(myNodes[1]);
            mySubNodes.push_back(myNodes[5]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[2]);
            break;
        }
        case 5: {
            mySubNodes.push_back(myNodes[2]);
            mySubNodes.push_back(myNodes[6]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[3]);
            break;
        }
        case 6: {
            mySubNodes.push_back(myNodes[3]);
            mySubNodes.push_back(myNodes[7]);
            mySubNodes.push_back(myNodes[4]);
            mySubNodes.push_back(myNodes[0]);
            break;
        }
    }

    //cout<< mySubNodes[0]<< " d "<< mySubNodes[1]<<endl;
    classElements *subElement2;
    subElement2 = new classLinearSquare(me, "2dLinearSquare", mySubNodes, 2, nod, bulkElement);
    subElements.push_back(subElement2);

}

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] tempGP Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] Point Point to be evaluated
*/
void classLinearHexahedra::defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim,
                                                const vector<double>& Point) {
    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(ndim * numNodes, 0);
    double x[8] = {}, y[8] = {}, z[8] = {}; // They are all of them triangles at this point
    for (int j = 0; j < numNodes; j++) {
        const vector<double>& xyz = nodes[myNodes.at(j)]->getXYZ();
        x[j] = xyz[0];
        y[j] = xyz[1];
        z[j] = xyz[2];
    }
    for (int k = 0; k < numNodes; k++) {
        xyzCoor[k] = x[k];
        xyzCoor[numNodes + k] = y[k];
        xyzCoor[2 * numNodes + k] = z[k];
    }
    vector<double> A(ndim * ndim, 0);
    vector<double> Ainv(ndim * ndim, 0);
    vector<double> X(ndim, 0);
    vector<double> Xres(ndim, 0);
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    x1 = nodes[myNodes[0]]->getXYZ()[0];
    x2 = nodes[myNodes[1]]->getXYZ()[0];
    x3 = nodes[myNodes[3]]->getXYZ()[0];
    x4 = nodes[myNodes[4]]->getXYZ()[0];
    y1 = nodes[myNodes[0]]->getXYZ()[1];
    y2 = nodes[myNodes[1]]->getXYZ()[1];
    y3 = nodes[myNodes[3]]->getXYZ()[1];
    y4 = nodes[myNodes[4]]->getXYZ()[1];
    z1 = nodes[myNodes[0]]->getXYZ()[2];
    z2 = nodes[myNodes[1]]->getXYZ()[2];
    z3 = nodes[myNodes[3]]->getXYZ()[2];
    z4 = nodes[myNodes[4]]->getXYZ()[2];
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
    xw[0][0] = 2*Xres[0] - 1;
    xw[0][1] = 2*Xres[1] - 1;
    xw[0][2] = 2*Xres[2] - 1;
    xw[0][3] = 1;//not used
    double phi[numNodes];     // basis functions
    vector<double> natDev(numNodes * ndim, 0);
    vector<double> Jacobian(ndim * ndim, 0);
    vector<double> invJacobian(ndim * ndim, 0);
    getInterpolationParameters(phi, natDev, xw, 0);
    multGeneralMatrices(xyzCoor, ndim, numNodes, natDev, numNodes, ndim, Jacobian);
    double J = determinantTensor(Jacobian, ndim);
    inverse(Jacobian, ndim, invJacobian);

    vector<double> xyzIn(ndim, 0);
    xyzIn[0] = Point[0];
    xyzIn[1] = Point[1];
    xyzIn[2] = Point[2];
    double weight = 8;
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
classLinearHexahedra::classLinearHexahedra(int me, string type,  const vector<int>& myNodes, int ndim,
                                             vector<classNodes *>& nod, int bulkElement) : classElements(me, type,
                                                                                                        myNodes, ndim,
                                                                                                        bulkElement) {
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
double classLinearHexahedra::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    const vector<double>& p0 = nodesVec[0];
    const vector<double>& p1 = nodesVec[1];
    const vector<double>& p2 = nodesVec[2];
    const vector<double>& p3 = nodesVec[3];
    const vector<double>& p4 = nodesVec[4];
    const vector<double>& p5 = nodesVec[5];
    const vector<double>& p6 = nodesVec[6];
    const vector<double>& p7 = nodesVec[7];
    vector<double> hexaNodes;
    hexaNodes.resize(ndim * 8);
    for (int i = 0; i < ndim; i++) {
        hexaNodes[i * 8] = p0[i];
        hexaNodes[i * 8 + 1] = p1[i];
        hexaNodes[i * 8 + 2] = p2[i];
        hexaNodes[i * 8 + 3] = p3[i];
        hexaNodes[i * 8 + 4] = p4[i];
        hexaNodes[i * 8 + 5] = p5[i];
        hexaNodes[i * 8 + 6] = p6[i];
        hexaNodes[i * 8 + 7] = p7[i];

    }
    double radius = inscribedSphereHex(hexaNodes, ndim);
    return 2 * radius;
};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
void classLinearHexahedra::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim) {

    int nGPs = 0;
    int order = 2;
    if (order == 1) {
        nGPs = 1;
    } else if (order == 2) {
        nGPs = 8;
    } else {
        ERROR("The integration order %i for Hexahedra is not available. Orders 1 and 2 are the ones available in MuPhiSim",
              order);
        exit(-1);
    }
    //double xi[nGPs],eta[nGPs],zeta[nGPs], weights[nGPs];
    double xw[8][4];//1gp, 3dim+w
    GPsHexahedron3DValues(order, nGPs, xw);//xi, eta, zeta, weights);

    vector<int> neighbours; // Temporal list of nodes in which a given GPs lie
    int numNodes = this->myNodes.size();
    neighbours = myNodes;
    vector<double> xyzCoor(ndim * numNodes, 0);
    //double x[4]={},y[4]={}, z[4]={}; // They are all of them triangles at this point

    //for(int j=0;j<numNodes;j++){
    //  xyz=nodes[myNodes[j]]->getXYZ();
    //  x[j]=xyz[0];
    //  y[j]=xyz[1];
    //  z[j]=xyz[2];
    //}
    vector<vector<double>> xyz;
    for (int i = 0; i < numNodes; i++) {
        vector<double> xyzTemp = nodes[myNodes[i]]->getXYZ();
        xyz.push_back(xyzTemp);
        xyzCoor[i] = xyzTemp[0];
        xyzCoor[numNodes + i] = xyzTemp[1];
        xyzCoor[2 * numNodes + i] = xyzTemp[2];
    }


    // https://people.fh-landshut.de/~maurer/numeth/node148.html
    //for(int k=0;k<numNodes;k++){
    //  xyzCoor[k]=x[k];
    //  xyzCoor[numNodes+k]=y[k];
    //  xyzCoor[2*numNodes+k]=z[k];
    //}

    //vector<double> Jacobian(ndim*ndim,0);

    //vector<double> invJacobian(ndim*ndim,0);
    //vector <double> natDev;
    //natDev.resize(myNodes.size()*ndim,0);

    // https://people.fh-landshut.de/~maurer/numeth/node148.html
	
    // Defining the indices of the GPs only for the initialisation
    bool flagInit = false;
    vector<int> listGPs = this->getMyGPs();
    if (listGPs.size() == 0) { flagInit = true; }; 

    for (int j = 0; j < nGPs; j++) {
        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * ndim, 0);
        vector<double> Jacobian(ndim * ndim, 0);
        vector<double> invJacobian(ndim * ndim, 0);
        getInterpolationParameters(phi, natDev, xw, j);
        multGeneralMatrices(xyzCoor, ndim, numNodes, natDev, numNodes, ndim, Jacobian);
        double J = determinantTensor(Jacobian, ndim);
        inverse(Jacobian, ndim, invJacobian);


        vector<double> xyzIn(localNdim, 0);
        for (int i = 0; i < ndim; i++) {
            for (int k = 0; k < numNodes; k++) {
                xyzIn[i] += xyz[k][i] * phi[k];
            }
        }

        double weight = 1*xw[j][3];// It is already divided by 1/6
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
        if (flagInit) {
            this->addMyGPs(GPs.size());
        } //This vector is a indices list from 0 to maximum number of GPs
        
        GPs.push_back(tempGP);
    }
};

void classLinearHexahedra::getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][4], int j) {
    //shape functions c++ style:
    phi[0] = 0.125 * (1 - xw[j][0]) * (1 - xw[j][1]) * (1 - xw[j][2]);
    phi[1] = 0.125 * (1 + xw[j][0]) * (1 - xw[j][1]) * (1 - xw[j][2]);
    phi[2] = 0.125 * (1 + xw[j][0]) * (1 + xw[j][1]) * (1 - xw[j][2]);
    phi[3] = 0.125 * (1 - xw[j][0]) * (1 + xw[j][1]) * (1 - xw[j][2]);
    phi[4] = 0.125 * (1 - xw[j][0]) * (1 - xw[j][1]) * (1 + xw[j][2]);
    phi[5] = 0.125 * (1 + xw[j][0]) * (1 - xw[j][1]) * (1 + xw[j][2]);
    phi[6] = 0.125 * (1 + xw[j][0]) * (1 + xw[j][1]) * (1 + xw[j][2]);
    phi[7] = 0.125 * (1 - xw[j][0]) * (1 + xw[j][1]) * (1 + xw[j][2]);

    //natural derivatives c++ style (i*ndim+j convention):
    natDev[0] = -0.125 * (1 - xw[j][1]) * (1 - xw[j][2]);
    natDev[3] = 0.125 * (1 - xw[j][1]) * (1 - xw[j][2]);
    natDev[6] = 0.125 * (1 + xw[j][1]) * (1 - xw[j][2]);
    natDev[9] = -0.125 * (1 + xw[j][1]) * (1 - xw[j][2]);
    natDev[12] = -0.125 * (1 - xw[j][1]) * (1 + xw[j][2]);
    natDev[15] = 0.125 * (1 - xw[j][1]) * (1 + xw[j][2]);
    natDev[18] = 0.125 * (1 + xw[j][1]) * (1 + xw[j][2]);
    natDev[21] = -0.125 * (1 + xw[j][1]) * (1 + xw[j][2]);
    natDev[1] = -0.125 * (1 - xw[j][0]) * (1 - xw[j][2]);
    natDev[4] = -0.125 * (1 + xw[j][0]) * (1 - xw[j][2]);
    natDev[7] = 0.125 * (1 + xw[j][0]) * (1 - xw[j][2]);
    natDev[10] = 0.125 * (1 - xw[j][0]) * (1 - xw[j][2]);
    natDev[13] = -0.125 * (1 - xw[j][0]) * (1 + xw[j][2]);
    natDev[16] = -0.125 * (1 + xw[j][0]) * (1 + xw[j][2]);
    natDev[19] = 0.125 * (1 + xw[j][0]) * (1 + xw[j][2]);
    natDev[22] = 0.125 * (1 - xw[j][0]) * (1 + xw[j][2]);
    natDev[2] = -0.125 * (1 - xw[j][0]) * (1 - xw[j][1]);
    natDev[5] = -0.125 * (1 + xw[j][0]) * (1 - xw[j][1]);
    natDev[8] = -0.125 * (1 + xw[j][0]) * (1 + xw[j][1]);
    natDev[11] = -0.125 * (1 - xw[j][0]) * (1 + xw[j][1]);
    natDev[14] = 0.125 * (1 - xw[j][0]) * (1 - xw[j][1]);
    natDev[17] = 0.125 * (1 + xw[j][0]) * (1 - xw[j][1]);
    natDev[20] = 0.125 * (1 + xw[j][0]) * (1 + xw[j][1]);
    natDev[23] = 0.125 * (1 - xw[j][0]) * (1 + xw[j][1]);

};

void classLinearHexahedra::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                      int bulkElement) {
    ERROR("should not be used");
    exit(-1);
};

/*! This function defines just the positions, weights and Jacobians of the GPs inside the element. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] orderd This value defines the number of GPs inside the element (dramatically important for MM simulations)
  @param[inout] mmGPs Array with all GPs in the MM domain
*/
void classLinearHexahedra::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                                int orderd, vector<classGPs *> &mmGPs) {

    int nGPs = 0;
    int order = 2;
    if (order == 1) {
        nGPs = 1;
    } else if (order == 2) {
        nGPs = 8;
    } else {
        ERROR("The integration order %i for Hexahedra is not available. Orders 1 and 2 are the ones available in MuPhiSim",
              order);
        exit(-1);
    }
    //double xi[nGPs],eta[nGPs],zeta[nGPs], weights[nGPs];
    double xw[8][4];//1gp, 3dim+w
    GPsHexahedron3DValues(order, nGPs, xw);//xi, eta, zeta, weights);
    classGPs *tempMMGPs;
    vector<double> xyzIn(ndim, 0);
    vector<int> neighbours; // Temporal list of nodes in which a given GPs lie
    int numNodes = this->myNodes.size();
    neighbours = myNodes;
    vector<double> xyzCoor(ndim * numNodes, 0);
    //double x[4]={},y[4]={}, z[4]={}; // They are all of them triangles at this point

    //for(int j=0;j<numNodes;j++){
    //  xyz=nodes[myNodes[j]]->getXYZ();
    //  x[j]=xyz[0];
    //  y[j]=xyz[1];
    //  z[j]=xyz[2];
    //}
    vector<vector<double>> xyz;
    for (int i = 0; i < numNodes; i++) {
        vector<double> xyzTemp = nodes[myNodes[i]]->getXYZ();
        xyz.push_back(xyzTemp);
        xyzCoor[i] = xyzTemp[0];
        xyzCoor[numNodes + i] = xyzTemp[1];
        xyzCoor[2 * numNodes + i] = xyzTemp[2];
    }


    // https://people.fh-landshut.de/~maurer/numeth/node148.html
    //for(int k=0;k<numNodes;k++){
    //  xyzCoor[k]=x[k];
    //  xyzCoor[numNodes+k]=y[k];
    //  xyzCoor[2*numNodes+k]=z[k];
    //}

    vector<double> Jacobian(ndim * ndim, 0);

    vector<double> invJacobian(ndim * ndim, 0);
    vector<double> natDev;
    natDev.resize(myNodes.size() * ndim, 0);


    for (int j = 0; j < nGPs; j++) {
        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * ndim, 0);
        vector<double> Jacobian(ndim * ndim, 0);
        vector<double> invJacobian(ndim * ndim, 0);
        getInterpolationParameters(phi, natDev, xw, j);
        multGeneralMatrices(xyzCoor, ndim, numNodes, natDev, numNodes, ndim, Jacobian);
        double J = determinantTensor(Jacobian, ndim);
        inverse(Jacobian, ndim, invJacobian);


        vector<double> xyzIn(localNdim, 0);
        for (int i = 0; i < ndim; i++) {
            for (int k = 0; k < numNodes; k++) {
                xyzIn[i] += xyz[k][i] * phi[k];
            }
        }

        double weight;
        weight = xw[j][3];// It is already divided by 1/6
        tempMMGPs = new classGPs(GPs.size(), xyzIn, weight, J, ndim);
        tempMMGPs->setCell(ind, this->getType());
        tempMMGPs->setConsModel(constitutiveModel);

        mmGPs.push_back(tempMMGPs);
        GPs.push_back(tempMMGPs);
    }

}
