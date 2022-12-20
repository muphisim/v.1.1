//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classLinearTetrahedra.cpp
  \brief This file contains the definition of the functions used of the classLinearTetrahedra
*/

#include "classLinearTetrahedra.h"

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classLinearTetrahedra::getn0(vector<classNodes *>& nod, int ndim) {
    ERROR("Tetrahedra cannot be a subelement for the application of the pressure BCs");
    exit(-1);
}

/*! This function calculates a subelement from a given element. This function should not be called for segments
  @param[in] me Label of the element
  @param[in] S label of the edge (2D) or surface (3D) that is the subelement
  @param[out] subElement Subelement
  @param[in] nod Array with all nodes in the domain
*/
void classLinearTetrahedra::getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod,
                                          int bulkElement, int surface_order) {
    int numSs = myNodes.size();
    if (S > numSs) {
        ERROR("You are trying to access to a subElement that does not exist in your original element");
        exit(-1);
    }
    vector<int> mySubNodes;
    if (S == 1) {
        mySubNodes.push_back(myNodes[S - 1]);
        mySubNodes.push_back(myNodes[S]);
        mySubNodes.push_back(myNodes[S + 1]);
    } else if ((S == 2) || (S == 3)) {
        mySubNodes.push_back(myNodes[S - 2]);
        mySubNodes.push_back(myNodes[3]);
        mySubNodes.push_back(myNodes[S - 1]);
    } else {
        mySubNodes.push_back(myNodes[0]);
        mySubNodes.push_back(myNodes[2]);
        mySubNodes.push_back(myNodes[S - 1]);
    }

    //cout<< mySubNodes[0]<< " d "<< mySubNodes[1]<<endl;
    classElements *subElement2;
    subElement2 = new classLinearTriangles(me, "2dLinearTriangle", mySubNodes, 2, nod, bulkElement);
    subElements.push_back(subElement2);

}

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] tempGP Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] Point Point to be evaluated
*/
void classLinearTetrahedra::defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim,
                                                const vector<double>& Point) {
    vector<int> neighbours = this->myNodes;
    int numNodes = myNodes.size();
    vector<double> xyzCoor(ndim * numNodes, 0);
    double x[4] = {}, y[4] = {}, z[4] = {}; // They are all of them triangles at this point
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
    double xw[1][4];
    xw[0][0] = 0;
    xw[0][1] = 0;
    xw[0][2] = 0;
    xw[0][3] = 0;//does not matter, natDev is constant
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

/*! Constructor for classElements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] nod Array with all nodes in the domain
  @param[in] bulkElement index of mother element
*/
classLinearTetrahedra::classLinearTetrahedra(int me, string type, const vector<int>& myNodes, int ndim,
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
double classLinearTetrahedra::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    const vector<double>& p0 = nodesVec[0];
    const vector<double>& p1 = nodesVec[1];
    const vector<double>& p2 = nodesVec[2];
    const vector<double>& p3 = nodesVec[3];
    vector<double> tetraNodes;
    tetraNodes.resize(ndim * 4);
    for (int i = 0; i < ndim; i++) {
        tetraNodes[i * 4] = p0[i];
        tetraNodes[i * 4 + 1] = p1[i];
        tetraNodes[i * 4 + 2] = p2[i];
        tetraNodes[i * 4 + 3] = p3[i];

    }
    double radius = inscribedSphere(tetraNodes, ndim);
    return 2 * radius;
};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
void classLinearTetrahedra::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim) {

    int nGPs = 0;
    int order = 1;
    if (order == 1) {
        nGPs = 1;
    } else if (order == 2) {
        nGPs = 4;
    } else {
        ERROR("The integration order %i for tetrahedra is not available. Orders 1 and 2 are the ones available in MuPhiSim",
              order);
        exit(-1);
    }
    //double xi[nGPs],eta[nGPs],zeta[nGPs], weights[nGPs];
    double xw[1][4];//1gp, 3dim+w
    GPsTetrahedron3DValues(order, nGPs, xw);//xi, eta, zeta, weights);

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
        if (flagInit) { 
            this->addMyGPs(GPs.size());
        } //This vector is a indices list from 0 to maximum number of GPs
        GPs.push_back(tempGP);
    }
};

void classLinearTetrahedra::getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][4], int j) {
    //shape functions c++ style:
    phi[0] = -1.0 * xw[j][0] - 1.0 * xw[j][1] - 1.0 * xw[j][2] + 1.0;
    phi[1] = 1.0 * xw[j][0];
    phi[2] = 1.0 * xw[j][1];
    phi[3] = 1.0 * xw[j][2];


    //natural derivatives c++ style (i*ndim+j convention):
    natDev[0] = -1.0 - 0.0 - 0.0 + 0.0;
    natDev[3] = 1.0;
    natDev[6] = 0.0;
    natDev[9] = 0.0;
    natDev[1] = -0.0 - 1.0 - 0.0 + 0.0;
    natDev[4] = 0.0;
    natDev[7] = 1.0;
    natDev[10] = 0.0;
    natDev[2] = -0.0 - 0.0 - 1.0 + 0.0;
    natDev[5] = 0.0;
    natDev[8] = 0.0;
    natDev[11] = 1.0;

};

void classLinearTetrahedra::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
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
void classLinearTetrahedra::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                                int orderd, vector<classGPs *> &mmGPs) {

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
    //double xi[nGPs],eta[nGPs],zeta[nGPs], weights[nGPs];
    double xw[4][4];//1gp, 3dim+w
    GPsTetrahedron3DValues(order, nGPs, xw);//xi, eta, zeta, weights);
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
