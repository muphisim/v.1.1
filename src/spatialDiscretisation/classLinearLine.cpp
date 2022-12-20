//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classLinearLine.cpp
  \brief This file contains the definition of the functions for classLinearLines
*/

#include "classLinearLine.h"

/*! This function calculates the normal of the segment based on the two points
  @param[in] nod Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
vector<double> classLinearLine::getn0(vector<classNodes *>& nod, int ndim) {
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
classLinearLine::classLinearLine(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
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
double classLinearLine::characteristicLength(int ndim, const vector<vector<double> >& nodesVec) {

    const vector<double>& xyz1 = nodesVec[0];
    const vector<double>& xyz2 = nodesVec[1];
    double dist = calcDist(xyz1, xyz2, ndim);
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
void classLinearLine::definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim,
                                          int orderd, vector<classGPs *> &mmGPs) {
//GPs of subelements are defined using defineGPs function (as in FEM)
    ERROR("How did you even get here?");
    exit(-1);
};

/*! This function fully defines the GPs inside the element (including the value of the shape functions in case FEM. ALL ELEMENTS NEW ELEMENTS SHOULD IMPLEMENT THIS FUNCTION
  @param[in] ind Label of th element
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
*/
void classLinearLine::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim) {
    // int localNdim=this->localNdim;
    double xw[1][2];//nr of gps, nr of abcissa+1weight/gp
    int order = 1;
    int nGPs = 0; // Number of GPs per elements. It is not the total number of GPs
    GPsLine1DValues(order, nGPs, xw); // Integration order, not elements order
    vector<int>& neighbours = this->myNodes;
    //vector<double> xyz;
    //vector<double> xyz2;
    int numNodes = neighbours.size();
    vector<vector<double> > xyz;
    vector<double> xyzCoorLocal(localNdim *numNodes,
    0);
    for (int i = 0; i < numNodes; i++) {
        const vector<double>& xyzTemp = nodes[myNodes[i]]->getXYZ();
        xyz.push_back(xyzTemp);
    }
    //define new coordinate basis on the vector formed by the two points for proper ndim reduction
    xyzCoorLocal[0] = 0;
    xyzCoorLocal[1] = calcDist(xyz[0], xyz[1], 2);
    //xyz=nodes[myNodes[0]]->getXYZ();
    //xyz2=nodes[myNodes[1]]->getXYZ();
    //double J=calcDist(xyz, xyz2, 2);
    vector<double> xyzIn(ndim, 0);
    for (int j = 0; j < nGPs; j++) {
        //double xi;
        //xi=xw[j][0];
        double phi[numNodes];     // basis functions
        vector<double> natDev(numNodes * localNdim, 0);
        vector<double> Jacobian(localNdim *localNdim,
        0);
        vector<double> invJacobian(localNdim *localNdim,
        0);
        getInterpolationParameters(phi, natDev, xw, j);
        multGeneralMatrices(xyzCoorLocal, localNdim, numNodes, natDev, numNodes, localNdim, Jacobian);
        double J = determinantTensor(Jacobian, localNdim);
        //phi[0]=0.5-0.5*xi;
        //phi[1]=0.5+0.5*xi;
        for (int i = 0; i < ndim; i++) {
            for (int k = 0; k < numNodes; k++) {
                xyzIn[i] += xyz[k][i] * phi[k];
            }
        }
        //xyzIn[0] = xyz[0]*phi[0]+xyz2[0]*phi[1];
        //xyzIn[1] = xyz[1]*phi[0]+xyz2[1]*phi[1];
        double weight = 0.5 * xw[j][1];
        classGPs *tempGP = new classGPs(GPs.size(), xyzIn, weight, J, ndim);
        int numNeighbours = neighbours.size();
        tempGP->setNeighbours(neighbours, numNeighbours);
        tempGP->setCell(ind, this->getType());
        tempGP->setConsModel(constitutiveModel);
        double dphi[localNdim * numNeighbours];   // first derivatives of basis functions
        double ddphi[localNdim * localNdim * numNeighbours];// second derivatives of basis functions
        //vector <double> natDev;
        //natDev.resize(myNodes.size()*ndim,0);
        //for(int i = 0; i<myNodes.size()*ndim; i++)
        //{
        //  natDev[i]=natDevTensor[j][i];//define natural derivative for current gauss point
        //}
        vector<double> result(numNeighbours * localNdim, 0);
        multGeneralMatrices(natDev, numNeighbours, localNdim, invJacobian, localNdim, localNdim, result);
        for (int i = 0; i < numNeighbours; i++) {
            for (int j = 0; j < localNdim; j++) {
                dphi[i * localNdim + j] = result[i * localNdim + j];
            }
        }
        //dphi[0] = natDev[0]/J;
        //dphi[1] = natDev[1]/J;

        tempGP->setShapeFunctions(phi, dphi, ddphi, ndim, numNeighbours);
        GPs.push_back(tempGP);
    }
};

void
classLinearLine::defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int _bulkElement) {
    int currentIndex = GPs.size();
    defineGPs(ind, GPs, nodes, ndim);
    for (int j = currentIndex; j < GPs.size(); j++) {
        GPs[j]->setBulkElement(_bulkElement);
    }
};


void classLinearLine::getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][2], int j) {
    //shape functions c++ style:
    phi[0] = 0.5 * xw[j][0] + 0.5;
    phi[1] = -0.5 * xw[j][0] + 0.5;


    //natural derivatives c++ style (i*ndim+j convention):
    natDev[0] = 0.5 + 0.0;
    natDev[1] = -0.5 + 0.0;

};
