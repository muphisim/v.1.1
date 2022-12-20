//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classLinearHexahedra_H
#define _classLinearHexahedra_H

#include "classElements.h"

/*! \brief  This class is for linear hexahedra */
class classLinearHexahedra : public classElements {

public:

    classLinearHexahedra(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod,
                          int bulkElement);

    virtual vector<double> getn0(vector<classNodes *>& nod, int ndim);

    virtual void getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod, int bulkElement,
                  int surface_order);

    virtual double characteristicLength(int ndim, const vector<vector<double> >& nodesVec);

    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim);

    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int bulkElement);

    virtual void getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][4], int j);

    virtual void definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int orderd,
                        vector<classGPs *> &mmGPs);

    virtual void defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim, const vector<double>& Point);

    ~classLinearHexahedra() {
    };
};


#endif
