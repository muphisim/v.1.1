//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classLinearSquare_H
#define _classLinearSquare_H

#include "classElements.h"

/*! \brief  This class is for linear squares */
class classLinearSquare : public classElements {

protected:
    vector<double> n0;


public:

    classLinearSquare(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod, int bulkElement);

    virtual vector<double> getn0(vector<classNodes *>& nod, int ndim);

    virtual void
    getSubElement(int me, int S, vector<classElements *> &subElements, vector<classNodes *>& nod, int bulkElement,
                  int surface_order);

    virtual void
    definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int orderd,
                        vector<classGPs *> &mmGPs);

    virtual double characteristicLength(int ndim, const vector<vector<double> >& nodesVec);

    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim);

    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int bulkElement);

    void getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][3], int j);

    virtual void
    defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim, const vector<double>& Point);

    ~classLinearSquare() {
    };
};

#endif
