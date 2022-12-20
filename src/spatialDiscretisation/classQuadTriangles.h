//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classQuadTriangles_H
#define _classQuadTriangles_H

#include "classElements.h"

/*! \brief  This class is for quadratic triangles */
class classQuadTriangles : public classElements {

protected:
    vector<double> n0;


public:

    classQuadTriangles(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod, int bulkElement);

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

    virtual void getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][3], int j);

    virtual void
    defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim, const vector<double>& Point);

    ~classQuadTriangles() {
    };
};

#endif
