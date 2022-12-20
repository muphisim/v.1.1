//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classQuadLine_H
#define _classQuadLine_H

#include "classElements.h"

/*! \brief  This class is for linear 1D line */
class classQuadLine : public classElements {

protected:
    vector<double> n0;


public:

    classQuadLine(int me, string type, const vector<int>& myNodes, int ndim, vector<classNodes *>& nod, int bulkElement);

    virtual vector<double> getn0(vector<classNodes *>& nod, int ndim);

    virtual void
    definePositionGPsMM(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int orderd,
                        vector<classGPs *> &mmGPs);

    virtual double characteristicLength(int ndim, const vector<vector<double> >& nodesVec);

    virtual void getInterpolationParameters(double phi[], vector<double> &natDev, double xw[][2], int j);

    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim);

    virtual void defineGPs(int ind, vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, int bulkElement);
    virtual void getSubElement(int me, int S, vector<classElements *>& subElements, vector<classNodes *>& nod, int bulkElement,
                  int surface_order)
    {
        ERROR("classQuadLine::getSubElement has not been defined");
        exit(-1);
    }
                  
    virtual void defineShapeOnAPoint(int ind, classGPs *&tempGP, vector<classNodes *> &nodes, int ndim, const vector<double>& Point)
    {
        ERROR("classQuadLine::defineShapeOnAPoint has not been defined");
        exit(-1);
    }
    ~classQuadLine() {
    };
};

#endif
