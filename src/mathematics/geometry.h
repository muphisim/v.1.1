//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _geometry_H_
#define _geometry_H_

#include "configuration.h"

extern void createOrthogonalBasisOnPlane(const vector<double> &a, const vector<double> &b, int ndim, vector<double> &map);

extern double signedAreaTriangle(const vector<double>& nodes, int ndim);

extern double areaTriangle(const vector<double>& nodes, int ndim);

extern double areaSquare(const vector<double>& nodes, int ndim);

extern double inscribedCircleQuad(const vector<double> &squareNodes, int ndim);

extern void crossProduct(const vector<double> &A, const vector<double> &B, vector<double> &C);

extern double inscribedCircle(const vector<double> &triangleNodes, int ndim);

extern double inscribedSphere(const vector<double> &tetraNodes, int ndim);

extern double calcDist(double sX, double sY, double eX, double eY);

extern double calcDist(double sX, double sY, double sZ, double eX, double eY, double eZ);

extern double calcDist(const vector<double>& s, const vector<double>& e, int ndim);

void deformedPositions(vector<double> &posXYZ, const vector<double>& disp, int idx, int ndim);

extern double inscribedSphereHex(const vector<double> &HexaNodes, int ndim);
#ifdef PARALLEL
extern void locDofToGlobDof(int &nNodGlobal, int &nNodLocal, int &dofLocal, int &dofGlobal, const vector<int> &nod_local, bool extra);
#endif

#endif
