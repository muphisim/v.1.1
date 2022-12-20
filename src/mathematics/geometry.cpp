//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*! \file geometry.cpp
\brief This file contains all functions related to geometrical calculi
*/

#include "geometry.h"
#include "shapeFunctions.h"
#include "solvers.h"
#include "maths.h"


/*! \brief This function calculates the distance between two points in 2D
  @param[in] sX x coordinate of point 1
  @param[in] sY y coordinate of point 1
  @param[in] eX x coordinate of point 2
  @param[in] eY y coordinate of point 2
  @return Distance between the points 
*/
extern double calcDist(double sX, double sY, double eX, double eY) {
    eX -= sX;
    eY -= sY;
    return sqrt(eX * eX + eY * eY);
};

/*! \brief This function calculates the distance between two points in 3D
  @param[in] sX x coordinate of point 1
  @param[in] sY y coordinate of point 1
  @param[in] sZ z coordinate of point 1
  @param[in] eX x coordinate of point 2
  @param[in] eY y coordinate of point 2
  @param[in] eZ z coordinate of point 2
  @return Distance between the points 
*/
extern double calcDist(double sX, double sY, double sZ, double eX, double eY, double eZ) {
    eX -= sX;
    eY -= sY;
    eZ -= sZ;
    return sqrt(eX * eX + eY * eY + eZ * eZ);
};

/*! \brief This function calculates the distance between two points
  @param[in] s Array with the coordinates of point 1
  @param[in] e Array with the coordinates of point 2
  @param[in] ndim Dimension of the problem
  @return Distance between the points 
*/
extern double calcDist(const vector<double>& s, const vector<double>& e, int ndim) {
    if (ndim == 2) {
        double eX = e[0];
        double eY = e[1];
        double sX = s[0];
        double sY = s[1];
        eX -= sX;
        eY -= sY;
        return sqrt(eX * eX + eY * eY);
    } else {
        double eX = e[0];
        double eY = e[1];
        double eZ = e[2];
        double sX = s[0];
        double sY = s[1];
        double sZ = s[2];
        eX -= sX;
        eY -= sY;
        eZ -= sZ;
        return sqrt(eX * eX + eY * eY + eZ * eZ);
    }
    return -1;
};

void deformedPositions(vector<double> &posXYZ, const vector<double>& disp, int idx, int ndim) {
    int nodeNum = disp.size() / ndim;
    for (int i = 0; i < ndim; i++) {
        posXYZ[i] += disp[nodeNum * i + idx];
    }
}

/*! \brief This calculates the area of a triangle in 2D simulations and 3D simulations
  @param[in] nodes Array with the coordinates of the vertices of the triangle (x1,x2,x3; y1,y2,y3; z1, z2, z3)
  @param[in] ndim Dimension of the problem
  @return Area of the triangle
*/
extern double areaTriangle(const vector<double>& nodes, int ndim) {

    double area = -1.0;
    vector<double> matrix(9, 0);
    int nNod = nodes.size();
    if (ndim == 2) { // && ndim == 2, split into two cases ndim=2 and ndim=3
        if (nNod != 6) {
            ERROR("A triangle should have 3 vertices. You provided %f", double(nNod)/2);
            exit(-1);
        } else {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 2; j++) {
                    matrix[i + j * 3] = nodes[i + j * 3];
                }
            }
            matrix[6] = 1;
            matrix[7] = 1;
            matrix[8] = 1;

            area = fabs(determinantTensor(matrix, 3) / 2);

        }
    } else if (ndim == 3) {
        if (nNod != 9) {
            ERROR("A triangle should have 3 vertices. You provided %f", double(nNod)/3);
            exit(-1);
        } else {
            //use |AxB|/2, first get A=v1-v3, B=v2-v3
            vector<double> A;
            vector<double> B;
            for (int i = 0; i < 3; i++) {
                double v1MINv3 = nodes[i * 3] - nodes[i * 3 + 2];
                double v2MINv3 = nodes[i * 3 + 1] - nodes[i * 3 + 2];
                A.push_back(v1MINv3);
                B.push_back(v2MINv3);
            }
            //compute cross product square sum of components
            vector<double> C;
            crossProduct(A, B, C);
            area = 0.5 * sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
        }

    }
    return area;


};

/*! \brief This calculates the area of a square in 2D simulations and 3D simulations
  @param[in] nodes Array with the coordinates of the vertices of the triangle (x1,x2,x3; y1,y2,y3; z1, z2, z3)
  @param[in] ndim Dimension of the problem
  @return Area of the square
*/
extern double areaSquare(const vector<double>& squareNodes, int ndim) {
    vector<double> xyz0, xyz1, xyz2, xyz3, triangle_1, triangle_2;
    xyz0.resize(ndim, 0);
    xyz1.resize(ndim, 0);
    xyz2.resize(ndim, 0);
    xyz3.resize(ndim, 0);
    triangle_1.resize(3*ndim, 0);
    triangle_2.resize(3*ndim, 0);
    for (int i = 0; i < ndim; i++) {
        xyz0[i] = squareNodes[i * 4];
        xyz1[i] = squareNodes[i * 4 + 1];
        xyz2[i] = squareNodes[i * 4 + 2];
        xyz3[i] = squareNodes[i * 4 + 3];
    }

    for (int i = 0; i < ndim; i++) {
        triangle_1[i * 3] = squareNodes[i * 4];
        triangle_1[i * 3 + 1] = squareNodes[i * 4 + 1];
        triangle_1[i * 3 + 2] = squareNodes[i * 4 + 2];

        triangle_2[i * 3] = squareNodes[i * 4];
        triangle_2[i * 3 + 1] = squareNodes[i * 4 + 2];
        triangle_2[i * 3 + 2] = squareNodes[i * 4 + 3];
    }
    double total_area;
    double area_1 = areaTriangle(triangle_1, ndim);
    double area_2 = areaTriangle(triangle_2, ndim);
    total_area = area_1 + area_2;
    return total_area;


};
/*! \brief This calculates the cross product of two vectors
  @param[in] A Array with the first vector
  @param[in] B Array with the second vector
  @param[out] C Array with output vector
*/
extern void crossProduct(const vector<double> &A, const vector<double> &B, vector<double> &C) {
    if (A.size() != 3 || B.size() != 3) {
        ERROR("Cross product is only defined for R3");
        exit(-1);
    } else {
        C.resize(3, 0);
        C[0] = A[1] * B[2] - A[2] * B[1];
        C[1] = A[2] * B[0] - A[0] * B[2];
        C[2] = A[0] * B[1] - A[1] * B[0];
    }
};

/*! \brief This calculates the radius of the circle inscribed in a triangle
  @param[in] nodes Array with the coordinates of the vertices of the triangle (x1,x2,x3; y1,y2,y3; z1, z2, z3)
  @param[in] ndim Dimension of the problem
  @return radius of the inscribed circle
  */
extern double inscribedCircle(const vector<double> &triangleNodes, int ndim) {
    if (triangleNodes.size() != 3 * ndim) {
        ERROR("Triangles have 3 nodes!");
        exit(-1);
    }
    vector<double> xyz0, xyz1, xyz2;
    xyz0.resize(ndim, 0);
    xyz1.resize(ndim, 0);
    xyz2.resize(ndim, 0);
    for (int i = 0; i < ndim; i++) {
        xyz0[i] = triangleNodes[i * 3];
        xyz1[i] = triangleNodes[i * 3 + 1];
        xyz2[i] = triangleNodes[i * 3 + 2];
    }
    double area = areaTriangle(triangleNodes, ndim);
    double perimeter = calcDist(xyz0, xyz1, ndim) + calcDist(xyz1, xyz2, ndim) + calcDist(xyz2, xyz0, ndim);
    double r = 2 * area / perimeter;
    return r;
};

/*! \brief This calculates the radius of the circle inscribed in a square
  @param[in] nodes Array with the coordinates of the vertices of the triangle (x1,x2,x3; y1,y2,y3; z1, z2, z3)
  @param[in] ndim Dimension of the problem
  @return radius of the inscribed circle
  */
extern double inscribedCircleQuad(const vector<double> &squareNodes, int ndim) {
    if (squareNodes.size() != 4 * ndim) {
        ERROR("Square have 4 nodes!");
        exit(-1);
    }
    vector<double> xyz0, xyz1, xyz2, xyz3, triangle_1, triangle_2;
    xyz0.resize(ndim, 0);
    xyz1.resize(ndim, 0);
    xyz2.resize(ndim, 0);
    xyz3.resize(ndim, 0);
    triangle_1.resize(3*ndim, 0);
    triangle_2.resize(3*ndim, 0);
    for (int i = 0; i < ndim; i++) {
        xyz0[i] = squareNodes[i * 4];
        xyz1[i] = squareNodes[i * 4 + 1];
        xyz2[i] = squareNodes[i * 4 + 2];
        xyz3[i] = squareNodes[i * 4 + 3];
    }

    for (int i = 0; i < ndim; i++) {
        triangle_1[i * 3] = squareNodes[i * 4];
        triangle_1[i * 3 + 1] = squareNodes[i * 4 + 1];
        triangle_1[i * 3 + 2] = squareNodes[i * 4 + 2];

        triangle_2[i * 3] = squareNodes[i * 4];
        triangle_2[i * 3 + 1] = squareNodes[i * 4 + 2];
        triangle_2[i * 3 + 2] = squareNodes[i * 4 + 3];
    }
    double total_area;
    double area_1 = areaTriangle(triangle_1, ndim);
    double area_2 = areaTriangle(triangle_2, ndim);
    total_area = area_1 + area_2;
    double perimeter = calcDist(xyz0, xyz1, ndim) + calcDist(xyz1, xyz2, ndim) + calcDist(xyz2, xyz3, ndim) + calcDist(xyz3, xyz0, ndim);
    double r = 2 * total_area / perimeter;
    return r;
};

/*! \brief This calculates the radius of the sphere inscribed in a tetrahedron
  @param[in] nodes Array with the coordinates of the vertices of the tetrahedra (x1,x2,x3,x4; y1,y2,y3,y4; z1, z2, z3,z4)
  @param[in] ndim Dimension of the problem
  @return radius of the inscribed circle
  */
extern double inscribedSphere(const vector<double> &tetraNodes, int ndim) {
    if (tetraNodes.size() != 4 * ndim) {
        ERROR("Tetrahedra have 4 nodes!");
        exit(-1);
    }
    vector<double> p0, p1, p2, p3;
    p0.resize(ndim, 0);
    p1.resize(ndim, 0);
    p2.resize(ndim, 0);
    p3.resize(ndim, 0);
    for (int i = 0; i < ndim; i++) {
        p0[i] = tetraNodes[i * 4];
        p1[i] = tetraNodes[i * 4 + 1];
        p2[i] = tetraNodes[i * 4 + 2];
        p3[i] = tetraNodes[i * 4 + 3];

    }
    double s1, s2, s3, s4;
    vector<double> areaArgument = {p0[0], p1[0], p2[0], p0[1], p1[1], p2[1], p0[2], p1[2], p2[2]};
    s1 = areaTriangle(areaArgument, ndim); //area of triangle formed by these nodes
    areaArgument = {p1[0], p2[0], p3[0], p1[1], p2[1], p3[1], p1[2], p2[2], p3[2]};
    s2 = areaTriangle(areaArgument, ndim);
    areaArgument = {p0[0], p2[0], p3[0], p0[1], p2[1], p3[1], p0[2], p2[2], p3[2]};
    s3 = areaTriangle(areaArgument, ndim);
    areaArgument = {p0[0], p1[0], p3[0], p0[1], p1[1], p3[1], p0[2], p1[2], p3[2]};
    s4 = areaTriangle(areaArgument, ndim);

    double A = fabs(s1) + fabs(s2) + fabs(s3) + fabs(s4);

    // Volume of subelement can be computed using 1/6*|det(p1-p0,p2-p0,p3-p0)|
    vector<double> tempMat = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2],
                              p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2],
                              p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};
    double V = 1.0 / 6.0 * fabs(determinantTensor(tempMat, ndim));
    double r = 3.0 * V / A;

    return r;
};

/*! \brief This calculates the radius of the sphere inscribed in a hexahedron
  @param[in] nodes Array with the coordinates of the vertices of the tetrahedra (x1,x2,x3,x4; y1,y2,y3,y4; z1, z2, z3,z4)
  @param[in] ndim Dimension of the problem
  @return radius of the inscribed circle
  */
extern double inscribedSphereHex(const vector<double> &HexaNodes, int ndim) {
    if (HexaNodes.size() != 8 * ndim) {
        ERROR("Hexahedra have 8 nodes!");
        exit(-1);
    }
    vector<double> p0, p1, p2, p3, p4, p5, p6, p7;
    p0.resize(ndim, 0);
    p1.resize(ndim, 0);
    p2.resize(ndim, 0);
    p3.resize(ndim, 0);
    p4.resize(ndim, 0);
    p5.resize(ndim, 0);
    p6.resize(ndim, 0);
    p7.resize(ndim, 0);
    for (int i = 0; i < ndim; i++) {
        p0[i] = HexaNodes[i * 8];
        p1[i] = HexaNodes[i * 8 + 1];
        p2[i] = HexaNodes[i * 8 + 2];
        p3[i] = HexaNodes[i * 8 + 3];
        p4[i] = HexaNodes[i * 8 + 4];
        p5[i] = HexaNodes[i * 8 + 5];
        p6[i] = HexaNodes[i * 8 + 6];
        p7[i] = HexaNodes[i * 8 + 7];
    }
    double s1, s2, s3, s4, s5, s6;
    vector<double> areaArgument = {p0[0], p1[0], p2[0], p3[0], p0[1], p1[1], p2[1], p3[1], p0[2], p1[2], p2[2], p3[2]};
    s1 = areaSquare(areaArgument, ndim); //area of triangle formed by these nodes
    areaArgument = {p4[0], p5[0], p6[0], p7[0], p4[1], p5[1], p6[1], p7[1], p4[2], p5[2], p6[2], p7[2]};
    s2 = areaSquare(areaArgument, ndim);
    areaArgument = {p0[0], p1[0], p5[0], p4[0], p0[1], p1[1], p5[1], p4[1], p0[2], p1[2], p5[2], p4[2]};
    s3 = areaSquare(areaArgument, ndim);
    areaArgument = {p1[0], p2[0], p6[0], p5[0], p1[1], p2[1], p6[1], p5[1], p1[2], p2[2], p6[2], p5[2]};
    s4 = areaSquare(areaArgument, ndim);
    areaArgument = {p2[0], p3[0], p7[0], p6[0], p2[1], p3[1], p7[1], p6[1], p2[2], p3[2], p7[2], p6[2]};
    s5 = areaSquare(areaArgument, ndim);
    areaArgument = {p0[0], p3[0], p7[0], p4[0], p0[1], p3[1], p7[1], p4[1], p0[2], p3[2], p7[2], p4[2]};
    s6 = areaSquare(areaArgument, ndim);
    double A = fabs(s1) + fabs(s2) + fabs(s3) + fabs(s4) + fabs(s5) + fabs(s6);

    // Volume of subelement can be computed using 1/6*|det(p1-p0,p2-p0,p3-p0)|
    vector<double> tempMat = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2],
                              p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2],
                              p4[0] - p0[0], p4[1] - p0[1], p4[2] - p0[2]};
    double V = fabs(determinantTensor(tempMat, ndim));
    double r = 3.0 * V / A;

    return r;
};

/*! \brief creates a coordinate basis on a plane spanned by two vectors, meant for triangle rotation
  @param[in] a first vector making up the plane
  @param[in] b second vector
  @param[in] ndim Dimension of the problem
  @param[out] mapping matrix out
*/
extern void createOrthogonalBasisOnPlane(const vector<double> &a, const vector<double> &b, int ndim, vector<double> &map) {
    // This function returns a mapping matrix instead of mapped vector ie M and not Vlocal=M*V
    // as it is inefficient to recalculate the mapping matrix for each time a vector is to mapped.
    // for quadratic triangles that is 6 times
    map.resize((ndim - 1) * ndim, 0);//size 3x2 as we redice dimension by 1
    if (ndim != 3) {
        INFO("Mapping from R3 to R2, you provided R%i", ndim);
        exit(-1);
    }
    //get normal
    vector<double> normal(ndim, 0);
    vector<double> xhat(ndim, 0);
    vector<double> yhat(ndim, 0);
    crossProduct(b, a, normal);
    //get x
    crossProduct(a, normal, xhat);
    double magA = norm2(a, ndim);//yhat is parallel to a
    double magX = norm2(xhat, ndim);
    for (int i = 0; i < ndim; i++) {
        yhat[i] = a[i] / magA;
        xhat[i] = xhat[i] / magX;
    }
    for (int i = 0; i < ndim; i++) {
        map[i] = xhat[i];
        map[ndim + i] = yhat[i];
    }

}

/*! \brief This calculates the signed area of a triangle in 2D simulations
  @param[in] nodes Array with the coordinates of the vertices of the triangle (x1,x2,x3; y1,y2,y3)
  @param[in] ndim Dimension of the problem
  @return Signed area of the triangle
*/
extern double signedAreaTriangle(const vector<double>& nodes, int ndim) {
    vector<double> matrix(9, 0);
    int nNod = nodes.size();
    if (nNod != 6) {
        ERROR("A triangle should have 3 vertices. You provided %f", double(nNod)/2);
        exit(-1);
    } else if (ndim == 3) {
        ERROR("This function is just valid for 2D simulations");
        exit(-1);
    }


    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            matrix[i + j * 3] = nodes[i + j * 3];
        }
    }
    matrix[6] = 1;
    matrix[7] = 1;
    matrix[8] = 1;
    double dete = determinantTensor(matrix, 3);
    return dete / fabs(dete);
};


/*! \brief Computation of the global index of a local node
  @param[in] nNodGlobal Number of nodes in the all domain
  @param[in] nNodLocal Number of nodes in the current processor part
  @param[in] dofLocal Index of the node in local part
  @param[in] nod_local List of global indexes for local nodes
  @param[in] extra Boolean, 1 if there is extra DoF
  @param[out] dofGlobal Index of the node in the global domain
*/
extern void
locDofToGlobDof(int &nNodGlobal, int &nNodLocal, int &dofLocal, int &dofGlobal,const vector<int> &nod_local, bool extra) {
    div_t divresult;
    int nod;
    int dim;
    int nodGlobal;
    divresult = div(dofLocal, nNodLocal);
    nod = divresult.rem;
    dim = divresult.quot + 1;
    nodGlobal = nod_local[nod];
    if (extra) {
        dofGlobal = nodGlobal;
    } else {
        dofGlobal = nNodGlobal * (dim - 1) + nodGlobal;
    }
}
