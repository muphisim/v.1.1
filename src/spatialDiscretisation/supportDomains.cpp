//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file supportDomains.cpp
  \brief This file contains all functions related to the calculation of the support domains as well as the functions to evaluate where each GP lies
*/
#include "supportDomains.h"
#include "commandLine.h"


template<typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}


/*! \brief This function defines the radius of the support domain of each node according to the distance between the nodes
  @param[in] nodesXYZ Array with the nodes coordinates
  @param[in] ndim Dimension of the domain
  @param[in] numNodes Total number of nodes
  @param[out] rmax Array with the radius of the support domain of each node
  @param[out] charlengthscale Typical lenght need in Sukumar's formulation
  @param[in] S0 number of nodes to average the distance
*/

extern void
supportDomain(double *nodesXYZ, int ndim, int numNodes, vector<double> &rmax, double &charlengthscale, int S0) {
    if (S0 <= 1E-16) {
        S0 = 3;// Can be changed, e.g., 2
    }
    double *me, *you;
    vector<double> euclidean;
    double temp = 0, temp2 = 0;
    double dot = 0;
    double dm = GeneralOptions::MMScalingFactor;// Can be changed, e.g., 2, introduce in comandlind with the key -MMScalingFactor followed by a number
    if (S0 > numNodes) { S0 = numNodes; }
    charlengthscale = 0.;
    for (int j = 0; j < numNodes; j++) {
        // For each node I will compute the rmax, which is the radious of the support domain. rmax=dmax*cI, where cI is the characteristic lenght of the neighbours (juli algorithm)
        // Compute the euclidean distance of all nodes respect to j (brute force until now, that should be optimised)
        // Euclidean distance of all nodes
        me = &nodesXYZ[j * ndim];
        for (int i = 0; i < numNodes; i++) {
            you = &nodesXYZ[i * ndim];
            dot = 0;
            for (int kk = 0; kk < ndim; kk++) {
                dot += (me[kk] - you[kk])*(me[kk] - you[kk]);
            }
            euclidean.push_back(sqrt(dot));
            // cout<< "j "<< j<< " i="<<i << " num "<< numNodes<< " rmax= "<< euclidean[i] << endl;
        }
        vector<size_t> idx;
        idx = sort_indexes(euclidean);
        // sum up the distance to the S0 closest neighbours and calculated the weighted average
        temp = 0;
        for (int i = 0; i < S0; i++) {
            temp = temp + euclidean[idx[i + 1]];    // that will be cI
        }

        temp2 = temp / double(S0);
        euclidean.clear();
        rmax[j] = dm * temp2;
        // This is the scale of the domain. Something to help the units independency
        charlengthscale = charlengthscale + rmax[j];
        //cout<< "i="<<j << "rmax= "<< rmax[j] << endl;
    }
    // This is the scale of the domain. Something to help the units independency
    charlengthscale = charlengthscale / numNodes;
}


/*! \brief This function finds the nodes in which the GP lies
  @param[in] GP Gauss point coordinates
  @param[in] me Label of the GP   
  @param[in] nodesXYZ Array with the nodes coordinates
  @param[in] ndim Dimension of the domain
  @param[in] numNodes Total number of nodes
  @param[out] neighbours Array with the neighboring nodes
  @param[in] rmax Array with the radius of the support domain of each node
*/
extern void findNeighbours(vector<double> GP, int me, double *nodesXYZ, int ndim, int numNodes, vector<int> &neighbours,
                           const vector<double> &rmax) {

    vector<double> A(ndim, 0);
    vector<double> I(ndim * ndim, 0);
    vector<double> result(ndim, 0);
    double dot = 0;
    vector<double> q(numNodes, 0);
    for (int i = 0; i < ndim; i++) {
        I[i * ndim + i] = 1;
    }
    for (int i = 0; i < numNodes; i++) {
        q[i] = 0;
        for (int j = 0; j < ndim; j++) {
            A[j] = nodesXYZ[i * ndim + j] - GP[j];
            result[j] = 0;
        }
        multTensorVector(I, ndim, ndim, A, ndim, result);
        dot =0;
        for (int j = 0; j < ndim; j++) 
        {
            dot += A[j]*result[j];
        }
        q[i] = sqrt(dot) / rmax[i];
        if (q[i] <= 1) {
            neighbours.push_back(i);
        }
    }
}

