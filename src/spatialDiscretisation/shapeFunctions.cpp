//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file shapeFunctions.cpp
  \brief Several functions to calculate the shape functions for MM
*/
#include "shapeFunctions.h"
#include "commandLine.h"

/*! \brief This function calculates the max-ent shape functions with the Sukumar approach (Sukumar and Wright 2007) at each GPs in the grid defined by "nodes"
  @param[inout] GPs Array with all GPs in the domain
  @param[in] nodes Array with all nodes in the domain
  @param[in] ndim Dimension of the domain
  @param[in] nodesMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[in] flagSupport It is true if the radius of the support domain of nodes has to be calculated
  @param[in] flagCurrCon It is true if the shape functions have to be calculated in the current configuration. Without remeshing, this flag is always false
*/

#ifdef MM           
extern void maxEntSFatGPs(vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, vector<int> &nodesMM, bool flagSupport,
              bool flagCurrCon) {

    int numNodes = nodesMM.size();
    vector<double> rmaxCurr(numNodes, 0);
    double charlengthscale = 1.;
    
    double nodesXYZCurr[numNodes * ndim];// Temporal array with the coordinates
    for (int j = 0; j < numNodes; j++) 
    {
        int ind = nodesMM[j];
        vector<double> xyz22;
        if (flagCurrCon) 
        {
            // current config
            xyz22 = (nodes[ind]->getXYZCurr());
        } 
        else 
        {
            // initial config
            xyz22 = (nodes[ind]->getXYZ());
        }
        for (int i=0; i< ndim; i++)
        {
            nodesXYZCurr[j * ndim+i] = xyz22[i];
        }
    }
    if (flagSupport) 
    {
        supportDomain(nodesXYZCurr, ndim, numNodes, rmaxCurr, charlengthscale, 0);
        for (int j = 0; j < nodesMM.size(); j++) 
        {
            int nod = nodesMM[j];
            nodes[nod]->setRmax(rmaxCurr[j]);
        }
    } 
    else 
    {
        charlengthscale = 0;
        for (int j = 0; j < nodesMM.size(); j++) 
        {
            int nod = nodesMM[j];
            rmaxCurr[j] = nodes[nod]->getRmax();
            charlengthscale += rmaxCurr[j];
            if (rmaxCurr[j] <= 0) {
                ERROR("The support domain of the nodes has not been properly defined");
                exit(-1);
            }
        }
        charlengthscale /= (double)numNodes;
    }

    
    char priorweight[80] = "cubic"; // uniform, cubic, quartic, gaussian, or gaussian-rbf
    char solnscheme[80] = "lbfgs";  // `newton', ``descent', 'lbfgs'
    char solnschemef[80];
    char priorweightf[80];
    ConvertToFortran(solnschemef, sizeof solnschemef, solnscheme);
    ConvertToFortran(priorweightf, sizeof priorweightf, priorweight);
    
    int maxit = 100;
    double eps = 1E-8, dettol = 1E-15;
    int printflag = 0;
    int onlyPhi = 0;
    int ierrorflag;
    int maxFails=100;
    int iter=0;
    bool flagShapeProblem = false;
    bool flagLocalProblem = false;
    while (!flagShapeProblem) 
    {
        iter++;
        if (iter > maxFails)
        {
            ERROR("The maximal number of iterations when creating shape functions is reached!");
            exit(-1);
        }
        // To calculate the shape functions for each evaluation point
        flagLocalProblem = false;
        for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
        {
            classGPs *tempGP = *ite;
            vector<double> xyz;
            if (flagCurrCon) 
            {
                xyz = tempGP->getXYZCurr();
            } 
            else 
            {
                xyz = tempGP->getXYZ();
            }
            vector<int> neighboursLocal(0);  // neighbours indexes in the MM list of nodes        
            findNeighbours(xyz, tempGP->getMe(), nodesXYZCurr, ndim, numNodes, neighboursLocal, rmaxCurr);
            vector<int> neighbours(neighboursLocal.size(),0);
            for (int j = 0; j < neighboursLocal.size(); j++) 
            {
                neighbours[j] = nodesMM[neighboursLocal[j]];
            }
            int numNeighbours = neighbours.size();
            double nodesXYZNeighbours[numNeighbours * ndim]; //nodal coordinates of the n neighbors
            vector<double> rmaxLocal(numNeighbours, 0);
            for (int j = 0; j < numNeighbours; j++) 
            {
                const vector<double>& xzyPoint= nodes[neighbours[j]]->getXYZ();
                for (int i=0; i< ndim; i++)
                {
                    nodesXYZNeighbours[j+i*numNeighbours] = xzyPoint[i];
                }
                rmaxLocal[j] = rmaxCurr[neighboursLocal[j]];
            }
            
            double D[ndim*ndim*numNeighbours]; //anisotropic metric = I when the isotropic case is used
            fill(D,D+ndim*ndim*numNeighbours,0.);
            for (int i=0; i< ndim; i++)
            {
                for (int j=0; j< numNeighbours; j++)
                {
                    D[i+ndim*i+ndim*ndim*j] = 1.;
                }
            }
            
            double phi[numNeighbours];     // basis functions
            double dphi[ndim * numNeighbours];   // first derivatives of basis functions
            double ddphi[ndim * ndim * numNeighbours];// second derivatives of basis functions
            fill(phi,phi+numNeighbours,0.);
            fill(dphi,dphi+ndim*numNeighbours,0.);
            fill(ddphi,ddphi+ndim*ndim*numNeighbours,0.);
            
            __maxent_MOD_drivermaxent(&numNeighbours, &ndim, &solnschemef[0], &priorweightf[0], nodesXYZNeighbours, &xyz[0],
                                      &rmaxLocal[0], D, &maxit, &eps, &printflag, &ierrorflag, &charlengthscale,
                                      &dettol, phi, dphi, ddphi);
            for (int i = 0; i < numNeighbours; i++) 
            {
                if ((!std::isfinite(phi[i])) or (ierrorflag != 0))
                {
                    double factorRmax = 1.1;
#ifdef PARALLEL
                    WARNING("rank =%d: Convex Hull problems GPs %i. Trying another RMAX, factor=%f", GeneralOptions::commRank, tempGP->getMe(), factorRmax);
#else
                    WARNING("Convex Hull problems GPs %i. Trying another RMAX, factor=%f", tempGP->getMe(), factorRmax);
#endif
                    if (!(GeneralOptions::MMAdaptiveSupportRadius))
                    {
                        ERROR("stop simulation!");
                        exit(-1);
                    }
                    for (int j = 0; j < nodesMM.size(); j++) {
                        int nod = nodesMM[j];
                        rmaxCurr[j] = rmaxCurr[j] * factorRmax;
                        nodes[nod]->setRmax(rmaxCurr[j]);
                    }
                    flagLocalProblem = true;
                    break;
                }
            }
            if (flagLocalProblem) 
            {
                break;
            }

            tempGP->setNeighbours(neighbours, numNeighbours);
            tempGP->setShapeFunctions(phi, dphi, ddphi, ndim, numNeighbours);
            
        }
        if (!flagLocalProblem) 
        {
            flagShapeProblem = true;
        }
    }
};

/*! \brief This function calculates the max-ent shape functions with the Sukumar approach (Sukumar and Wright 2007) at each node in the grid defined by "nodes2" (themselves)
  @param[inout] nodes Array with grid
  @param[in] nodes2 Array with the targets
  @param[in] ndim Dimension of the domain
  @param[in] nodesMM Array with the labels that correspond to the MM nodes in the nodes general array
  @param[in] flagSupport It is true if the radius of the support domain of nodes has to be calculated
  @param[in] flagCurrCon It is true if the shape functions have to be calculated in the current configuration. Without remeshing, this flag is always false
*/
extern void maxEntSFatNod(vector<classNodes *> &nodes, vector<classNodes *> &nodes2, int ndim, vector<int> &nodesMM,
                          bool flagSupport, bool flagCurrCon) {

    double factor = 1;
    int numNodes = nodesMM.size();
    vector<double> rmaxCurr(numNodes, 0);
    char priorweight[80] = "cubic"; // uniform, cubic, quartic, gaussian, or gaussian-rbf
    char solnscheme[80] = "lbfgs";  // `newton', ``descent', 'lbfgs'
    char solnschemef[80];
    char priorweightf[80];
    ConvertToFortran(solnschemef, sizeof solnschemef, solnscheme);
    ConvertToFortran(priorweightf, sizeof priorweightf, priorweight);
            
    int maxit = 100;
    double eps = 1E-6, charlengthscale = 1, dettol = 1E-15;
    int printflag = 0;
    int onlyPhi = 1;
    int ierrorflag;
    
    double nodesXYZCurr[numNodes * ndim];// Temporal array with the coordinates
    for (int j = 0; j < numNodes; j++) 
    {
        int ind = nodesMM[j];
        vector<double> xyz22;
        if (flagCurrCon) 
        {
            // current config
            xyz22 = (nodes[ind]->getXYZCurr());
        } 
        else 
        {
            // initial config
            xyz22 = (nodes[ind]->getXYZ());
        }
        for (int i=0; i< ndim; i++)
        {
            nodesXYZCurr[j * ndim+i] = xyz22[i];
        }
    }
    if (flagSupport) 
    {
        supportDomain(nodesXYZCurr, ndim, numNodes, rmaxCurr, charlengthscale, 0);
        for (int j = 0; j < nodesMM.size(); j++) 
        {
            int nod = nodesMM[j];
            nodes[nod]->setRmax(rmaxCurr[j]);
        }
    } 
    else 
    {
        charlengthscale = 0;
        for (int j = 0; j < nodesMM.size(); j++) 
        {
            int nod = nodesMM[j];
            rmaxCurr[j] = nodes[nod]->getRmax();
            charlengthscale += rmaxCurr[j];
            if (rmaxCurr[j] <= 0) {
                ERROR("The support domain of the nodes has not been properly defined");
                exit(-1);
            }
        }
        charlengthscale/= double(numNodes);
    }

    for (int j = 0; j < nodesMM.size(); j++) 
    {
        int ind = nodesMM[j];
        vector<double> xyz;
        if (flagCurrCon) 
        {
            xyz = nodes[ind]->getXYZCurr();
        } 
        else 
        {
            xyz = nodes[ind]->getXYZ();
        }
        vector<int> neighboursLocal(0);  // neighbours indexes in the MM list of nodes        
        findNeighbours(xyz, nodes[ind]->getMe(), nodesXYZCurr, ndim, numNodes, neighboursLocal, rmaxCurr);
        vector<int> neighbours(neighboursLocal.size(),0);
        for (int j = 0; j < neighboursLocal.size(); j++) 
        {
            neighbours[j] = nodesMM[neighboursLocal[j]];
        }
        int numNeighbours = neighbours.size();
        double nodesXYZNeighbours[numNeighbours * ndim]; //nodal coordinates of the n neighbors
        vector<double> rmaxLocal(numNeighbours, 0);
        for (int j = 0; j < numNeighbours; j++) 
        {
            const vector<double>& xzyPoint= nodes[neighbours[j]]->getXYZ();
            for (int i=0; i< ndim; i++)
            {
                nodesXYZNeighbours[j+i*numNeighbours] = xzyPoint[i];
            }
            rmaxLocal[j] = rmaxCurr[neighboursLocal[j]];
        }
        
        double D[ndim*ndim*numNeighbours]; //anisotropic metric = I when the isotropic case is used
        fill(D,D+ndim*ndim*numNeighbours,0.);
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< numNeighbours; j++)
            {
                D[i+ndim*i+ndim*ndim*j] = 1.;
            }
        }
        
        double phi[numNeighbours];     // basis functions
        double dphi[ndim * numNeighbours];   // first derivatives of basis functions
        double ddphi[ndim * ndim * numNeighbours];// second derivatives of basis functions
        fill(phi,phi+numNeighbours,0.);
        fill(dphi,dphi+ndim*numNeighbours,0.);
        fill(ddphi,ddphi+ndim*ndim*numNeighbours,0.);
        __maxent_MOD_drivermaxent(&numNeighbours, &ndim, &solnschemef[0], &priorweightf[0], nodesXYZNeighbours, &xyz[0],
                                  &rmaxLocal[0], D, &maxit,
                                  &eps, &printflag, &ierrorflag, &charlengthscale, &dettol, phi, dphi, ddphi);
        for (int i = 0; i < numNeighbours; i++) 
        {
            if (!std::isfinite(phi[i])) {
                WARNING("Cannot get shape function at node %i  numNeighbours=%d", nodes[ind]->getMe(),numNeighbours);
                printVector(xyz,"coordinates");
                printVector(neighbours,"neighboursTRUE");
                WARNING("approximate by fictitious dofs at this point");
                for (int j = 0; j < numNeighbours; j++)
                {
                    if (neighbours[j] == nodes[ind]->getMe())
                    {
                        phi[j]=1;
                    }
                    else
                    {
                        phi[j]=0;
                    }
                }
                break;
            }
        }
        // Be careful with the format data of phi, dphi and ddphi. What I am going to do in order to avoid any transformation, is that dphi has ndim files with numNeighbours columns and ddphi has ndim rows, ndims columns and for each node who is the last one. Careful with the way of accessing to the data in those arrays.
        nodes[ind]->setNeighbours(neighbours, numNeighbours);
        nodes[ind]->setShapeFunctions(phi, dphi, ddphi, ndim, numNeighbours);
    }
};



#else
extern void maxEntSFatGPs(vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, vector<int> &nodesMM, bool flagSupport,
              bool flagCurrCon)
{
  ERROR("the code with be compiled with MM");
  exit(-1);
}

extern void maxEntSFatNod(vector<classNodes *> &nodes, vector<classNodes *> &nodes2, int ndim, vector<int> &nodesMM,
                          bool flagSupport, bool flagCurrCon)
{
  ERROR("the code with be compiled with MM");
  exit(-1);
}

#endif //MM
