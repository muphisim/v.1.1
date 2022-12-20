//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file boundaryConditions.cpp
  \brief This file contains all functions related to the application of the boundary conditions: Dirichlet and Neumann
*/
#include "boundaryConditions.h"
#include "GPsDistribution.h"
#include "classNeumannBCs.h"
#include "constitutiveList.h"

static void BCsNodalForces(vector<classNodes *> &nodes, vector<classElements *> &elements,
                           vector<classNeumann *> &NeumannNodes, int &ndim, vector<double> &force, double timeRun,
                           double dt);
static void BCsTractionRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classTractionRamp *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt);

static void BCsTractionInst(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classTractionInst *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt);

static void BCsStochasticPressureRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStochasticPressureRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt);

static void BCsPressureRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classPressureRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt,
                            bool tangentEstimation, map<pair<int,int>,double>& stiffBC);

static void BCsHertzianRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classHertzianRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt, 
                             bool tangentEstimation, map<pair<int,int>,double>& stiffBC);

static void BCsStochasticHertzianRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStochasticHertzianRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt);

static void BCsPressureInst(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classPressureInst *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt);

static void BCsGravity(vector<classNodes *> &nodes, vector<classElements *> &elements, int &ndim, vector<double> &force,
                       classGravity *&BCsGrav);
                       
static void BCsVolHeatFlux(vector<classGPs *> &GPs,  vector<classNodes *> &nodes, classVolHeatFlux *&BCsVolHeatFlux, int &ndim,
                       vector<double>& force, double timeRun, double dt);
                       
static void BCsVolHeatFluxGaussian(vector<classGPs *> &GPs,  vector<classNodes *> &nodes, classVolHeatFluxGaussian *&BCsVolHeatFlux, int &ndim,
                       vector<double>& force, double timeRun, double dt);
                       
static void BCsCurrentInst(vector<classGPs *> &GPs, vector<classNodes *> &nodes, classCurrentInst *&BCsCurrent, int &ndim,
                        vector<double>& force, double timeRun, double dt);
                        
static void BCsHeatFluxInst(const vector<double> &uK, vector<double>& force,  vector<classNodes *> &nodes,
                            vector<classElements *> &elements, classSurfaceHeatFluxInst *&BCsHeat, int &ndim, double timeRun,
                            double dt) ;
static void BCsConvection(const vector<double> &uK, vector<double>& force, vector<classNodes *> &nodes,
                          vector<classElements *> &elements, classConvection *&BCsConvection, int &ndim, double timeRun,
                          double dt,  bool tangentEstimation, map<pair<int,int>,double>& stiffBC);
static void BCsRadiation(const vector<double> &uK, vector<double>& force, vector<classNodes *> &nodes,
                         vector<classElements *> &elements, classRadiation *&BCsRadiation, int &ndim, double timeRun,
                         double dt,  bool tangentEstimation, map<pair<int,int>,double>& stiffBC);

/*! \brief This function manages all possible Neumann boundary conditions
  @param[in] uK Array with the displacements
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] NeumannBCs Array with all Neumann BCs (all types)
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
extern void NeumannBCManagement(const vector<double> &uK, vector<classGPs *> &GPs, vector<classNodes *> &nodes, vector<classElements *> &elements,
                                vector<classNeumannBCs *> &NeumannBCs, int &ndim, vector<double> &force, double timeRun,
                                double dt, bool tangentEstimation, map<pair<int,int>,double>& stiffBC) {
    fill(force.begin(), force.end(), 0); //reset to zero at each time step

    for (int i = 0; i < NeumannBCs.size(); i++) {
        if (NeumannBCs[i]->getType() == "FORCE") {
            classNodalForces *noda = static_cast<classNodalForces *>(NeumannBCs[i]);
            vector<classNeumann *> lal = noda->getListOfNodes();
            BCsNodalForces(nodes, elements, lal, ndim, force, timeRun, dt);
        } else if (NeumannBCs[i]->getType() == "TRACTION RAMP") {
            classTractionRamp *nodas = static_cast<classTractionRamp *>(NeumannBCs[i]);
            BCsTractionRamp(uK, nodes, elements, nodas, ndim, force, timeRun, dt);
        } else if (NeumannBCs[i]->getType() == "TRACTION INST") {
            classTractionInst *nodas = static_cast<classTractionInst *>(NeumannBCs[i]);
            BCsTractionInst(uK, nodes, elements, nodas, ndim, force, timeRun, dt);
        } else if (NeumannBCs[i]->getType() == "PRESSURE RAMP") {
            classPressureRamp *nodas = static_cast<classPressureRamp *>(NeumannBCs[i]);
            BCsPressureRamp(uK, nodes, elements, nodas, ndim, force, timeRun, dt, tangentEstimation,stiffBC);
        } else if (NeumannBCs[i]->getType() == "PRESSURE INST") {
            classPressureInst *nodas = static_cast<classPressureInst *>(NeumannBCs[i]);
            BCsPressureInst(uK, nodes, elements, nodas, ndim, force, timeRun, dt);
        } else if (NeumannBCs[i]->getType() == "GRAVITY") {
            classGravity *nodas = static_cast<classGravity *>(NeumannBCs[i]);
            BCsGravity(nodes, elements, ndim, force, nodas);
        } else if (NeumannBCs[i]->getType() == "HERTZIAN RAMP") {
            classHertzianRamp *nodas = static_cast<classHertzianRamp *>(NeumannBCs[i]);
            BCsHertzianRamp(uK, nodes, elements, nodas, ndim, force, timeRun, dt, tangentEstimation, stiffBC);
        }
        else if (NeumannBCs[i]->getType() == "STOCHASTIC HERTZIAN RAMP") {
            classStochasticHertzianRamp *nodas = static_cast<classStochasticHertzianRamp *>(NeumannBCs[i]);
            BCsStochasticHertzianRamp(uK, nodes, elements, nodas, ndim, force, timeRun, dt);
        } else if (NeumannBCs[i]->getType() == "STOCHASTIC PRESSURE RAMP") {
            classStochasticPressureRamp *nodas = static_cast<classStochasticPressureRamp *>(NeumannBCs[i]);
            BCsStochasticPressureRamp(uK, nodes, elements, nodas, ndim, force, timeRun, dt);
        }
        else if (NeumannBCs[i]->getType() == "VOLUMETRIC HEAT FLUX"){
            classVolHeatFlux *nodas = static_cast<classVolHeatFlux *>(NeumannBCs[i]);
            BCsVolHeatFlux(GPs, nodes, nodas, ndim, force, timeRun, dt);
        }
        else if (NeumannBCs[i]->getType() == "GAUSSIAN VOLUMETRIC HEAT FLUX"){
            classVolHeatFluxGaussian *nodas = static_cast<classVolHeatFluxGaussian *>(NeumannBCs[i]);
            BCsVolHeatFluxGaussian(GPs, nodes, nodas, ndim, force, timeRun, dt);
        }
        else if (NeumannBCs[i]->getType() == "CURRENT INST") {
            classCurrentInst *nodas = static_cast<classCurrentInst *>(NeumannBCs[i]);
            BCsCurrentInst(GPs, nodes, nodas, ndim, force, timeRun, dt);
        }
        else if (NeumannBCs[i]->getType() == "HEAT INST") {
            classSurfaceHeatFluxInst *nodas = static_cast<classSurfaceHeatFluxInst *>(NeumannBCs[i]);
            BCsHeatFluxInst(uK, force, nodes, elements, nodas, ndim, timeRun, dt);
        }
        else if (NeumannBCs[i]->getType() == "CONVECTION") {
            classConvection *nodas = static_cast<classConvection *>(NeumannBCs[i]);
            BCsConvection(uK, force, nodes, elements, nodas, ndim, timeRun, dt, tangentEstimation, stiffBC);
        } 
        else if (NeumannBCs[i]->getType() == "RADIATION") {
            classRadiation *nodas = static_cast<classRadiation *>(NeumannBCs[i]);
            BCsRadiation(uK, force, nodes, elements, nodas, ndim, timeRun, dt, tangentEstimation, stiffBC);
        }
        
        
    }
}

/*! \brief This function deletes all possible Neumann boundary conditions
  @param[in] NeumannBCs Array with all Neumann BCs (all types)
*/
extern void deleteNeumannBCs(vector<classNeumannBCs *> &NeumannBCs) {
    for (int i = 0; i < NeumannBCs.size(); i++) {
        if (NeumannBCs[i]->getType() == "FORCE") {
            classNodalForces *noda = static_cast<classNodalForces *>(NeumannBCs[i]);
            vector<classNeumann *> lal = noda->getListOfNodes();
            for (int j = 0; j < lal.size(); j++) {
                delete lal[j];
            }
        } else if (NeumannBCs[i]->getType() == "TRACTION RAMP") {
            classTractionRamp *nodas = static_cast<classTractionRamp *>(NeumannBCs[i]);
            vector<classElements *> lal = nodas->getListOfElements();
            for (int j = 0; j < lal.size(); j++) {
                delete lal[j];
            }
        } else if (NeumannBCs[i]->getType() == "PRESSURE RAMP") {
            classPressureRamp *nodas = static_cast<classPressureRamp *>(NeumannBCs[i]);
            vector<classElements *> lal = nodas->getListOfElements();
            for (int j = 0; j < lal.size(); j++) {
                delete lal[j];
            }
        } else if (NeumannBCs[i]->getType() == "PRESSURE INST") {
            classPressureInst *nodas = static_cast<classPressureInst *>(NeumannBCs[i]);
            vector<classElements *> lal = nodas->getListOfElements();
            for (int j = 0; j < lal.size(); j++) {
                delete lal[j];
            }
        }  // No sub Elements to delete for extradof flux
        delete NeumannBCs[i];
    }
    NeumannBCs.erase(NeumannBCs.begin(), NeumannBCs.end());
}


/*! \brief This function finds all dof that are active during the simulation
  @param[in] DirichletNodes Dirichlet boundary nodes
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of degrees of freedom
  @param[out] activeDof Array with the labels of the active dof
*/
extern void findActiveDof(vector<classDirichlet *> DirichletNodes, int ndim, int numDof, vector<int> &activeDof) {
    activeSearching(DirichletNodes, ndim, numDof, activeDof);
}

/*! \brief This function finds all dof that are active during the simulation
  @param[in] DirichletNodes Dirichlet boundary nodes
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of degrees of freedom
  @param[out] activeDof Array with the labels of the active dof
*/
extern void activeSearching(vector<classDirichlet *> DirichletNodes, int ndim, int numDof, vector<int> &activeDof) {
    vector<int> givenDof;
    int numGivenDof = DirichletNodes.size();
    for (int i = 0; i < numGivenDof; i++) {
        givenDof.push_back(DirichletNodes[i]->getMe());
    }

    for (int i = 0; i < numDof; i++) {
        vector<int>::iterator it;
        it = find(givenDof.begin(), givenDof.end(), i);
        if (it == givenDof.end()) {
            activeDof.push_back(i);
        }
    }
};

/*! \brief This function applies the traction (Neumann boundary condition)
  @param[in] uK Array with the displacements
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[inout] BCsTraction Neumann BCs traction
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsTractionInst(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classTractionInst *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt) {
    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsTraction->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    std::vector<double> Traction(3,0.);
    BCsTraction->getTraction(timeRun, dt,Traction[0],Traction[1],Traction[2]);
    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        for (int a = 0; a < neighbours.size(); a++) {
            for (int i = 0; i < ndim; i++) {
                int indi11 = neighbours[a] + nodes.size() * i;
                force[indi11] += weight * std::abs(J) * phi[a] * Traction[i];
            }
        }
        delete tempGP;
    }
}


/*! \brief This function applies the traction (Neumann boundary condition)
  @param[in] uK Array with the displacements
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[inout] BCsTraction Neumann BCs traction
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsTractionRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classTractionRamp *&BCsTraction, int &ndim, vector<double> &force, double timeRun, double dt) {
    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsTraction->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    std::vector<double> Traction(3,0.);
    BCsTraction->getTraction(timeRun, dt,Traction[0],Traction[1],Traction[2]);
    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        for (int a = 0; a < neighbours.size(); a++) {
            for (int i = 0; i < ndim; i++) {
                int indi11 = neighbours[a] + nodes.size() * i;
                force[indi11] += weight * std::abs(J) * phi[a] * Traction[i];
            }
        }
        delete tempGP;
    }
}

/*! \brief This function applies the nodal forces (Neumann boundary condition)
  @param[in] uK Array with the displacements
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[inout] BCsPress Neumann BCs pressure
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsStochasticPressureRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStochasticPressureRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsPress->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    vector<double> P = BCsPress->getPressure(timeRun, dt);
    vector<constitutiveModels *> elementsCons = elements[0]->getConsModel();
    int nbr_dof = 0 ;
    nbr_dof = (elementsCons[0]->getNbrDofConsMod());
    int nstoch = 0;
    nstoch = nbr_dof/ndim ;
    vector<double> C_bis(nstoch*nstoch*nstoch,0);
        const classStochasticHyperElasticStVenantKirchhoff* firstLaw = dynamic_cast<const classStochasticHyperElasticStVenantKirchhoff*>(elementsCons[0]);
    if (firstLaw == NULL)
    {
        ERROR("classStochasticHyperElasticStVenantKirchhoff must be considered");
        exit(-1);
    }
    C_bis = firstLaw->getC_Stoch();

    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        vector<double> n0 = elems[tempGP->getCell()]->getn0(nodes, ndim);// Be carefull with this call
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim * nstoch, 0);
        vector<double> FFinv(ndim * ndim * nstoch, 0);
        vector<double> FFinvT(ndim * ndim * nstoch, 0);
        vector<double> n(ndim * nstoch, 0);
        vector<double> n0_stoch(ndim * nstoch, 0);
        for (int i=0; i<ndim; i++){
            n0_stoch[i] = n0[i];
        }

        tempGPF->getF(uK, nodes.size(), ndim, FF);
        vector <double> Jacobian(nstoch,0);
        vector <double> Jacobian_P(nstoch,0);
        DeterminantStochasticTensor(FF, ndim, nstoch, Jacobian, C_bis);
        InverseTensorStochastic(FF, ndim, nstoch, C_bis, FFinv);
/*        vector<double> FF_simplified(ndim*ndim, 0);
        vector<double> FF_inv_simplified(ndim*ndim, 0);
        FF_simplified = {FF.begin(), FF.begin() + ndim*ndim};
        double Jacobian = determinantTensor(FF_simplified, ndim);
        inverse(FF_simplified, ndim, FF_inv_simplified);
        for (int i=0;i<ndim*ndim;i++){
            FFinv[i] = FF_inv_simplified[i];
        }*/
        transposeStochastic(FFinv, ndim, nstoch, FFinvT);
        multTensorVectorStochastic(FFinvT, ndim, nstoch, n0_stoch, C_bis, n);
        multiRandomVariableRandomVariableStochastic(P, Jacobian, C_bis, Jacobian_P);
        for (int a = 0; a < neighbours.size(); a++) {
            for (int j = 0; j < nstoch; j++) {
                for (int i = 0; i < ndim; i++) {
                    int indi11 = neighbours[a] + nodes.size() * i +  j * ndim * nodes.size();
                    force[indi11] += weight * std::abs(J) * phi[a] * Jacobian_P[j] * n[i];
                }
            }
        }
        delete tempGPF;
        delete tempGP;
    }

}



/*! \brief This function applies the nodal forces (Neumann boundary condition)
  @param[in] uK Array with the displacements
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[inout] BCsPress Neumann BCs pressure
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsPressureRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classPressureRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt,
                            bool tangentEstimation, map<pair<int,int>,double>& stiffBC) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsPress->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    double P = BCsPress->getPressure(timeRun, dt);

    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        vector<double> n0 = elems[tempGP->getCell()]->getn0(nodes, ndim);// Be carefull with this call
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        const vector<double>& phi = tempGPF->getPhi();
        const vector<double>& dphi = tempGPF->getDphi();
        const vector<int>& neighbours = tempGPF->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim, 0);
        vector<double> FFinv(ndim * ndim, 0);
        vector<double> FFinvT(ndim * ndim, 0);
        vector<double> n(ndim, 0);
        tempGPF->getF(uK, nodes.size(), ndim, FF);
        double Jacobian = determinantTensor(FF, ndim); // Volume
        inverse(FF, ndim, FFinv);
        transpose(FFinv, ndim, FFinvT);
        multTensorVector(FFinvT, ndim, ndim, n0, ndim, n);
        
        double DnJacobianDF[3][3][3];
        if (tangentEstimation)
        {
            for (int j=0; j< ndim; j++)
            {
                for (int k=0; k< ndim; k++)
                {
                    for (int l=0; l< ndim; l++)
                    {
                        DnJacobianDF[j][k][l] = 0.;
                        for (int i=0; i< ndim; i++)
                        {
                            DnJacobianDF[j][k][l] += n0[i]*(Jacobian*FFinv[i*ndim+j]*FFinv[l*ndim+k] -Jacobian*FFinv[i*ndim+k]*FFinv[l*ndim+j]);
                        }
                    }
                }
            }
        }
        for (int a = 0; a < neighbours.size(); a++) {
            for (int i = 0; i < ndim; i++) {
                int indi11 = neighbours[a] + nodes.size() * i;
                force[indi11] += weight * std::abs(J) * phi[a] * Jacobian * P * n[i];
                if (tangentEstimation)
                {
                    for (int b = 0; b < neighbours.size(); b++) 
                    {
                        for (int j = 0; j < ndim; j++) 
                        { 
                            int indi12 = neighbours[b] + nodes.size() * j;
                            for (int kk=0; kk< ndim; kk++)
                            {
                                double DforceDu = weight * std::abs(J) * phi[a]  * P * DnJacobianDF[i][j][kk]*dphi[b*ndim + kk];
                                if (stiffBC.find(pair<int,int>(indi11,indi12)) == stiffBC.end())
                                {
                                    stiffBC[pair<int,int>(indi11,indi12)] = DforceDu;
                                }
                                else
                                {
                                    stiffBC[pair<int,int>(indi11,indi12)] += DforceDu;
                                }
                            }
                        }
                    }
                    
                }
            }
        }
        delete tempGPF;
        delete tempGP;
    }

}

static void BCsPressureInst(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classPressureInst *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsPress->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    double P = BCsPress->getPressure(timeRun, dt);

    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        vector<double> n0 = elems[tempGP->getCell()]->getn0(nodes, ndim);// Be carefull with this call
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim, 0);
        vector<double> FFinv(ndim * ndim, 0);
        vector<double> FFinvT(ndim * ndim, 0);
        vector<double> n(ndim, 0);
        tempGPF->getF(uK, nodes.size(), ndim, FF);
        double Jacobian = determinantTensor(FF, ndim); // Volume
        inverse(FF, ndim, FFinv);
        transpose(FFinv, ndim, FFinvT);
        multTensorVector(FFinvT, ndim, ndim, n0, ndim, n);
        for (int a = 0; a < neighbours.size(); a++) {
            for (int i = 0; i < ndim; i++) {
                int indi11 = neighbours[a] + nodes.size() * i;
                force[indi11] += weight * std::abs(J) * phi[a] * Jacobian * P * n[i];
            }
        }
        delete tempGPF;
        delete tempGP;
    }

}



/*! \brief This function applies the nodal forces (Neumann boundary condition)
  @param[in] uK Array with the displacements
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[inout] BCsPress Neumann BCs pressure
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsHertzianRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classHertzianRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt,
                             bool tangentEstimation, map<pair<int,int>,double>& stiffBC) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsPress->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    double P = BCsPress->getNormalPressure(timeRun, dt);
    double Q = BCsPress->getTangentialPressure(timeRun, dt);
    double rad = BCsPress->getRadius();
    double c = BCsPress->getRadiusStick();
    vector<double> center = BCsPress->getCenter();

    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        bool isInContact = false;
        bool isInStick = false;

        classGPs *tempGP = *ite;
        vector<double> n0 = elems[tempGP->getCell()]->getn0(nodes, ndim);// Be carefull with this call
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        vector<double> phi = tempGP->getPhi();
        vector<double> dphi = tempGP->getDphi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        isInCircle(isInContact, n0, center, Point, rad);
        isInCircle(isInStick, n0, center, Point, c);
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim, 0);
        vector<double> FFinv(ndim * ndim, 0);
        vector<double> FFinvT(ndim * ndim, 0);
        vector<double> n(ndim, 0);
        vector<double> n_test(ndim, 0);
        vector<double> n_tang(ndim, 0);
        n_test[0] = 1;
        double r;
        r = distance(center, Point);

        tempGPF->getF(uK, nodes.size(), ndim, FF);
        double Jacobian = determinantTensor(FF, ndim); // Volume
        inverse(FF, ndim, FFinv);
        transpose(FFinv, ndim, FFinvT);
        multTensorVector(FFinvT, ndim, ndim, n0, ndim, n);
        multTensorVector(FFinvT, ndim, ndim, n_test, ndim, n_tang);
        if (isInContact) {
            for (int a = 0; a < neighbours.size(); a++) {
                for (int i = 0; i < ndim; i++) {
                    int indi11 = neighbours[a] + nodes.size() * i;
                    force[indi11] += weight * std::abs(J) * phi[a] * Jacobian * P * n[i] * sqrt(1 - pow(r / rad, 2));
                    force[indi11] += weight * std::abs(J) * phi[a] * Jacobian * Q * n_tang[i] * sqrt(1 - pow(r / rad, 2));
                    for (int b = 0; b < neighbours.size(); b++) {
                        for (int k = 0; k < ndim; k++) {
                            double value = 0;
                            int indi22 = neighbours[b] + nodes.size() * k;
                            for (int J = 0; J < ndim; J++) {
                                for (int L = 0; L < ndim; L++) { // the 0 deactivates the geometrical contribution of the ext force
                                    value -= 0 * weight * std::abs(J) * phi[a] * Jacobian * P * sqrt(1 - pow(r / rad, 2)) *
                                             dphi[b * ndim + L] * (FFinvT[i * ndim + J] * n0[J] * FFinv[L * ndim + k] -
                                                                   FFinv[J * ndim + k] * FFinv[L * ndim + i] *
                                                                     n0[J]);
                                    value -= 0 * weight * std::abs(J) * phi[a] * Jacobian * Q * sqrt(1 - pow(r / rad, 2)) *
                                             dphi[b * ndim + L] * (FFinvT[i * ndim + J] * n_test[J] * FFinv[L * ndim + k] -
                                                                   FFinv[J * ndim + k] * FFinv[L * ndim + i] * n_test[J]);
                                }
                            }
                            if (stiffBC.find(pair<int,int>(indi11,indi22)) == stiffBC.end())
                            {
                                stiffBC[pair<int,int>(indi11,indi22)] = value;
                            }
                            else
                            {
                                stiffBC[pair<int,int>(indi11,indi22)] += value;
                            }
                        }
                    }
                }
            }
            if (isInStick) {
                for (int a = 0; a < neighbours.size(); a++) {
                    for (int i = 0; i < ndim; i++) {
                        int indi11 = neighbours[a] + nodes.size() * i;
                        force[indi11] += -weight * std::abs(J) * phi[a] * Jacobian * n_tang[i] * Q * (c / rad) *
                                         sqrt(1 - pow(r / c, 2));
                        for (int b = 0; b < neighbours.size(); b++) {
                            for (int k = 0; k < ndim; k++) {
                                double value = 0;
                                int indi22 = neighbours[b] + nodes.size() * k;
                                for (int J = 0; J < ndim; J++) {
                                    for (int L = 0; L < ndim; L++) {
                                        value -= -0 * weight * std::abs(J) * phi[a] * Jacobian * Q * (c / rad) *
                                                 sqrt(1 - pow(r / c, 2)) * dphi[b * ndim + L] *
                                                 (FFinvT[i * ndim + J] * n_test[J] * FFinv[L * ndim + k] -
                                                  FFinv[J * ndim + k] * FFinv[L * ndim + i] *
                                                  n_test[J]);
                                    }
                                }
                                if (stiffBC.find(pair<int,int>(indi11,indi22)) == stiffBC.end())
                                {
                                    stiffBC[pair<int,int>(indi11,indi22)] = value;
                                }
                                else
                                {
                                    stiffBC[pair<int,int>(indi11,indi22)] += value;
                                }
                            }
                        }

                    }
                }
            }
        }
        delete tempGPF;
        tempGPF = NULL;
        delete tempGP;
        tempGPF = NULL;
    }
}

/*! \brief This function applies the nodal forces (Neumann boundary condition)
  @param[in] uK Array with the displacements
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[inout] BCsPress Neumann BCs pressure
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsStochasticHertzianRamp(const vector<double> &uK, vector<classNodes *> &nodes, vector<classElements *> &elements,
                            classStochasticHertzianRamp *&BCsPress, int &ndim, vector<double> &force, double timeRun, double dt) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsPress->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    double P = BCsPress->getNormalPressure(timeRun, dt);
    vector<double> Q = BCsPress->getTangentialPressure(timeRun, dt);
    vector<double> mu = BCsPress->getFriction();
    double rad = BCsPress->getRadius();
    vector<double> c = BCsPress->getRadiusStick();
    vector<double> center = BCsPress->getCenter();
    vector<constitutiveModels *> elementsCons = elements[0]->getConsModel();
    int nbr_dof  = (elementsCons[0]->getNbrDofConsMod());
    int nstoch  = nbr_dof/ndim;
    vector<double> C_bis(nstoch*nstoch*nstoch,0);
    const classStochasticHyperElasticStVenantKirchhoff* firstLaw = dynamic_cast<const classStochasticHyperElasticStVenantKirchhoff*>(elementsCons[0]);
    if (firstLaw == NULL)
    {
        ERROR("classStochasticHyperElasticStVenantKirchhoff must be considered");
        exit(-1);
    }
    C_bis = firstLaw->getC_Stoch();
    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        bool isInContact = false;
        classGPs *tempGP = *ite;
        vector<double> no = elems[tempGP->getCell()]->getn0(nodes, ndim);// Be carefull with this call
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        isInCircle(isInContact,no,center,Point,rad);
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim* nstoch, 0);
        vector<double> FFinv(ndim * ndim * nstoch, 0);
        vector<double> FFinvT(ndim * ndim * nstoch, 0);
        vector<double> n(ndim *nstoch, 0);
        vector<double> n_test(ndim * nstoch, 0);
        vector<double> n_tang(ndim * nstoch, 0);
        vector<double> nobis(ndim*nstoch,0);
        vector<double> Jacobian(nstoch, 0);
        vector<double> relative_distance(nstoch, 0);
        vector<double> sigmoid_value(nstoch, 0);
        vector<double> sigmoid_complementary_value(nstoch, 0);
        relative_distance = c;

        for(int i=0; i<ndim;i++){
            nobis[i] = no[i];
        }

        n_test[0] = 1;
        double r = 0;
        r = distance (center,Point);
 
        tempGPF->getF(uK, nodes.size(), ndim, FF);
        DeterminantStochasticTensor(FF, ndim, nstoch, Jacobian,C_bis);
        InverseTensorStochastic(FF, ndim, nstoch, C_bis, FFinv);
        transposeStochastic(FFinv, ndim, nstoch, FFinvT);

        multTensorVectorStochastic(FFinvT, ndim, nstoch, nobis, C_bis, n);


        multTensorVectorStochastic(FFinvT, ndim, nstoch, n_test, C_bis, n_tang);

        double alpha = 1; //Numerical parameter
        vector<double> dist(nstoch,0);
        vector<double> ratio(nstoch,0);
        dist[0] = r*alpha;

        vector<double> approximation(nstoch,0);
        DivisionRandomVariableRandomVariableStochastic(dist, c, C_bis,
                                                       ratio);

        piecewisefit(r, c, C_bis, approximation);
        vector<double>r_over_c(nstoch,0);
        r_over_c = ratio;
        vector<double> square_ratio(nstoch,0);
        multiRandomVariableRandomVariableStochastic(ratio, ratio, C_bis, square_ratio);

        linearCombinationVariableStochastic(square_ratio,-1,1);

        //weight_2 = ratio;
        vector<double> c_over_rad(nstoch,0);
        c_over_rad = c;
        linearCombinationVariableStochastic(c_over_rad,1/rad,0);
        vector<double>sqrt_1minusx(nstoch,0);
        vector<double>slip_tang(nstoch,0);
        vector<double>stick_tang(nstoch,0);
        slip_tang = Q;
        double coeff_1 = 0;
        if((r/rad)<1) {
            coeff_1 = sqrt(1 - pow( r / rad,2));
        }
        linearCombinationVariableStochastic(slip_tang,coeff_1,0);
        vector<double>sqrt_r_over_c(nstoch,0);
        if(r<1.65){
            squareRootIntegral(square_ratio, sqrt_r_over_c, C_bis);
        }
        if((r>1.65)&&(r<2.1634)){
            InterpolationSqrt(r_over_c,sqrt_r_over_c,C_bis);
        }

        vector<double> prod_tang(nstoch,0);
        multiRandomVariableRandomVariableStochastic(sqrt_r_over_c, c_over_rad, C_bis, prod_tang);
        vector<double> tang_force_stick (ndim*nstoch,0);
        vector<double> Q_n (ndim*nstoch,0);
        vector<double> approximation_force(ndim*nstoch, 0);
        multiSVectorRandomVariableStochastic(n_tang, ndim, Q, nstoch, C_bis,
                                             Q_n);
        multiSVectorRandomVariableStochastic(Q_n, ndim, prod_tang, nstoch, C_bis,
                                             tang_force_stick);
        multiSVectorRandomVariableStochastic(Q_n, ndim, approximation, nstoch, C_bis,
                                             approximation_force);
        vector<double> tang_force_slip (ndim*nstoch,0);
        multiSVectorRandomVariableStochastic(n_tang, ndim, slip_tang, nstoch, C_bis,
                                             tang_force_slip);
        vector<double> J_n(ndim*nstoch,0);
        vector<double> J_tang_slip(ndim*nstoch,0);
        vector<double> J_tang_stick(ndim*nstoch,0);
        vector<double> J_tang_stick_weighted(ndim * nstoch, 0);
        vector<double> J_approximation_force(ndim * nstoch, 0);
        multiSVectorRandomVariableStochastic(n, ndim, Jacobian, nstoch, C_bis,
                                             J_n);
        multiSVectorRandomVariableStochastic(tang_force_slip, ndim, Jacobian, nstoch, C_bis,
                                             J_tang_slip);
        multiSVectorRandomVariableStochastic(tang_force_stick, ndim, Jacobian, nstoch, C_bis,
                                             J_tang_stick);
        multiSVectorRandomVariableStochastic(approximation_force, ndim, Jacobian, nstoch, C_bis,
                                            J_approximation_force);
        //multiSVectorRandomVariableStochastic(tang_force_stick, ndim, Jacobian, nstoch, C_bis,
        // J_tang_stick);
        //multiSVectorRandomVariableStochastic(J_tang_stick, ndim, sigmoid_complementary_value, nstoch, C_bis,
        // J_tang_stick_weighted);

        if(isInContact) {
            for (int a = 0; a < neighbours.size(); a++) {
                for (int j = 0; j < nstoch; j++) {
                    for (int i = 0; i < ndim; i++) {
                        int indi11 = neighbours[a] + nodes.size() * i +  j * ndim * nodes.size();
                        //force[indi11] +=
                                //weight * std::abs(J) * phi[a] * Jacobian * P * n[j*ndim + i] * sqrt(1 - pow(r / rad,2));
                        //force[indi11] +=
                                //weight * std::abs(J) * phi[a] * Jacobian * ( tang_force_slip[j*ndim + i] - tang_force_stick[j*ndim + i]) ;
                        force[indi11] += weight * std::abs(J) * phi[a] * P * J_n[j*ndim + i] * sqrt(1 - pow(r / rad,2));
                        force[indi11] += weight * std::abs(J) * phi[a] * (J_tang_slip[j*ndim + i] - J_tang_stick[j*ndim + i]);                      
                    }
                }
            }
        }
         delete tempGPF;
         delete tempGP;
    }

}




/*! \brief This function applies a volumetric heat flux on Extra dof (Neumann boundary condition) 
  @param[in] GPs all GPs
  @param[in] nodes Array with all nodes in the domain
  @param[in] BCsCurrent Neumann BCs Volumetric heat flux
  @param[in] ndim Dimension of the domain
  @param[out] force external force
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsVolHeatFluxGaussian(vector<classGPs *> &GPs, vector<classNodes *> &nodes, classVolHeatFluxGaussian *&BCsVolHeatFlux, int &ndim,
                           vector<double>& force, double timeRun, double dt) {
    int numNodes = nodes.size();
    int fieldIndex = BCsVolHeatFlux->getFieldIndex();
    
    double t0 = BCsVolHeatFlux->get_t0();
    double tf = BCsVolHeatFlux->get_tf();
    if ((timeRun < t0) or (timeRun > tf))
        return;
    
    for (int igp = 0; igp < GPs.size(); igp++)
    {
        classGPs *GaussPoint = GPs[igp];
        if (GaussPoint->getActivate() ==0)
            continue;
        
        vector<double> XYZ =  GaussPoint->getXYZ();
        for (int i=XYZ.size(); i< 3; i++)
        {
            // if 2D is used, the last coordinate = 0
            XYZ.push_back(0);
        }
        double r = BCsVolHeatFlux->getVolHeatFlux(XYZ[0],XYZ[1],XYZ[2],timeRun, dt);
        //
        const vector<int>& neighbours = GaussPoint->getNeighbours(); 
        int numNeighbours = neighbours.size();
        const vector<double>& phi = GaussPoint->getPhi();
        int numMechDofsPerNode = GaussPoint->getCurrentState().getNumMechanicalDofsPerNode();
        double wJ = GaussPoint->getWeightJ();
        for (int a = 0; a < numNeighbours; a++) 
        {
            int id1 = neighbours[a] + (numMechDofsPerNode+ fieldIndex) * numNodes;
            force[id1] += wJ * r * phi[a];
        }    
    }
}

/*! \brief This function applies a volumetric heat flux on Extra dof (Neumann boundary condition) 
  @param[in] GPs all GPs
  @param[in] nodes Array with all nodes in the domain
  @param[in] BCsCurrent Neumann BCs Volumetric heat flux
  @param[in] ndim Dimension of the domain
  @param[out] force external force
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsVolHeatFlux(vector<classGPs *> &GPs, vector<classNodes *> &nodes, classVolHeatFlux *&BCsVolHeatFlux, int &ndim,
                           vector<double>& force, double timeRun, double dt) {

    vector<classElements *>& elems = BCsVolHeatFlux->getListOfElements();
    double r = BCsVolHeatFlux->getVolHeatFlux(timeRun, dt);
    if (fabs(r) < 1e-12) return;
    int numNodes = nodes.size();
    int fieldIndex = BCsVolHeatFlux->getFieldIndex();
    for (vector<classElements *>::iterator ite = elems.begin(); ite != elems.end(); ++ite) 
    {
        // Catching all the GPs refered to tempElem
        classElements *tempElem = *ite;
        const vector<int>& listIndGPs = tempElem->getMyGPs();
        for (int j=0; j< listIndGPs.size(); j++)
        {
            classGPs *GaussPoint = GPs[listIndGPs[j]];
            if (GaussPoint->getActivate() ==0)
                continue;
            const vector<int>& neighbours = GaussPoint->getNeighbours(); 
            int numNeighbours = neighbours.size();
            const vector<double>& phi = GaussPoint->getPhi();
            int numMechDofsPerNode = GaussPoint->getCurrentState().getNumMechanicalDofsPerNode();
            double wJ = GaussPoint->getWeightJ();
            for (int a = 0; a < numNeighbours; a++) 
            {
                int id1 = neighbours[a] + (numMechDofsPerNode+ fieldIndex) * numNodes;
                force[id1] += wJ * r * phi[a];
            }
        }
    }
}


/*! \brief This function applies the current on Extra dof (Neumann boundary condition)
  @param[in] GPs all GPs
  @param[in] nodes Array with all nodes in the domain
  @param[in] BCsCurrent Neumann BCs current
  @param[in] ndim Dimension of the domain
  @param[out] force external force
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
void BCsCurrentInst(vector<classGPs *> &GPs, vector<classNodes *> &nodes, classCurrentInst *&BCsCurrent, int &ndim,
                        vector<double>& force, double timeRun, double dt)

{
    vector<classElements *> elems = BCsCurrent->getListOfElements();
    double r = BCsCurrent->getCurrent(timeRun, dt);
    int fieldIndex = BCsCurrent->getFieldIndex();
    int numNodes = nodes.size();
    for (vector<classElements *>::iterator ite = elems.begin(); ite != elems.end(); ++ite) {
        // Catching all the GPs refered to tempElem
        classElements *tempElem = *ite;
        vector<int> listIndGPs = tempElem->getMyGPs();

        for (int indGPlocal = 0; indGPlocal < listIndGPs.size(); indGPlocal++) 
        {
            classGPs *GaussPoint = GPs[listIndGPs[indGPlocal]];
            if (GaussPoint->getActivate() ==0)
                continue;
                
            const vector<int>& neighbours = GaussPoint->getNeighbours(); 
            int numNeighbours = neighbours.size();
            const vector<double>& phi = GaussPoint->getPhi();
            int numMechDofsPerNode = GaussPoint->getCurrentState().getNumMechanicalDofsPerNode();
            double wJ = GaussPoint->getWeightJ();
            for (int a = 0; a < numNeighbours; a++) 
            {
                int id1 = neighbours[a] + (numMechDofsPerNode+ fieldIndex) * numNodes;
                force[id1] += wJ * r * phi[a];
            }
        }
    }
}

static void BCsHeatFluxInst(const vector<double> &uK, vector<double>& force,  vector<classNodes *> &nodes,
                            vector<classElements *> &elements, classSurfaceHeatFluxInst *&BCsHeat, int &ndim, double timeRun,
                            double dt) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsHeat->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    vector<double> q = BCsHeat->getHeat(timeRun, dt);
    int fieldIndex = BCsHeat->getFieldIndex();
    int numNodes = nodes.size();
    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        vector<double> n0 = elems[tempGP->getCell()]->getn0(nodes, ndim);// Be carefull with this call
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim, 0);
        vector<double> FFinv(ndim * ndim, 0);
        vector<double> FFinvT(ndim * ndim, 0);
        vector<double> n(ndim, 0);
        tempGPF->getF(uK, nodes.size(), ndim, FF);
        double Jacobian = determinantTensor(FF, ndim); // Volume
        inverse(FF, ndim, FFinv);
        transpose(FFinv, ndim, FFinvT);
        multTensorVector(FFinvT, ndim, ndim, n0, ndim, n);
        int numMechDofsPerNode = tempGP->getCurrentState().getNumMechanicalDofsPerNode();
        double wJ = tempGP->getWeightJ();
        for (int a = 0; a < neighbours.size(); a++) {
            int indi11 = neighbours[a] + (numMechDofsPerNode+ fieldIndex) * numNodes;
            double val = 0;
            for (int i = 0; i < ndim; i++) {
                val += wJ * phi[a] * Jacobian * q[i] * n[i];
            }
        }
        delete tempGPF;
        delete tempGP;
    }
};

static void BCsConvection(const vector<double> &uK, vector<double>& force, vector<classNodes *> &nodes,
                          vector<classElements *> &elements, classConvection *&BCsConvection, int &ndim, double timeRun,
                          double dt, bool tangentEstimation, map<pair<int,int>,double>& stiffBC) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsConvection->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    int fieldIndex = BCsConvection->getFieldIndex();
    int numNodes = nodes.size();
    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim, 0);
        tempGPF->getF(uK, nodes.size(), ndim, FF);
        double Jacobian = determinantTensor(FF, ndim); // Volume

        int numMechDofsPerNode = tempGP->getCurrentState().getNumMechanicalDofsPerNode();
        double wJ = tempGP->getWeightJ();
        
        double theta = 0;
        int numNeighbours = neighbours.size();
        for (int a = 0; a < numNeighbours; a++) {
            int ind = neighbours[a]+(numMechDofsPerNode+ fieldIndex) * numNodes;
            theta += uK[ind] * phi[a];
        }
        double q = BCsConvection->getConvection(timeRun, dt, theta);
        double h = BCsConvection->getConvCons(timeRun, dt); //Convection constant h
        
        for (int a = 0; a < neighbours.size(); a++) {
            int indi11 = neighbours[a]+(numMechDofsPerNode+ fieldIndex) * numNodes;
            double val  = wJ * phi[a] * Jacobian * q; //Convection flux
            force[indi11]+= val;
        }
        
        // add to stiffness matrix
        if (tangentEstimation)
        {
            for (int a = 0; a < neighbours.size(); a++) {
                int indi11 = neighbours[a]+(numMechDofsPerNode+ fieldIndex) * numNodes;
                for (int b = 0; b < neighbours.size(); b++) {
                    double val = 0;
                    int indi12 = neighbours[b]+(numMechDofsPerNode+ fieldIndex) * numNodes;
                    val += wJ * phi[a] * Jacobian * h * phi[b];
                    
                    if (stiffBC.find(pair<int,int>(indi11,indi12)) == stiffBC.end())
                    {
                        stiffBC[pair<int,int>(indi11,indi12)] = val;
                    }
                    else
                    {
                        stiffBC[pair<int,int>(indi11,indi12)] += val;
                    }
                }
            }
            
        }

        delete tempGPF;
        delete tempGP;
    }

}

static void BCsRadiation(const vector<double> &uK, vector<double>& force, vector<classNodes *> &nodes,
                         vector<classElements *> &elements, classRadiation *&BCsRadiation, int &ndim, double timeRun,
                         double dt, bool tangentEstimation, map<pair<int,int>,double>& stiffBC) {

    vector<classGPs *> boundaryGPs;
    vector<classElements *> elems = BCsRadiation->getListOfElements();
    defineBCGPs(boundaryGPs, elems, nodes, ndim, true);
    int fieldIndex = BCsRadiation->getFieldIndex();
    int numNodes = nodes.size();
    for (vector<classGPs *>::iterator ite = boundaryGPs.begin(); ite != boundaryGPs.end(); ++ite) {
        classGPs *tempGP = *ite;
        classGPs *tempGPF;
        vector<double> Point = tempGP->getXYZ();
        elements[tempGP->getBulkElement()]->defineShapeOnAPoint(tempGP->getCell(), tempGPF, nodes, ndim, Point);
        vector<double> phi = tempGP->getPhi();
        vector<int> neighbours = tempGP->getNeighbours();
        double weight = tempGP->getW();
        double J = tempGP->getJ();
        // For linear triangles and tetrahedra, the deformation gradient is constant for the whole element
        vector<double> FF(ndim * ndim, 0);
        tempGPF->getF(uK, nodes.size(), ndim, FF);
        double Jacobian = determinantTensor(FF, ndim); // Volume
        int numMechDofsPerNode = tempGP->getCurrentState().getNumMechanicalDofsPerNode();
        
        double theta = 0;
        int numNeighbours = neighbours.size();
        for (int a = 0; a < numNeighbours; a++) {
            int ind = neighbours[a]+(numMechDofsPerNode+ fieldIndex) * numNodes;
            theta += uK[ind] * phi[a];
        }

        double q = BCsRadiation->getRadiation(timeRun, dt, theta);
        double A = BCsRadiation->getRadiaCons(timeRun, dt); //Radiation constant A
        
        double wJ = tempGP->getWeightJ();
        for (int a = 0; a < neighbours.size(); a++) {
            int indi11 = neighbours[a]+(numMechDofsPerNode+ fieldIndex) * numNodes;
            double val = wJ * phi[a] * Jacobian * q;
            force[indi11]+= val; //Radiation flux
        }
        
        // add to stiffness matrix
        if (tangentEstimation)
        {
            for (int a = 0; a < neighbours.size(); a++) {
                int indi11 = neighbours[a]+(numMechDofsPerNode+ fieldIndex) * numNodes;
                for (int b = 0; b < neighbours.size(); b++) {
                    double val = 0;
                    int indi12  = neighbours[b]+(numMechDofsPerNode+ fieldIndex) * numNodes;
                    val += wJ * phi[a] * Jacobian * 4 * A * pow(theta, 3) * phi[b];
                    if (stiffBC.find(pair<int,int>(indi11,indi12)) == stiffBC.end())
                    {
                        stiffBC[pair<int,int>(indi11,indi12)] = val;
                    }
                    else
                    {
                        stiffBC[pair<int,int>(indi11,indi12)] += val;
                    }
                }
            }
        }

        delete tempGPF;
        delete tempGP;
    }

}



/*! \brief This function applies the nodal forces (Neumann boundary condition)
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] NeumannNodes Neumann boundary nodes
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsNodalForces(vector<classNodes *> &nodes, vector<classElements *> &elements,
                           vector<classNeumann *> &NeumannNodes, int &ndim, vector<double> &force, double timeRun,
                           double dt) {

    int numNodes = nodes.size();
    int numNNodes = NeumannNodes.size();
    vector<double> F;
    int me = 0;
    for (int j = 0; j < numNNodes; j++) {
        F = NeumannNodes[j]->getTracAmount(timeRun, dt);
        me = NeumannNodes[j]->getMe();
        for (int k = 0; k < ndim; k++) {
            force[me + k * numNodes] += F[k];
        }
    }
}

/*! \brief This function applies the gravity forces (Neumann boundary condition)
  @param[in] nodes Array with all nodes in the domain
  @param[in] elements Array with all elements (or background cells) in the domain
  @param[in] NeumannNodes Neumann boundary nodes
  @param[in] ndim Dimension of the domain
  @param[inout] force Array with the external force applied until timeRun
  @param[in] timeRun Current simulation time
  @param[in] dt Time step
*/
static void BCsGravity(vector<classNodes *> &nodes, vector<classElements *> &elements, int &ndim, vector<double> &force,
                       classGravity *&BCsGrav) {
    int numNodes = nodes.size();
    vector<classGPs *> GPs;
    double rho;
    double g;
    g = BCsGrav->getGravity();
    for (int i = 0; i < elements.size(); i++) {
        elements[i]->defineGPs(i, GPs, nodes, ndim);
    }

    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) {

        classGPs *GaussPoint = *ite;
        double w = GaussPoint->getW();
        double Jacobian = GaussPoint->getJ();
        vector<double> phi = GaussPoint->getPhi();
        vector<constitutiveModels *> GaussPointCons = GaussPoint->getConstitutiveManager().getConsModel();
        rho = GaussPointCons[0]->getRho();
        vector<int> neighbours = GaussPoint->getNeighbours(); // GaussPointoral list of nodes in which a given GPs lies
        int numNeighbours = neighbours.size();
        for (int a = 0; a < numNeighbours; a++) {
            int id1 = neighbours[a] + (ndim - 1) * numNodes;
            force[id1] += -w * rho * Jacobian * g * phi[a];
        }
    }
}


#ifdef PARALLEL
/*! \brief This function finds all dof that are active during the simulation
  @param[in] DirichletNodes Dirichlet boundary nodes
  @param[in] ndim Dimension of the domain
  @param[in] numDof Number of degrees of freedom
  @param[out] activeDof Array with the labels of the active dof
*/
extern void findActiveDof(vector<classDirichlet *> DirichletNodes, int ndim, int numDof, vector<int> &activeDof,
                          vector<int> nod_local, vector<int> &nDof_local, int nNod) {

    activeSearching(DirichletNodes, ndim, numDof, activeDof);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < nod_local.size(); j++) {
            nDof_local.push_back((i * nNod) + j);
        }
    }

}
#endif
