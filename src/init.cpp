//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file init.cpp
  \brief All functions related to the initialisation of the problem will be placed in this file
*/
#include "init.h"
#include "geometry.h"
#include "commandLine.h"

/*! \brief This function initialises the position, weights and contributions of the GPs based on background cells for MM and on the elements for FEM
    @param[inout] GPs Array with all GPs in the domain
    @param[in] nod Array with all nodes in the domain
    @param[in] elem Array with all elements (or background cells) in the domain
    @param[in] ndim Dimension of the domain
    @param[in] sim Type of simulation executed
    @param[in] nodFEM Array with the labels that correspond to the FEM nodes in the nodes general array
    @param[in] nodMM Array with the labels that correspond to the MM nodes in the nodes general array
    @param[in] nodShare Array with the labels that correspond to the shared nodes in the nodes general array
    @param[in] elemFEM Array with the labels that correspond to the FEM elements in the elements general array
    @param[in] elemMM Array with the labels that correspond to the MM background cells in the elements general array
    @param[out] consMod Array with all constitutive models in the domain
*/
extern void
initGPsDefinition(vector<classGPs *> &GPs, vector<classNodes *> &nod, vector<classElements *> &elem, int ndim, int sim,
                  vector<int> &nodFEM, vector<int> &nodMM, vector<int> &nodShare, vector<int> &elemFEM,
                  vector<int> &elemMM, vector<constitutiveModels *> &consMod) {

    vector<classGPs *> GPsMM;
    if (sim == 2) {
        GPsDistributionFEM(GPs, nod, elem, ndim, nodFEM, elemFEM);
    } else {
        if (sim == 1) {
            // The treatment of the GPs at MM regions is a bit different. I have to place them first and then I will calculate the shape functions. The following function is to use the elements given in the input file as integration cells for MM
            GPsDistributionMM(GPs, nod, elem, ndim, nodMM, elemMM, GPsMM, GeneralOptions::MMIntegrationOrder);
            maxEntSFatGPs(GPsMM, nod, ndim, nodMM, true, false);
        } else if (sim == 3) {
            // The finite element region should be called first. Then, the MM region is called into the two different steps
            GPsDistributionMM(GPs, nod, elem, ndim, nodMM, elemMM, GPsMM, GeneralOptions::MMIntegrationOrder);
            GPsDistributionFEM(GPs, nod, elem, ndim, nodFEM, elemFEM);
            maxEntSFatGPs(GPsMM, nod, ndim, nodMM, true, false);
        }
    }
};

/*! \brief This function initialises the internal variables of the constitutive model assigned to each GP
    @param[inout] GPs Array with all GPs in the domain
    @param[in] ndim Dimension of the domain
*/
extern void initInternalVariables(vector<classGPs *> &GPs, vector<classElements *> &elem) {

    for (vector<classGPs *>::iterator ite = GPs.begin(); ite != GPs.end(); ++ite) 
    {
        classGPs *temp = *ite;
        const vector<int>& intVarsPosition = elem[temp->getCell()]->get_IntVarsPosition(); //To get the position of internal variables stored in the element to which the GP belongs
        const vector<double>& intVarsVal = elem[temp->getCell()]->get_IntVarsValue();  //To get the value of internal variables stored in the element to which the GP belongs
        // the internal variables are created by push_back into a vector
        // it must start from a zero-size vector
        vector<double> intVars(0);
        temp->getConstitutiveManager().initIntVars(intVars);
        for (int j = 0; j < intVarsPosition.size(); j++) 
        { //Loop to assign the value, defined in the inputFile
            if (intVarsPosition[j] > (intVars.size() - 1)) 
            {//error if the positions introduced by the users are greater than the size of the internal variables vector
                ERROR("Position of an internal variable not valid");
                exit(-1);
            } else 
            {
                intVars[intVarsPosition[j]] = intVarsVal[j]; //if not, we assign the value of the internal variable
            }
        }
        temp->getCurrentState().getInternalVariables() = (intVars);
        temp->getPreviousState().getInternalVariables() = (intVars);
        temp->getInitialState().getInternalVariables() = (intVars);
    }
};
