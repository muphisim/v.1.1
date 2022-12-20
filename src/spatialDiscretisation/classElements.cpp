//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file classElements.cpp
  \brief This file contains the definition of the constructors of the classElements
*/

#include "classElements.h"

/*! Constructor for classElements. It should not be called.

*/
classElements::classElements() {
    ERROR("This constructor should not be called");
    exit(-1);
};

/*! Constructor for classElements in case of subelements.
  @param[in] me Number of element
  @param[in] type Type of element
  @param[in] myNodes Nodes that build the element
  @param[in] ndim Dimension of the domain
  @param[in] bulkElement number of parent element to subelement
*/
classElements::classElements(int me, string type, const vector<int>& myNodes, int ndim, int bulkElement) {

    this->bulkElement = bulkElement;
    this->me = me;
    this->type = type;
    this->myNodes = myNodes;
    this->Cauchy.assign(ndim * ndim, 0);
    this->E.assign(ndim * ndim, 0);
    this->VMS.assign(1, 0);
    this->volumetricStrain.assign(1, 0);
    this->equivalentStrain.assign(1, 0);
    this->n0.assign(3, 0);
    this->localNdim = ndim;
};
