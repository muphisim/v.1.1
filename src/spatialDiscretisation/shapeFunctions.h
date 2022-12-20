//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _shapeFunctions_H_
#define _shapeFunctions_H_


/*!\file shapeFunctions.h
  \brief Several functions to calculate the shape functions for MM
*/
#include "configuration.h"
#include "supportDomains.h"
#include "classNodes.h"
#include "classGPs.h"
#include "init.h"
#include "maths.h"

extern void maxEntSFatGPs(vector<classGPs *> &GPs, vector<classNodes *> &nodes, int ndim, vector<int> &nodesMM, bool flagSupport,
              bool flagCurrCon);

extern void maxEntSFatNod(vector<classNodes *> &nodes, vector<classNodes *> &nodes2, int ndim, vector<int> &nodesMM,
                          bool flagSupport, bool flagCurrCon);

#endif

