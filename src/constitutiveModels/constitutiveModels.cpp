//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file constitutiveModels.cpp
  \brief This file contains all functions related to the base constitutive model class. In this wrapper, there are two type of functions: the virtual ones that MUST be implemented by the user for all new constitutive models, and the particular ones that are going to be called somewhere else in the program. The user must not modify those functions.
*/
#include "constitutiveModels.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
*/
constitutiveModels::constitutiveModels(int ndim, string name, string name_cons, double rho) : 
    _name(name),_name_cons(name_cons), _ndim(ndim), _rho(rho), _extraFieldIndexes(), _isInitialised(false),
    _term(NULL)
{
 
};

constitutiveModels::~constitutiveModels()
{
    if (_term) delete _term;
    _term = NULL;
}

