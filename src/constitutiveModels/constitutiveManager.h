//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _CONSTITUTIVEMANAGER_H_
#define _CONSTITUTIVEMANAGER_H_

#include "constitutiveModels.h"

class constitutiveManager 
{
    protected:
        std::vector<constitutiveModels*> _allLaws;
        
    public:
        constitutiveManager();
        ~constitutiveManager();
        void clear();
        void addLaw(constitutiveModels* law);
        std::vector<constitutiveModels*>& getConsModel() {return _allLaws;};
        const std::vector<constitutiveModels*>& getConsModel() const {return _allLaws;};
        int getTotalNumDofPerNode() const;
        int getNumMechDofPerNode() const;
        int getNumExtraDofPerNode() const;
        void initIntVars(vector<double>& invars) const;
        void checkActivation(classGPs *GaussPoint, double timeRun) const;
        void constitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) const;
        void computeKinematicVariables(classGPs *GaussPoint, const vector<double>& uK, int numNodes) const;
        
        void computeInertialForce(int numNodes, const vector<double>& aK, classGPs *GaussPoint, vector<int>& vDofs, mVector& vForce) const;
        void computeInternalForce(int numNodes, classGPs *GaussPoint, vector<int>& vDofs, mVector& vForce) const;
        void computeStiffnessMatrix(double timeRun, double dt, vector<double>& uK, int numNodes, 
                                    classGPs *GaussPoint, vector<int>& vDofs, mMatrix& stiffness, bool byPerturbation, double tol) const;
        void computeMassMatrix(int numNodes, classGPs *GaussPoint, vector<int>& vDofs, mMatrix& mass) const;
        void computeLumpedMassVector(int numNodes, classGPs *GaussPoint, vector<int>& vDofs, mVector& mass) const;
        void computeInternalForceExplicit(int numNodes, classGPs *GaussPoint, vector<int>& vDofs, mVector& vForce, double timeRun, double dt) const;
};

#endif