//
//
// File authors:  see Authors.txt
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file viscoElastic.cpp
  \brief This file contains all functions related to viscoElastic constitutive model.  The model is the original formulation proposed by Simo et al (nonlinear viscoelasticity, Simo J.C., Computational Inelasticity, Springer). There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/

#include "viscoElastic.h"
#include "maths.h"
#include "classGPs.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] K Bulk modulus
  @param[in] N Number of internal variables. Number of branches in the generalised Maxwell model
  @param[in] mus Array of mus for N branches. Its size is N+1, because the initial position mus[0] is the for the branch with the spring alone (Einfi)
  @param[in] etas Array of etas for N branches. Its size is N 
*/
classViscoElastic::classViscoElastic(int ndim, string name, string name_cons, double rho, double K, int N,
                                     const vector<double>& mus, const vector<double>& etas) : 
                                     constitutiveModels(ndim, name,  name_cons, rho) {

    this->_mu = mus;
    this->_eta = etas;
    this->_K = K;
    double sum_mu = 0.;

    for (int i = 0; i < _mu.size(); ++i) {
        sum_mu += _mu[i];
    }
    this->_E = 9. * _K * sum_mu / (3. * _K + sum_mu);
    this->_nu = (3. * _K - 2. * sum_mu) / 2. / (3. * _K + sum_mu);

    for (int i = 0; i < _eta.size(); ++i) {
        _tau.push_back(_eta.at(i) / _mu.at(i + 1));
    }

    this->_N = N;
}

classViscoElastic::~classViscoElastic()
{
   
}

void classViscoElastic::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term;
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classViscoElastic::soundSpeed() const{
    double factornu = (1. - _nu) / ((1. + _nu) * (1. - 2. * _nu));
    return sqrt(_E * factornu / _rho);
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classViscoElastic::initIntVars(vector<double> &intVars) {
    // This constitutive model has as many internal variables as maxwell branches. Then, for each branch there are two tensors Hn and Sn. Then the internal variables array will be built as follows:
    /*
      ndim=3;
      TotalSize=N*2*ndim*ndim; That reads: N internal variables with two ndim order tensors
     */
    int ndim = getDim();
    vector<double> Hn(ndim *ndim, 0.);
    vector<double> Sn(ndim *ndim, 0.);
    int finalSize = _N * (Hn.size() + Sn.size());

    this->_finalSize = finalSize;
    // It should be initialised to zero, these loops are just to illustrate how to access to the internal variables;
    for (int i = 0; i < _finalSize; i++) {
        intVars.push_back(0.);
    }


}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] Fn Current deformation gradient tensor
  @param[in] F0 Current deformation gradient tensor
  @param[in] intVarsPrev Array of the previous state of the internal variables
  @param[out] intVarsCurr Array of the current state of the internal variables
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] stiff Flag to calculate the tangent modulus
  @param[out] P First Piola Kirchhoff stress tensor
  @param[out] PK1Prev First Piola Kirchhoff stress tensor, previous state
  @param[out] tanModuli Tangent moduli,  Fourth order tensor
*/
void classViscoElastic::updateConstitutive(const vector<double> &Fn, const vector<double> &F0, const vector<double> &intVarsPrev,
                                           vector<double> &intVarsCurr, double dt, double timeRun, bool stiff,
                                           vector<double> &P, const vector<double> &PK1Prev, classTensor4 *tanModuli) {

    this->_timeStep = dt;
    this->update(F0, Fn, P, intVarsCurr, intVarsPrev, tanModuli, stiff);
}


/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classViscoElastic::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &Fn = currState.getDeformationGradient();
    const vector<double> &F0 = prevState.getDeformationGradient();
    vector<double> &E = currState.getGL();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    int ndim = getDim();
    computeGL(Fn, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    vector<double> &P = currState.getFirstPiolaKirchhoffStress();
    classTensor4 *tanModuli = currState.getTangent();
    this->_timeStep = dt;
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();
    this->update(F0, Fn, P, intVarsCurr, intVarsPrev, tanModuli, flagTanMod);

}

/*! \brief This function calculates the initial elastic Kirchhoff stress tensor
  @param[in] Fn Current deformation gradient tensor
  @param[out] Kirchhoff_ Initial Kirchooff stress tensor
*/
void classViscoElastic::initialKirchhoffStressTensor(const vector<double> &Fn, vector<double> &Kirchhoff_) const {
    int ndim = getDim();
    double Jn = determinantTensor(Fn, ndim);
    // Juli 2d
    //double frac = 1./3.;
    double frac = 1. / double(ndim);
    // STensor3 Fn_ = pow(Jn, -frac) * Fn; // volume-preserving deformation gradient
    vector<double> Fn_(ndim *ndim,
    0.);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Fn_[i * ndim + j] = pow(Jn, -frac) * Fn[i * ndim + j]; // defining the kronecker delta
        }
    }
    double sum_mu = 0.;

    for (int i = 0; i < _mu.size(); ++i) {
        sum_mu += _mu.at(i);
    }

    double derivative = 1. / 2. * sum_mu;
    vector<double> total(ndim *ndim,
    0.);
    multSTensor3SecondTranspose(Fn_, ndim, ndim, Fn_, ndim, total);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            total[i * ndim + j] = 2.0 * derivative * total[i * ndim + j];
        }
    }

    Kirchhoff_ = devSTensor3(total, ndim); // T = dev[mu * F_*F_']

}




/*! \brief This function updates the internal variables of the constitutive model
  @param[in] F1 Current deformation gradient tensor
  @param[in] kirchhoff1_ Initial Kirchooff stress tensor
  @param[out] IPvis Array of the current state of the internal variables
  @param[in] IPvisprev Array of the previous state of the internal variables
  @param[in] g Equal to r_inf = mu0/sumMu
  @param[inout] hn Algoritmic tensor
  @param[out] r Array defined r(i) = mu(i)/sumMu
  @param[out] tanPartD Tangent moduli corresponding to the viscoelastic part
  @param[out] ICinv Tensors needed to reconstruct the total tangent moduli
  @param[out] CinvTensorCinv Tensors needed to reconstruct the total tangent moduli
  @param[in] C Cauchy-Green tensor 
  @param[in] Cinv Cauchy-Green inverse tensor 
*/
// update internal variables
void classViscoElastic::updateIP(const vector<double> &F1, const vector<double> &kirchhoff1_, vector<double> &IPvis,
                                 const vector<double> &IPvisprev, double &g, vector<double> &hn, vector<double> &r,
                                 classTensor4 *&tanPartD, classTensor4 *&ICinv, classTensor4 *&CinvTensorCinv,
                                 const vector<double> &C, const vector<double> &Cinv) const {
    // cout<< "tau "<<_tau[0]<<endl;
    // Juli 2d
    int ndim = getDim();
    double frac = 1. / double(ndim);
    double dt = _timeStep;
    double Jn1 = determinantTensor(F1, ndim);
    vector<double> F1_(ndim *ndim,
    0.);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            F1_[i * ndim + j] = pow(Jn1, -frac) * F1[i * ndim + j];
        }
    }
    vector<double> Inv_F1_(ndim *ndim,
    0.);
    inverse(F1_, ndim, Inv_F1_);
    vector<double> inter1(ndim *ndim,
    0.);
    fill(IPvis.begin(), IPvis.end(), 0);

    classTensor4 *tanPartC2, *tanPartC3, *tanPartC4;
    tanPartC2 = new classTensor4(ndim);
    tanPartC3 = new classTensor4(ndim);
    tanPartC4 = new classTensor4(ndim);

    for (int i = 0; i < _N; i++) { // for all the internal variables
        vector<double> Hn_c(ndim *ndim, 0.);
        vector<double> Hn_cDEV(ndim *ndim, 0.);
        vector<double> Hn_nplus1(ndim *ndim, 0.);
        vector<double> Sn_nplus1(ndim *ndim, 0.);
        vector<double> _Hn;
        vector<double> _Sn;
        _Hn.assign(&IPvisprev[i * (ndim * ndim + ndim * ndim)],
                   &IPvisprev[i * (ndim * ndim + ndim * ndim)] + ndim * ndim);
        _Sn.assign(&IPvisprev[i * (ndim * ndim + ndim * ndim) + ndim * ndim],
                   &IPvisprev[i * (ndim * ndim + ndim * ndim) + ndim * ndim] + ndim * ndim);

        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                Hn_c[j * ndim + k] =
                        exp(-dt / _tau.at(i)) * _Hn[j * ndim + k] + (-exp(-dt / (2. * _tau.at(i))) * _Sn[j * ndim + k]);
            }
        }
        multTensorTensor3(Inv_F1_, ndim, ndim, kirchhoff1_, ndim, inter1);
        multSTensor3SecondTranspose(inter1, ndim, ndim, Inv_F1_, ndim, Sn_nplus1);
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                Hn_nplus1[j * ndim + k] = Hn_c[j * ndim + k] + exp(-dt / (2. * _tau.at(i))) * Sn_nplus1[j * ndim + k];
            }
        }
        g = g + r.at(i) * exp(-dt / (2 * _tau.at(i))); //g = r_inf + sum(r(i)*exp(-dt/(2*tau(i)))
        vector<double> total(ndim *ndim, 0.);
        vector<double> intertotal(ndim *ndim, 0.);
        vector<double> temp(ndim *ndim, 0.);
        multTensorTensor3(F1_, ndim, ndim, Hn_c, ndim, intertotal);
        multSTensor3SecondTranspose(intertotal, ndim, ndim, F1_, ndim, total);
        temp = devSTensor3(total, ndim);
        for (int j = 0; j < ndim; j++) { // for all the internal variables
            for (int k = 0; k < ndim; k++) { // for all the internal variables
                hn[j * ndim + k] = hn[j * ndim + k] + r.at(i) * temp[j * ndim + k];
            }
        }
        Hn_cDEV = DEVSTensor3(Hn_c, ndim, C, Cinv);
        dyadicProductSecondTensorSecondTensor(Hn_cDEV, ndim, ndim, Cinv, ndim, ndim, tanPartC2);
        dyadicProductSecondTensorSecondTensor(Cinv, ndim, ndim, Hn_cDEV, ndim, ndim, tanPartC3);
        // double contraction
        double tempi = doubleContractionSecondTensors(Hn_c, C, ndim);
        double temp5 = 0;
        //double factorDon=1./double(ndim)-1;
        //double factorDon=-2./double(ndim);
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) { // for all the internal variables
                for (int l = 0; l < ndim; l++) { // for all the internal variables
                    for (int m = 0; m < ndim; m++) { // for all the internal variables
                        tanPartC2->setValues(j, k, l, m, pow(Jn1, -2. / double(ndim)) * r.at(i) * (-2. / double(ndim)) *
                                                         tanPartC2->get(j, k, l, m));
                        tanPartC3->setValues(j, k, l, m, pow(Jn1, -2. / double(ndim)) * r.at(i) * (-2. / double(ndim)) *
                                                         tanPartC3->get(j, k, l, m));
                        tanPartC4->setValues(j, k, l, m,
                                             pow(Jn1, -2. / double(ndim)) * r.at(i) * (2. / double(ndim)) * tempi *
                                             (ICinv->get(j, k, l, m) -
                                              1. / double(ndim) * CinvTensorCinv->get(j, k, l, m)));

                        //tanPartC2->setValues(j,k,l,m,pow(Jn1, factorDon)*r.at(i)*(factorDon)*tanPartC2->get(j,k,l,m));
                        // tanPartC3->setValues(j,k,l,m,pow(Jn1, factorDon)*r.at(i)*(factorDon)*tanPartC3->get(j,k,l,m));
                        // tanPartC4->setValues(j,k,l,m,pow(Jn1, factorDon)*r.at(i)*(-factorDon)*tempi*(ICinv->get(j,k,l,m)-1./double(ndim)*CinvTensorCinv->get(j,k,l,m)));
                        temp5 = tanPartD->get(j, k, l, m) + tanPartC2->get(i, j, k, l) + tanPartC3->get(i, j, k, l) +
                                tanPartC4->get(i, j, k, l);
                        tanPartD->setValues(j, k, l, m, temp5);
                    }
                }
            }
        }
        //     hn = hn +r.at(i)* devSTensor3(total); //hn = sum(r(i)*dev[F_*Hn_c(i)*F_'])
        for (int j = 0; j < ndim; j++) { // for all the internal variables
            for (int k = 0; k < ndim; k++) { // for all the internal variables
                IPvis[i * (ndim * ndim + ndim * ndim) + j * ndim + k] = Hn_nplus1[j * ndim + k];
                IPvis[i * (ndim * ndim + ndim * ndim) + ndim * ndim + j * ndim + k] = Sn_nplus1[j * ndim + k];
            }
        }
    }
    delete tanPartC2;
    delete tanPartC3;
    delete tanPartC4;
}

/*! \brief This function updates the constitutive response
  @param[in] F0 Previous deformation gradient tensor
  @param[in] F Current deformation gradient tensor
  @param[out] P First Piola Kirchhoff stress tensor
  @param[out] IPvis Array of the current state of the internal variables
  @param[in] IPvisprev Array of the previous state of the internal variables
  @param[out] tanModuli Tangent moduli
  @param[in] stiff Flag to calculate (or not) the tangent moduli
*/
void classViscoElastic::update(const vector<double> &F0, const vector<double> &F, vector<double> &P, vector<double> &IPvis,
                               const vector<double> &IPvisprev, classTensor4 *&tanModuli, bool stiff) const {

    int ndim = getDim();
    vector<double> Finv(ndim *ndim, 0.);
    inverse(F, ndim, Finv);
    vector<double> T(ndim *ndim, 0.); // Kirchhoff stress
    vector<double> I(ndim *ndim, 0.);
    classTensor4 *tanModuliNoNominal; // composed by parts a, b and c. See Julian's handnotes
    tanModuliNoNominal = new classTensor4(ndim);
    for (int i = 0; i < ndim; i++) {
        I[i * ndim + i] = 1.; // defining the kronecker delta
    }
    this->stress(T, F, IPvis, IPvisprev, tanModuliNoNominal); // calculate Kirchhoff stress
    //  This tanModuli will be dS/dE. It should be transformed to dP/dF that is what is required by MuPhiSim
    // cout<< "TanModuliNoNominal"<<endl;
    // for (int i = 0; i < ndim; i++){
    //   for (int j = 0; j < ndim; j++){
    //     for (int k = 0; k < ndim; k++){
    // 	for (int l = 0; l < ndim; l++){
    // 	  cout<<"i="<< i<<j<<k<<l<< " "<<tanModuliNoNominal->get(i,j,k,l)<<endl;
    // 	}
    //     }
    //   }
    // }
    // exit(0);
    multSTensor3SecondTranspose(T, ndim, ndim, Finv, ndim, P);




    if (stiff)
    { 

        vector<double> S(ndim *ndim,
        0.);
        multTensorTensor3(Finv, ndim, ndim, P, ndim, S);
        classTensor4 *finalTanPartA, *finalTanPartB, *finalTanPartBA, *finalTanPartBB, *finalTanPartBBFinal, *tempFour, *tempFour2;
        finalTanPartA = new classTensor4(ndim);
        finalTanPartB = new classTensor4(ndim);
        finalTanPartBA = new classTensor4(ndim);
        finalTanPartBB = new classTensor4(ndim);
        finalTanPartBBFinal = new classTensor4(ndim);
        tempFour = new classTensor4(ndim);
        tempFour2 = new classTensor4(ndim);

        // This tensor product is different. JULIAN please take care of this
        tensorProductSecondTensorSecondTensor(I, ndim, ndim, S, ndim, ndim, finalTanPartA);
        tensorProductSecondTensorSecondTensor(F, ndim, ndim, I, ndim, ndim, finalTanPartBA);
        fourthTensorTranspose(ndim, finalTanPartBA, finalTanPartBB);

        // cout<< "finalA"<<endl;
        // conta2=0;
        // for (int i = 0; i < ndim; i++){
        //   for (int j = 0; j < ndim; j++){
        //     for (int k = 0; k < ndim; k++){
        // 	for (int l = 0; l < ndim; l++){
        // 	  cout<<"i="<< i<<j<<k<<l<< " "<<finalTanPartA->get(i,j,k,l)<<endl;
        // 	  conta2++;
        // 	}
        //     }
        //   }
        // }

        // cout<< "finalBA"<<endl;
        // conta2=0;
        // for (int i = 0; i < ndim; i++){
        //   for (int j = 0; j < ndim; j++){
        //     for (int k = 0; k < ndim; k++){
        // 	for (int l = 0; l < ndim; l++){
        // 	  cout<<"i="<< i<<j<<k<<l<< " "<<finalTanPartBA->get(i,j,k,l)<<endl;
        // 	  conta2++;
        // 	}
        //     }
        //   }
        // }

        // cout<< "finalBB"<<endl;
        // conta2=0;
        // for (int i = 0; i < ndim; i++){
        //   for (int j = 0; j < ndim; j++){
        //     for (int k = 0; k < ndim; k++){
        // 	for (int l = 0; l < ndim; l++){
        // 	  cout<<"i="<< i<<j<<k<<l<< " "<<finalTanPartBB->get(i,j,k,l)<<endl;
        // 	  conta2++;
        // 	}
        //     }
        //   }
        // }



        // cout<< "Second Piola"<<endl;
        // for(int i=0;i<ndim;i++){
        //   for(int j=0;j<ndim;j++){
        // 	cout<<S[i*ndim+j]<<" ";
        //   }
        //   cout<<endl;
        // }
        multFourthTensorFourthTensor(ndim, tanModuliNoNominal, finalTanPartBB, tempFour);
        // cout<< "tempFour"<<endl;
        // conta2=0;
        // for (int i = 0; i < ndim; i++){
        //   for (int j = 0; j < ndim; j++){
        //     for (int k = 0; k < ndim; k++){
        // 	for (int l = 0; l < ndim; l++){
        // 	  cout<<"i="<< i<<j<<k<<l<< " "<<tempFour->get(i,j,k,l)<<endl;
        // 	  conta2++;
        // 	}
        //     }
        //   }
        // }
        multFourthTensorFourthTensor(ndim, finalTanPartBA, tempFour, tempFour2);

        // cout<< "tempFour2"<<endl;
        // conta2=0;
        // for (int i = 0; i < ndim; i++){
        //   for (int j = 0; j < ndim; j++){
        //     for (int k = 0; k < ndim; k++){
        // 	for (int l = 0; l < ndim; l++){
        // 	  cout<<"i="<< i<<j<<k<<l<< " "<<tempFour2->get(i,j,k,l)<<endl;
        // 	  conta2++;
        // 	}
        //     }
        //   }
        // }

        double tempo = 0;
        //cout<< "TanModuli nominal"<<endl;
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
                for (int k = 0; k < ndim; k++) {
                    for (int l = 0; l < ndim; l++) {
                        tempo = finalTanPartA->get(i, j, k, l) + tempFour2->get(i, j, k, l);
                        tanModuli->setValues(i, j, k, l, tempo);
                        //cout<<"i="<< i<<j<<k<<l<< " "<<tanModuli->get(i,j,k,l)<<endl;
                        // conta++;
                    }
                }
            }
        }

        // cout<< "Piola"<<endl;
        //  for(int i=0;i<ndim;i++){
        //    for(int j=0;j<ndim;j++){
        //    	cout<<P[i*ndim+j]<<" ";
        //     }
        //     cout<<endl;
        //  }

        delete finalTanPartA;
        delete finalTanPartB;
        delete finalTanPartBA;
        delete finalTanPartBB;
        delete finalTanPartBBFinal;
        delete tempFour;
        delete tempFour2;
    }
    delete tanModuliNoNominal;
}

/*! \brief This function calculates the Kirchhooff stress tensor
  @param[out] Sig_ Kirchhoff stress tensor
  @param[in]  Fn Current deformation gradient tensor
  @param[out] IPvis Array of the current state of the internal variables
  @param[in] IPvisprev Array of the previous state of the internal variables
  @param[out] tanModuli Tangent moduli
*/
void classViscoElastic::stress(vector<double> &Sig_, const vector<double> &Fn, vector<double> &IPvis, const vector<double> &IPvisprev,
                          classTensor4 *&tanModuli) const {
    int ndim = getDim();
    double Jn = determinantTensor(Fn, ndim);
    // volumetric part ( this is for a free energy function Psi=1/2 k (1/2 J^2 -1) - ln(J) )
    double P = 0.5 * _K * (Jn - 1. / Jn);
    double Up = P; // Julian's style
    double Upp = 0.5 * _K * (1 + 1. / (Jn * Jn));
    vector<double> I(ndim *ndim,
    0.);
    for (int i = 0; i < ndim; i++) {
        I[i * ndim + i] = 1.;
    }
    vector<double> C(ndim *ndim, 0.);
    vector<double> Cinv(ndim *ndim, 0.);
    multSTensor3FirstTranspose(Fn, ndim, ndim, Fn, ndim, C);
    inverse(C, ndim, Cinv);
    // tanModuliNoNominal is composed by parts a, b and c. See Julian's handnotes
    classTensor4 *tanPartA, *tanPartB, *tanPartC, *CinvTensorCinv, *ICinv, *ICinv2;
    classTensor4 *tanPartC1, *tanPartC12, *tanPartD;
    classTensor4 *tanPartC122, *tanPartC123, *tanPartC124;
    tanPartA = new classTensor4(ndim);
    tanPartB = new classTensor4(ndim);
    tanPartC = new classTensor4(ndim);
    CinvTensorCinv = new classTensor4(ndim);
    ICinv = new classTensor4(ndim);
    ICinv2 = new classTensor4(ndim);
    tanPartC1 = new classTensor4(ndim);
    tanPartC12 = new classTensor4(ndim);
    tanPartD = new classTensor4(ndim);
    tanPartC122 = new classTensor4(ndim);
    tanPartC123 = new classTensor4(ndim);
    tanPartC124 = new classTensor4(ndim);
    // tanModuli dS/dE (it is not the nominal dP/dF)
    // Check Julian's Notes 08/03/2016
    // tanModuli= 4*(tanPartA +tanPartB) + tanPartC
    // Value needed in several parts
    dyadicProductSecondTensorSecondTensor(Cinv, ndim, ndim, Cinv, ndim, ndim, CinvTensorCinv);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                for (int l = 0; l < ndim; l++) {
                    // tanPartA
                    tanPartA->setValues(i, j, k, l, Upp * Jn * 0.5 * CinvTensorCinv->get(i, j, k, l));
                    // Value needed in several parts
                    //ICinv->setValues(i,j,k,l,-Cinv[i*ndim+k]*Cinv[l*ndim+j]);
                    ICinv->setValues(i, j, k, l, 0.5 * (Cinv[i * ndim + k] * Cinv[j * ndim + l] +
                                                        Cinv[j * ndim + k] * Cinv[i * ndim + l]));
                    // ICinv2->setValues(i,j,k,l,-Cinv[i*ndim+k]*Cinv[l*ndim+j]);
                    tanPartB->setValues(i, j, k, l, Up * Jn * 0.5 *
                                                    (0.5 * CinvTensorCinv->get(i, j, k, l) - ICinv->get(i, j, k, l)));
                }
            }
        }
    }

    //STensor3 Part1 = Jn*P * I; //volumetric part = JPI
    vector<double> Part1(ndim *ndim,
    0.);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Part1[i * ndim + j] = Jn * P * I[i * ndim + j];
        }
    }
    double SumMu = 0.;
    vector<double> r;
    for (int i = 0; i < _mu.size(); i++) {
        SumMu += _mu.at(i); //_mu.at(0) = mu0 for t=inf
    }
    for (int i = 0; i < _N; i++) {
        r.push_back(_mu.at(i + 1) / SumMu); //r(i) = mu(i)/sumMu
    }
    double r_inf = _mu.at(0) / SumMu; //r_inf = mu0/sumMu
    double g = r_inf;
    vector<double> IniKirchhoff_(ndim *ndim, 0.);
    this->initialKirchhoffStressTensor(Fn, IniKirchhoff_);
    vector<double> hn(ndim *ndim, 0.);
    double derivative = 1. / 2. * SumMu;

    updateIP(Fn, IniKirchhoff_, IPvis, IPvisprev, g, hn, r, tanPartD, ICinv, CinvTensorCinv, C, Cinv);


    vector<double> Part2(ndim *ndim, 0.);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Part2[i * ndim + j] = g * IniKirchhoff_[i * ndim + j];
        }
    }
    // Sig_ = Part1 + Part2 + hn; //Kirchhoff stress
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            Sig_[i * ndim + j] = Part1[i * ndim + j] + Part2[i * ndim + j] + hn[i * ndim + j];
        }
    }
    // Tangent moduli no nominal
    // tanPartC1 and tanPartC12
    vector<double> dWdC(ndim *ndim, 0.);
    vector<double> S0Bar(ndim *ndim, 0.);
    // double factorDon=1./double(ndim)-1;
    //double factorDon=-2./double(ndim);


    for (int j = 0; j < ndim; j++) { // for all the internal variables
        dWdC[j * ndim + j] = 2. * derivative * I[j * ndim + j];
    }
    double temp23 = doubleContractionSecondTensors(dWdC, C, ndim);
    S0Bar = DEVSTensor3(dWdC, ndim, C, Cinv);
    for (int j = 0; j < ndim; j++) { // for all the internal variables
        for (int k = 0; k < ndim; k++) { // for all the internal variables
            S0Bar[j * ndim + k] = pow(Jn, -2. / double(ndim)) * S0Bar[j * ndim + k];
            //S0Bar[j*ndim+k]=pow(Jn, factorDon)*S0Bar[j*ndim+k];
        }
    }
    //  tanPart121=0; The second derivative of the stored energety function is zero
    dyadicProductSecondTensorSecondTensor(S0Bar, ndim, ndim, Cinv, ndim, ndim, tanPartC122);
    dyadicProductSecondTensorSecondTensor(Cinv, ndim, ndim, S0Bar, ndim, ndim, tanPartC123);
    //cout<< "TanPartC"<<endl;
    for (int j = 0; j < ndim; j++) {
        for (int k = 0; k < ndim; k++) {
            for (int l = 0; l < ndim; l++) {
                for (int m = 0; m < ndim; m++) {
                    tanPartC122->setValues(j, k, l, m, -2. / double(ndim) * tanPartC122->get(j, k, l, m));
                    tanPartC123->setValues(j, k, l, m, -2. / double(ndim) * tanPartC123->get(j, k, l, m));
                    tanPartC124->setValues(j, k, l, m, pow(Jn, 2. / double(ndim)) * (2. / double(ndim)) * temp23 *
                                                       (ICinv->get(j, k, l, m) -
                                                        1. / double(ndim) * CinvTensorCinv->get(j, k, l, m)));
                    tanPartC12->setValues(j, k, l, m, tanPartC122->get(j, k, l, m) + tanPartC123->get(j, k, l, m) +
                                                      tanPartC124->get(j, k, l, m));
                    tanPartC1->setValues(j, k, l, m, g * tanPartC12->get(j, k, l, m));
                    tanPartC->setValues(j, k, l, m, tanPartC1->get(j, k, l, m) + tanPartD->get(j, k, l, m));
                    tanModuli->setValues(j, k, l, m, 4 * (tanPartA->get(j, k, l, m) + tanPartB->get(j, k, l, m)) +
                                                     tanPartC->get(j, k, l, m));
                }
            }
        }
    }

    delete tanPartA;
    delete tanPartB;
    delete tanPartC;
    delete CinvTensorCinv;
    delete ICinv;
    delete ICinv2;
    delete tanPartC1;
    delete tanPartC12;
    delete tanPartD;
    delete tanPartC122;
    delete tanPartC123;
    delete tanPartC124;

}


