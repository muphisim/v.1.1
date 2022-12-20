//
//
// File authors:  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//
/*!\file hyperElasticNeoHookean.cpp
  \brief This file contains all functions related to the hyperelastic neoHookean constitutive model. There are two type of functions: the virtual ones inheritated from the general constitutive law (they MUST be implemented by the user), and the particular ones that the user wants to implement for his convenience.
*/
#include "hyperElasticNeoHookean.h"
#include "maths.h"
#include "classGPs.h"
#include "TensorOperations.h"

/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
*/
classHyperElasticNeoHookean::classHyperElasticNeoHookean(int ndim, string name, string name_cons, double rho,
                                                         double Young, double nu) : 
                        constitutiveModels(ndim, name, name_cons, rho),
                        _Young(Young),
                        _nu(nu),
                        _lambda((Young * nu) / ((1. + nu) * (1. - 2. * nu))),
                        _mu(0.5 * Young / (1. + nu)),
                        _K1(Young / (3 * (1 - 2 * nu))) 
{
    
}

classHyperElasticNeoHookean::~classHyperElasticNeoHookean()
{
   
}

void classHyperElasticNeoHookean::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalTerm();
};

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classHyperElasticNeoHookean::soundSpeed() const
{
    return sqrt((_lambda + 2 * _mu) / _rho);
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classHyperElasticNeoHookean::initIntVars(vector<double> &intVars) {
    // This constitutive model does not have internal variables

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classHyperElasticNeoHookean::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);

    classTensor2 Ftensor(ndim, FCurr); // construct tensor from vector
    static classTensor2 B(ndim), invF(ndim);
    
    classTensor4* DBDFp = NULL;
    classTensor4* DFinvDFp = NULL;
    classTensor2* DJacobianDFp=NULL;
    double Jacobian = 1;
    if (flagTanMod)
    {
        static classTensor4 DBDF(ndim), DFinvDF(ndim);
        static classTensor2 DJacobianDF(ndim);
        TensorOperations::leftCauchy(Ftensor,B, &DBDF);
        Jacobian = TensorOperations::determinant(Ftensor,&DJacobianDF);
        TensorOperations::inverseTensor2(Ftensor,invF,&DFinvDF,false);
        DBDFp = &DBDF;
        DFinvDFp = &DFinvDF;
        DJacobianDFp = &DJacobianDF;
    }
    else
    {
        TensorOperations::leftCauchy(Ftensor,B);
        Jacobian = TensorOperations::determinant(Ftensor);
        TensorOperations::inverseTensor2(Ftensor,invF);
    }
    double traceB = B.trace();
    static classTensor2 I(ndim,1.); // unit tensor
    
    static classTensor2 Cauchy(ndim);
    Cauchy.setAll(0);
    double fact1 = pow(Jacobian, 5.0/double(ndim));
    double fact2 = _mu/fact1;
    Cauchy.axpy(B,fact2);
    double fact3 = -(_mu / double(ndim)) * traceB/fact1+ _K1 * (Jacobian - 1.0); 
    Cauchy.axpy(I,fact3);
    
    static classTensor2 PKtensor(ndim);
    TensorOperations::dotSecondTranspose(Cauchy,invF,PKtensor,Jacobian); // PK = J*sigma*F^-T
    PK1 = PKtensor.getData();
    cauchy = Cauchy.getData();
    computeVMS(cauchy, ndim, VMS);
    if (flagTanMod) 
    {
        double Dfact1DJacobian =  (5.0 / double(ndim))*pow(Jacobian, 5.0 / double(ndim)-1.);
        double Dfact2DJacobian = -_mu*Dfact1DJacobian/(fact1*fact1);
        double Dfact3DtraceB = -(_mu/double(ndim))/fact1;
        double Dfact3DJacobian = (_mu/double(ndim))*traceB*Dfact1DJacobian/(fact1*fact1)+ _K1; 
    
        static classTensor4 DCauchyDF(ndim);
        DCauchyDF = (*DBDFp);
        DCauchyDF.scale(fact2);
        TensorOperations::dyadicProductAdd(B,*DJacobianDFp,DCauchyDF,Dfact2DJacobian);
        
        static classTensor2 DtraceBDF(ndim);
        TensorOperations::doubleDot(I,*DBDFp,DtraceBDF); // DtraceB/DB= I
        TensorOperations::dyadicProductAdd(I,DtraceBDF,DCauchyDF,Dfact3DtraceB);
        TensorOperations::dyadicProductAdd(I,*DJacobianDFp,DCauchyDF,Dfact3DJacobian);
        
        classTensor4 *tanModuli = currState.getTangent();
        TensorOperations::dyadicProduct(PKtensor,*DJacobianDFp,*tanModuli,1./Jacobian);
        for (int i = 0; i < ndim; i++) 
        {
            for (int j = 0; j < ndim; j++) 
            {
                for (int k = 0; k < ndim; k++) 
                {
                    for (int l = 0; l < ndim; l++) 
                    {
                        for (int p=0; p< ndim; p++)
                        {
                            (*tanModuli)(i,j,k,l) += Jacobian*(DCauchyDF(i,p,k,l)*invF(j,p) + Cauchy(i,p)*DFinvDFp->operator()(j,p,k,l));
                        }
                    }
                }
            }
        }
    }
};



/*! \brief Constructor of the class
  @param[in] ndim Dimension of the domain
  @param[in] name Name of the constitutive model, defined by the user in the input file
  @param[in] rho Density
  @param[in] Young Young modulus
  @param[in] nu Poisson ratio
*/
classHyperElastic3DCompressibleNeoHookean::classHyperElastic3DCompressibleNeoHookean(int ndim, string name, string name_cons, double rho,
                                                         double Young, double nu) : 
                        constitutiveModels(ndim, name, name_cons, rho),
                        _Young(Young),
                        _nu(nu),
                        _lambda((Young * nu) / ((1. + nu) * (1. - 2. * nu))),
                        _mu(0.5 * Young / (1. + nu))
{
    
}

classHyperElastic3DCompressibleNeoHookean::~classHyperElastic3DCompressibleNeoHookean()
{
   
}

void classHyperElastic3DCompressibleNeoHookean::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalTerm();
};

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double classHyperElastic3DCompressibleNeoHookean::soundSpeed() const
{
    return sqrt((_lambda + 2 * _mu) / _rho);
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void classHyperElastic3DCompressibleNeoHookean::initIntVars(vector<double> &intVars) {
    // This constitutive model does not have internal variables

}

/*! \brief This function updates the constitutive response. This function is mandatory and it must be implemented by the user.
  @param[in] GaussPoint Gauss point
  @param[in] dt Time step
  @param[in] timeRun Current time
  @param[in] flagTanMod Flag to calculate the tangent modulus
*/
void classHyperElastic3DCompressibleNeoHookean::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod) 
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    classTensor2 Ftensor(ndim, FCurr); // construct tensor from vector
    classTensor2 C(ndim), invC(ndim);
    classTensor4 DinvCDC(ndim);
    double J=TensorOperations::determinant(Ftensor);
    TensorOperations::rightCauchy(Ftensor,C);
    if (flagTanMod)
    {
        TensorOperations::inverseTensor2(C,invC,&DinvCDC,true);
    }
    else
    {
        TensorOperations::inverseTensor2(C,invC);
    };
    classTensor2 PK2(ndim);
    //PK2=mu*(I-invC)+lambda*(logJ)*invC
    PK2.setAll(0);
    PK2.addDiagonal(_mu);
    PK2.axpy(invC,_lambda*log(J)-_mu);
    
    
    classTensor2 PKtensor(ndim);
    TensorOperations::dot(Ftensor,PK2,PKtensor); // PK1 = F*PK2
    PK1 = PKtensor.getData();
    computeCauchy(FCurr, ndim, PK1, cauchy);
    computeVMS(cauchy, ndim, VMS);
    if (flagTanMod) 
    {
        
        static classTensor4 DPK2DC(ndim);
        DPK2DC = (DinvCDC);
        DPK2DC.scale(_lambda*log(J)-_mu);
        TensorOperations::dyadicProductAdd(invC,invC,DPK2DC,0.5*_lambda);
        //
        classTensor4 *tanModuli = currState.getTangent();
        for (int i = 0; i < ndim; i++) {
            for (int J = 0; J < ndim; J++) {
                for (int k = 0; k < ndim; k++) {
                    for (int L = 0; L < ndim; L++) {
                        double val = 0;
                        for (int K = 0; K < ndim; K++) {
                            for (int I = 0; I < ndim; I++) {
                                val += 2.*DPK2DC(I, J, K, L) * FCurr[i * ndim + I] * FCurr[k * ndim + K];
                            }
                        }
                        if (i == k) {
                            val += PK2(J,L);
                        }
                        tanModuli->setValues(i, J, k, L, val);
                    }
                }
            }
        }

    }
};

double classHyperElastic3DCompressibleNeoHookean::defoEnergy(classGPs *GaussPoint) const
{
    // W = 0.5*lambda*log(J)**2 - mu*log(J)+ 0.5*mu*trace(FT*F)
    const classStates& currState = GaussPoint->getCurrentState();
    const vector<double> &FCurr = currState.getDeformationGradient();
    int ndim = getDim();
    classTensor2 Ftensor(ndim, FCurr); // construct tensor from vector
    classTensor2 C(ndim);
    
    double J=TensorOperations::determinant(Ftensor);
    TensorOperations::rightCauchy(Ftensor,C);
    double traceC = C.trace();
    return 0.5*_lambda*log(J)*log(J) - _mu*log(J) + 0.5*_mu*traceC;
}
