//
//
// File author(s):  see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//


#include "J2plasticityFiniteStrain.h"
#include "classGPs.h"
#include "classStates.h"
#include "TensorOperations.h"

IsotropicHardening* IsotropicHardening::allocate(const char what[], const std::vector<double>& data)
{
    INFO("allocation hardening law %s",what);
    printVector(data,"haardening data");
    if (strcmp(what,"LinISO") == 0 )
    {
        if (data.size() <2 )
        {
            ERROR("The number of parameters must be equal to 2 with linearIsotropicHardening");
            exit(-1);
        }
        return new linearIsotropicHardening(data[0],data[1]);
    }
    else if (strcmp(what,"SwiftISO") == 0 )
    {
        if (data.size() <3 )
        {
            ERROR("The number of parameters must be equal to 3 with SwiftIsotropicHardening");
            exit(-1);
        }
        return new SwiftIsotropicHardening(data[0],data[1],data[2]);
    } 
    else if (strcmp(what,"VoceISO") == 0 )
    {
        if (data.size() <3 )
        {
            ERROR("The number of parameters must be equal to 3 with VoceIsotropicHardening");
            exit(-1);
        }
        return new VoceIsotropicHardening(data[0],data[1],data[2]);
    } 
    else
    {
        ERROR("IsotropicHardening cannot be allocated with %s",what);
        exit(-1);
    }
};


J2plasticityFiniteStrain::J2plasticityFiniteStrain(int ndim, string name, string name_cons, double rho, 
                                 double E, double nu, const IsotropicHardening& hardenLaw, int order):
                                 constitutiveModels(ndim, name, name_cons, rho),
                                 _E(E),_nu(nu), _K(E/3./(1.-2.*nu)), _mu(0.5*E/(1.+nu)),
                                 _intVarsSize(0), _order(order), _Cel(3), _I4dev(4)
{
    INFO("young = %g, nu = %g, order = %d",E,nu,order);
    _isoHardening = hardenLaw.clone();
    
    double lambda = (E*nu)/(1.+nu)/(1.-2.*nu);

    _Cel.setAll(0);
    _Cel(0,0,0,0) = lambda + 2*_mu;
    _Cel(1,1,0,0) = lambda;
    _Cel(2,2,0,0) = lambda;
    _Cel(0,0,1,1) = lambda;
    _Cel(1,1,1,1) = lambda + 2*_mu;
    _Cel(2,2,1,1) = lambda;
    _Cel(0,0,2,2) = lambda;
    _Cel(1,1,2,2) = lambda;
    _Cel(2,2,2,2) = lambda + 2*_mu;

    _Cel(1,0,1,0) = _mu;
    _Cel(2,0,2,0) = _mu;
    _Cel(0,1,0,1) = _mu;
    _Cel(2,1,2,1) = _mu;
    _Cel(0,2,0,2) = _mu;
    _Cel(1,2,1,2) = _mu;

    _Cel(0,1,1,0) = _mu;
    _Cel(0,2,2,0) = _mu;
    _Cel(1,0,0,1) = _mu;
    _Cel(1,2,2,1) = _mu;
    _Cel(2,0,0,2) = _mu;
    _Cel(2,1,1,2) = _mu;
    
    classTensor2 I(3,1.);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            for (int k=0; k<3; k++)
            {
                for (int l=0; l<3; l++)
                {
                    _I4dev(i,j,k,l) = 0.5*(I(i,k)*I(j,l)+I(i,l)*I(j,k)) - I(i,j)*I(k,l)/3.;
                }
            }
        }
    }
};

J2plasticityFiniteStrain::~J2plasticityFiniteStrain()
{
    if (_isoHardening) 
    {
        delete _isoHardening;
        _isoHardening = NULL;
    }
};


void J2plasticityFiniteStrain::initLawFromOptions()
{
    if (_isInitialised ) return;
    _isInitialised = true;
    if (_term) delete _term; 
    _term = new MechanicalTerm();
}

/*! \brief This function provides the sounds velocity in the material to calculate the critical time step
  @return Sound speed in the material
*/
double J2plasticityFiniteStrain::soundSpeed()  const
{
    double factornu = (1.-_nu)/((1.+_nu)*(1.-2.*_nu));
    double rho = this->getRho();
    double sound = sqrt(_E *factornu/ rho);
    return sound;
}

/*! \brief This function initialise the internal variables this constitutive model. This function is mandatory and it must be implemented by the user.
  @param[inout] intVars Array of the internal variables
*/
void J2plasticityFiniteStrain::initIntVars(vector<double> &intVars) 
{
    _intVarsSize = intVars.size();
    intVars.resize(_intVarsSize + 10, 0.); //  1 - plastic deformation 9 components of - plastic strain tensor 
    classTensor2 I(3,1.);
    for (int i=0; i<9; i++)
    {
        intVars[_intVarsSize+1+i] = I.getData()[i];
    }
};


void J2plasticityFiniteStrain::updateConstitutive(classGPs *GaussPoint, double dt, double timeRun, bool flagTanMod)
{
    classStates& currState = GaussPoint->getCurrentState();
    const classStates& prevState = GaussPoint->getPreviousState();
    //
    const vector<double> &FCurr = currState.getDeformationGradient();
    const vector<double> &FPrev = prevState.getDeformationGradient();
    vector<double> &PK1 = currState.getFirstPiolaKirchhoffStress();
    vector<double> &E = currState.getGL();
    vector<double> &Cauchy = currState.getCauchy();
    vector<double> &volumetricStrain = currState.getVolumetricStrain();
    vector<double> &equivalentStrain = currState.getEquivalentStrain();
    vector<double> &VMS = currState.getVMS();
    int ndim = getDim();
    computeGL(FCurr, ndim, E);
    computeVolumetricStrain(E, ndim, volumetricStrain);
    computeEquivalentStrain(E, ndim, equivalentStrain);
    vector<double>& intVarsCurr = currState.getInternalVariables();
    const vector<double>& intVarsPrev = prevState.getInternalVariables();

    classTensor2 F1tensor(ndim, currState.getDeformationGradient());
    classTensor2 F1(3,1.);
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            F1(i,j) = F1tensor(i,j);
        }
    }
    
    double p0 = intVarsPrev[_intVarsSize];
    classTensor2 Fp0(3, vector<double>(intVarsPrev.begin()+_intVarsSize+1, intVarsPrev.begin()+_intVarsSize+10)); 
    
    // elastic predictor
    classTensor2 Fp1(Fp0);
    double p1 = p0;
    static classTensor2 invFp0(3);
    static classTensor2 Fepr(3); // predictor of elastic strain
    TensorOperations::inverseTensor2(Fp0,invFp0);
    TensorOperations::dot(F1,invFp0,Fepr);
    static classTensor2 Fpinv(3), Fe(3), Ce(3), Ee(3);
    Fpinv = invFp0;
    Fe = Fepr;
    TensorOperations::rightCauchy(Fe,Ce);
    // compute strain by logarithmic operator
    static classTensor4 Lpr(3); // dlogCepr/dCepr
    static classTensor6 dLDCe(3,3); // ddlogCepr/ddCepr
    static classTensor2 Eepr(3); // 0.5*log(Cepr)
    TensorOperations::logClassTensor2(Ce,_order,Eepr,&Lpr,&dLDCe);
    Eepr.scale(0.5); // elastic strain
    Ee = Eepr;
    static classTensor4 DexpA(3); // use for compute tangent later
    static classTensor4 L(3);
    L = (Lpr);
    // compute corotational Kirchhoff stress
    double pres = _K*Eepr.trace();
    static classTensor2 corKirDev(3); // deviatoric part
    Eepr.dev(corKirDev);
    corKirDev.scale(2.*_mu);
    
    double Deps = 0.;
    double R=0, H=0;
    _isoHardening->hardening(p0,R,H);
    double R0 = _isoHardening->getR0();
    double sigVMpr = sqrt(1.5*TensorOperations::doubleDot(corKirDev,corKirDev));
    static classTensor2 N(3); // plastic normal
    N = corKirDev;
    N.scale(1.5/sigVMpr);
    double VMCr = sigVMpr-R;
    double tol = 1e-6;
    if (VMCr > tol*R0)
    {
        // plastic occurs
        int ite = 0;
        int maxInt = 1000;
        while (fabs(VMCr) > tol*R0)
        {
            ite++;
            double sol = VMCr/(3*_mu+H);
            if (Deps+sol <0)
            {
                Deps *=0.5;
            }
            else
            {
                Deps += sol;
            }
            _isoHardening->hardening(p0+Deps,R,H);
            VMCr = sigVMpr - 3.*_mu*Deps -R;
            
            if (ite > maxInt)
            {
                ERROR("plastic corrector does not convergce");
                PK1[0] = sqrt(-1); // to force the solver reduce time step
                return;
            };
        };
        
        // update plastic strain
        p1  = Deps+p0;        
        // update corotational kirchhoff stress
        static classTensor2 DepsN(3);
        DepsN = N;
        DepsN.scale(Deps);
          
        // update corKirDev
        corKirDev.axpy(DepsN,-2.*_mu);

        // update Fp
        classTensor2 expDepsN(3);
        TensorOperations::expClassTensor2(DepsN,_order,expDepsN,&DexpA,NULL);
        TensorOperations::dot(expDepsN,Fp0,Fp1);
        
        // update elastic strain
        TensorOperations::inverseTensor2(Fp1,Fpinv);
        TensorOperations::dot(F1,Fpinv,Fe);
        TensorOperations::rightCauchy(Fe,Ce);
        TensorOperations::logClassTensor2(Ce,_order,Ee,&L,&dLDCe);
        Ee.scale(0.5);
    };
    
    // estimation of PK stress
    classTensor2 corKir(3);
    corKir = corKirDev;
    corKir.addDiagonal(pres);
    
        
    classTensor2 S(3), FeS(3);
    TensorOperations::doubleDot(corKir,L,S);
    TensorOperations::dot(Fe,S,FeS);
    classTensor2 P1(3); 
    TensorOperations::dotSecondTranspose(FeS,Fpinv,P1);
    // update
    classTensor2 P1dim(ndim,0);
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            P1dim(i,j) = P1(i,j);
        }
    }
    PK1 = P1dim.getData();
    intVarsCurr[_intVarsSize] = p1;
    computeCauchy(FCurr, ndim, PK1, Cauchy);
    computeVMS(Cauchy, ndim, VMS);
    for (int i=0; i<9; i++)
    {
        intVarsCurr[_intVarsSize+i+1] =  Fp1.getData()[i];
    }
        
    if (flagTanMod)
    {
        // tangent computation
        static classTensor4 DcorKirDEepr(3);
        DcorKirDEepr =  (_Cel);
        static classTensor4 dFpDEepr(3);
        if (Deps > 0)
        {
            static classTensor2 devEepr(3);
            Eepr.dev(devEepr);
            double EeprEq = sqrt(TensorOperations::doubleDot(devEepr,devEepr)/1.5);
            static classTensor4 BE(3);
            TensorOperations::dyadicProduct(N,N,BE);
            BE.scale(2.*_mu/(3.*_mu+H)- 2.*Deps/3./EeprEq);
            double val=Deps/EeprEq;
            BE.axpy(_I4dev, val);

            // update with plasticity corKir
            DcorKirDEepr.axpy(BE,-2.*_mu);

            static classTensor4 DexpABE(3);
            TensorOperations::doubleDot(DexpA,BE,DexpABE);
            // update with plasticity Fp
            for (int i=0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                {
                    for (int k=0; k<3; k++)
                    {
                        for (int l=0; l<3; l++)
                        {
                            dFpDEepr(i,j,k,l) = 0.;
                            for (int p=0; p<3; p++)
                            {
                                dFpDEepr(i,j,k,l) += DexpABE(i,p,k,l)*Fp0(p,j);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            dFpDEepr.setAll(0);
        }

        static classTensor4 EprToF(3);
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                for (int k=0; k<3; k++){
                    for (int l=0; l<3; l++){
                        EprToF(i,j,k,l) = 0.;
                        for (int p=0; p<3; p++){
                            for (int q=0; q<3; q++){
                                EprToF(i,j,k,l) += Lpr(i,j,p,q)*Fepr(k,p)*invFp0(l,q);
                            }
                        }
                    }
                }
            }
        }

        static classTensor4 DcorKirDF(3), dFpdF(3);
        TensorOperations::doubleDot(DcorKirDEepr,EprToF,DcorKirDF);
        if (Deps>0.)
        {
            TensorOperations::doubleDot(dFpDEepr,EprToF,dFpdF);
        }
        else
        {
            dFpdF.setAll(0);
        }

        // done DcorKirDF, DcorKirDT, DFpDF, DFpDT

        static classTensor4 DinvFpdF(3), dFedF(3);
        if (Deps >0){
            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    for (int k=0; k<3; k++){
                        for (int l=0; l<3; l++){
                            DinvFpdF(i,j,k,l) = 0.;
                            for (int p=0; p<3; p++){
                                for (int q=0; q<3; q++){
                                    DinvFpdF(i,j,k,l) -= Fpinv(i,p)*dFpdF(p,q,k,l)*Fpinv(q,j);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            DinvFpdF.setAll(0);
        }

        dFedF.setAll(0);
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                for (int k=0; k<3; k++){
                    dFedF(i,j,i,k) += Fpinv(k,j);
                    if (Deps> 0){
                        for (int l=0; l<3; l++){
                            for (int p=0; p<3; p++){
                                dFedF(i,j,k,l) += F1(i,p)*DinvFpdF(p,j,k,l);
                            }
                        }
                    }
                }
            }
        }

        static classTensor6 DLDF(3,3);
        if (_order != 1)
        {
            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    for (int k=0; k<3; k++){
                        for (int l=0; l<3; l++){
                            for (int p=0; p<3; p++){
                                for (int q=0; q<3; q++){
                                    DLDF(i,j,k,l,p,q) = 0.;
                                    for (int r=0; r<3; r++){
                                        for (int s=0; s<3; s++){
                                            for (int a=0; a<3; a++){
                                                DLDF(i,j,k,l,p,q) += dLDCe(i,j,k,l,r,s)*2.*Fe(a,r)*dFedF(a,s,p,q);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            DLDF.setAll(0);
        }
        //

        static classTensor4 DSDF(3); // S = corKir:L
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                for (int k=0; k<3; k++){
                    for (int l=0; l<3; l++){
                        DSDF(i,j,k,l) = 0.;
                        for (int r=0; r<3; r++){
                            for (int s=0; s<3; s++){
                                DSDF(i,j,k,l) += DcorKirDF(r,s,k,l)*L(r,s,i,j) ;
                                if (_order != 1){
                                    DSDF(i,j,k,l) += corKir(r,s)*DLDF(r,s,i,j,k,l);
                                }
                            }
                        }
                    }
                }
            }
        }


        classTensor4& Tangent = *(currState.getTangent());
        Tangent.setAll(0);
        for (int i=0; i<ndim; i++){
            for (int j=0; j<ndim; j++){
                for (int k=0; k<3; k++){
                    for (int l=0; l<3; l++){
                        for (int p=0; p<ndim; p++){
                            for (int q=0; q<ndim; q++){
                                Tangent(i,j,p,q) += (dFedF(i,k,p,q)*S(k,l)*Fpinv(j,l) + Fe(i,k)*DSDF(k,l,p,q)*Fpinv(j,l)+Fe(i,k)*S(k,l)*DinvFpdF(j,l,p,q));
                            }
                        }
                    }
                }
            }
        }
    }
};