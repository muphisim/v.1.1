//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _classTensorOperations_H_
#define _classTensorOperations_H_

#include "maths.h"

namespace TensorOperations
{

/*! \brief inverse of a tensor
  @param[in] A input
  @param[out] invA output
  @param[out] DinvADA output if this pointer is non-null
 */
inline void inverseTensor2(const classTensor2& A, classTensor2& invA, classTensor4* DinvADA =NULL, bool sym = true)
{
    int ndim = A.getDim();
    inverse(A.getData(),ndim,invA.getData());
    if (DinvADA !=NULL)
    {
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< ndim; j++)
            {
                for (int k=0; k< ndim; k++)
                {
                    for (int l=0; l< ndim; l++)
                    {
                        if (sym)
                        {
                            (*DinvADA)(i,j,k,l) = -0.5*(invA(i,k)*invA(j,l)+ invA(i,l)*invA(j,k));
                        }
                        else
                        {
                            (*DinvADA)(i,j,k,l)  = -invA(i,k)*invA(l,j);
                        }
                    }
                }
            }
        }
    }
};

/*! \brief Determinant J=det(A)
  @param[in] A input
  @param[out] DdetA/DA output if this pointer is non-null
   
 */
inline double determinant(const classTensor2& A, classTensor2* DJDA=NULL)
{
    int ndim = A.getDim();
    double detA =  determinantTensor(A.getData(),ndim);
    if (DJDA != NULL)
    {
        classTensor2 invA(ndim);
        TensorOperations::inverseTensor2(A, invA);
        // DdetADA = detA*A^{-T}
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< ndim; j++)
            {
                (*DJDA)(i,j) = detA*invA(j,i);
            }
        }
    }
    return detA;
};


/*! \brief dyadic of two tensors
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dyadicProduct(const classTensor1& A, const classTensor1& B, classTensor2& C, double alpha=1.)
{
    int ndim = A.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            C(i,j) = A(i)*B(j)*alpha;
            
        };
    };
};

/*! \brief dyadic of two tensors
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dyadicProduct(const classTensor2& A, const classTensor2& B, classTensor4& C, double alpha=1.)
{
    int ndim = A.getDim();
    C.setAll(0);
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                      C(i,j,k,l) = A(i,j)*B(k,l)*alpha;
                }
            }
        };
    };
};

/*! \brief dyadic of two tensors
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dyadicProductAdd(const classTensor2& A, const classTensor2& B, classTensor4& C, double alpha=1.)
{
    int ndim = A.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                      C(i,j,k,l) += A(i,j)*B(k,l)*alpha;
                }
            }
        };
    };
};

/*! \brief double dot of two tensors
  @param[in] A input
  @param[in] B input
 */
inline double doubleDot(const classTensor2& A, const classTensor2& B, double fact = 1.)
{
    double v = 0;
    int ndim = A.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            v += A(i,j)*B(i,j)*fact;
            
        };
    };
    return v;
};

/*! \brief double dot of two tensors
  @param[in] A input
  @param[in] B input
  @param[in] C output
 */
inline void doubleDot(const classTensor4& A, const classTensor4& B, classTensor4& C, double fact = 1.)
{
    int ndim = A.getDim();
    C.setAll(0);
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    for (int p=0; p< ndim; p++)
                    {
                        for (int q=0; q< ndim; q++)
                        {
                            C(i,j,k,l) += A(i,j,p,q)*B(p,q,k,l)*fact;
                        }
                    }
                    
                }
            }
        }
    }
};

/*! \brief transpose of a tensor
  @param[in] A input
  @param[out] AT output
 */
inline void transpose(const classTensor2& A, classTensor2& AT)
{
    int ndim = A.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            AT(i,j) = A(j,i);
        }
    }
};



/*! \brief Multiplication C = A.B*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dot(const classTensor2& A, const classTensor2& B, classTensor2& C, double fact =1.)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                C(i,j) += A(i,k)*B(k,j)*fact;
            }
        }
    }
};

/*! \brief Multiplication C = A.B*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
  @param[in] indexFirst
  @param[in] indexSecond indexSecond
 */
inline void dot(const classTensor2& A, const classTensor4& B, classTensor4& C, int indexFirst, int indexSecond, double fact =1.)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    for (int p=0; p< ndim; p++)
                    {
                        if (indexFirst == 0)
                        {
                            if (indexSecond == 0)
                            {
                                C(i,j,k,l) += A(p,i)*B(p,j,k,l)*fact;
                            }
                            else if (indexSecond == 1)
                            {
                                C(i,j,k,l) += A(p,i)*B(j,p,k,l)*fact;
                            }
                            else if (indexSecond == 2)
                            {
                                C(i,j,k,l) += A(p,i)*B(j,k,p,l)*fact;
                            }
                            else if (indexSecond == 3)
                            {
                                C(i,j,k,l) += A(p,i)*B(j,k,l,p)*fact;
                            }
                        }
                        else if (indexFirst == 1)
                        {
                            if (indexSecond == 0)
                            {
                                C(i,j,k,l) += A(i,p)*B(p,j,k,l)*fact;
                            }
                            else if (indexSecond == 1)
                            {
                                C(i,j,k,l) += A(i,p)*B(j,p,k,l)*fact;
                            }
                            else if (indexSecond == 2)
                            {
                                C(i,j,k,l) += A(i,p)*B(j,k,p,l)*fact;
                            }
                            else if (indexSecond == 3)
                            {
                                C(i,j,k,l) += A(i,p)*B(j,k,l,p)*fact;
                            }
                        }
                         
                    }
                }
               
            }
        }
    }
};

/*! \brief Multiplication D += A.B.C*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dotThreeTensorsAdd(const classTensor2& A, const classTensor2& B, const classTensor2& C, classTensor2& D, double fact =1.)
{
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    D(i,j) += A(i,k)*B(k,l)*C(l,j)*fact;
                }
            }
        }
    }
};

/*! \brief Multiplication D = A.B.C*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dotThreeTensors(const classTensor2& A, const classTensor2& B, const classTensor2& C, classTensor2& D, double fact =1.)
{
    D.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    D(i,j) += A(i,k)*B(k,l)*C(l,j)*fact;
                }
            }
        }
    }
};

/*! \brief Multiplication D += A.B.C*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dotThreeTensorsAdd_firstTwoIndexesFourthTensor(const classTensor2& A, const classTensor4& B, const classTensor2& C, classTensor4& D, double fact =1.)
{
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    for (int p=0; p< ndim; p++)
                    {
                        for (int q=0; q< ndim; q++)
                        {
                            D(i,j,p,q) += A(i,k)*B(k,l,p,q)*C(l,j)*fact;
                        }
                    }
                }
            }
        }
    }
};

/*! \brief Multiplication D += A.B.C*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dotThreeTensors_firstTwoIndexesFourthTensor(const classTensor2& A, const classTensor4& B, const classTensor2& C, classTensor4& D, double fact =1.)
{
    D.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    for (int p=0; p< ndim; p++)
                    {
                        for (int q=0; q< ndim; q++)
                        {
                            D(i,j,p,q) += A(i,k)*B(k,l,p,q)*C(l,j)*fact;
                        }
                    }
                }
            }
        }
    }
};

/*! \brief Multiplication C = A.B*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dot(const classTensor2& A, const classTensor1& B, classTensor1& C, double fact = 1.)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            C(i) += A(i,j)*B(j)*fact;
        }
    }
};


/*! \brief Multiplication C = A:B
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void doubleDot(const classTensor4& A, const classTensor2& B, classTensor2& C, double fact = 1.)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    C(i,j) += A(i,j,k,l)*B(k,l)*fact;
                }
            }
        }
    }
};

/*! \brief Multiplication C = A:B
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void doubleDot(const classTensor2& A, const classTensor4& B, classTensor2& C, double fact=1.)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                for (int l=0; l< ndim; l++)
                {
                    C(k,l) += A(i,j)*B(i,j,k,l)*fact;
                }
            }
        }
    }
};


/*! \brief Multiplication C = FT.F
  @param[in] A input
  @param[out] C output
  @param[out] DCDA output if this pointer is non-null
 */
inline void rightCauchy(const classTensor2& F, classTensor2& C, classTensor4* DCDF=NULL)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                C(i,j) += F(k,i)*F(k,j);
            }
        }
    }
    if (DCDF != NULL)
    {
        static classTensor2 I(ndim, 1.);
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< ndim; j++)
            {
                 for (int p=0; p< ndim; p++)
                {
                    for (int q=0; q< ndim; q++)
                    {
                        (*DCDF)(i,j,p,q) = F(p,j)*I(i,q) + F(p,i)*I(j,q);
                    }
                }
            }
        }
                
    }
};

/*! \brief Multiplication C = F.FT
  @param[in] A input
  @param[out] C output
  @param[out] DCDA output if this pointer is non-null
 */
inline void leftCauchy(const classTensor2& F, classTensor2& C, classTensor4* DCDF=NULL)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                C(i,j) += F(i,k)*F(j,k);
            }
        }
    }
    if (DCDF != NULL)
    {
        static classTensor2 I(ndim, 1.);
        for (int i=0; i< ndim; i++)
        {
            for (int j=0; j< ndim; j++)
            {
                 for (int p=0; p< ndim; p++)
                {
                    for (int q=0; q< ndim; q++)
                    {
                        (*DCDF)(i,j,p,q) = F(j,q)*I(i,p) + F(i,q)*I(j,p);
                    }
                }
            }
        }
                
    }
};

/*! \brief Multiplication C = AT.B*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dotFirstTranspose(const classTensor2& A, const classTensor2& B, classTensor2& C, double fact = 1.)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                C(i,j) += A(k,i)*B(k,j)*fact;
            }
        }
    }
};

/*! \brief Multiplication C = A.BT*fact
  @param[in] A input
  @param[in] B input
  @param[out] C output
 */
inline void dotSecondTranspose(const classTensor2& A, const classTensor2& B, classTensor2& C, double fact =1.)
{
    C.setAll(0);
    int ndim = C.getDim();
    for (int i=0; i< ndim; i++)
    {
        for (int j=0; j< ndim; j++)
        {
            for (int k=0; k< ndim; k++)
            {
                C(i,j) += A(i,k)*B(j,k)*fact;
            }
        }
    }
};


inline bool polynomial(const classTensor2& a, const int order, const vector<double>& coeffs, 
                classTensor2& F, classTensor4* dF, classTensor6* ddF){
    if (order+1 != coeffs.size()) 
    {
        printVector(coeffs, "coeffs");
        ERROR("approximation is wrongly defined with the order = %d",order);
        exit(-1);
    };
    int ndim = a.getDim();
    const classTensor2 I(ndim,1.);
    const classTensor4 I4(ndim,1.,true);

    if (a.isnan())
    {
        WARNING("polynomial of a NAN tensor cannot be performed");
        F.setAll(0);
        if (dF)
        {
            (*dF) = I4;
            if (ddF)
            {
                ddF->setAll(0);
            }
        }
        return false;
    }
    
    classTensor2 Bn(ndim);
    classTensor4 dBnDa(ndim);
    classTensor6 ddBnDDa(ndim,ndim);
      
    for (int j=0; j<ndim; j++)
    {
        for (int k=0; k<ndim; k++)
        {
          Bn(j,k) = I(j,k)*coeffs[order];
        }
    }
    if (dF != NULL)
    {
        dBnDa.setAll(0);
        if (ddF!= NULL)
        {
          ddBnDDa.setAll(0);
        }
    }
  
    for (int i=0; i< order; i++)
    {
        F.setAll(0);
        if (dF != NULL)
        {
            dF->setAll(0);
            if (ddF != NULL)
            {
                ddF->setAll(0);
            }
        }
        
        int n =order-i;
        for (int j=0; j<ndim; j++)
        {
            for (int k=0; k<ndim; k++)
            {
                F(j,k) = coeffs[n-1]*I(j,k);
                for (int l=0; l<ndim; l++)
                {
                    F(j,k) +=  Bn(j,l)*a(l,k);
                    if (dF != NULL)
                    {
                        for (int r=0; r<ndim; r++)
                        {
                            for (int s=0; s<ndim; s++)
                            {
                                (*dF)(j,k,r,s) +=  dBnDa(j,l,r,s)*a(l,k) + Bn(j,l)*I4(l,k,r,s);
                                if (ddF != NULL)
                                {
                                    for (int p=0; p<ndim; p++)
                                    {
                                        for (int q=0; q<ndim; q++)
                                        {
                                            (*ddF)(j,k,r,s,p,q) +=  ddBnDDa(j,l,r,s,p,q)*a(l,k)+dBnDa(j,l,r,s)*I4(l,k,p,q) + dBnDa(j,l,p,q)*I4(l,k,r,s);
                                        }
                                    }
                                }
                            }
                        }              
                    }
                }
            }
        }
        Bn = F;
        if (dF != NULL)
        {
            dBnDa = (*dF);
            if (ddF != NULL)
            {
                ddBnDDa = (*ddF);
            }
        }
    }
    return true;
};

inline void logClassTensor2(const classTensor2& a, const int order, classTensor2& loga, classTensor4* dloga, classTensor6* ddloga)
{
    int ndim = a.getDim();
    const classTensor2 I(ndim,1.);
    const classTensor4 I4(ndim,1.,true);
    
    classTensor2 ami(ndim,-1.);
    ami.axpy(a,1.);

    // Choose approximation order
    if (order == 1)
    {
        // Order 1 : log(1+A) = A
        loga = ami;
        // First derivative
        if(dloga!=NULL)
        {
          (*dloga) = I4;
        }
        // Second derivative
        if(ddloga!=NULL)
        {
            ddloga->setAll(0);            
        }
    }
    else if (order > 1)
    {
        vector<double> coeffs(order+1);
        // logarithmic approximation
        for (int i=0; i< order+1; i++)
        {
            if (i==0)
            {
                coeffs[i] = 0.;
            }
            else if (i%2 == 0)
            {
                coeffs[i] = -1./double(i);
            }
            else
            {
                coeffs[i] = 1./double(i);
            }
        }
    
        bool ok =  polynomial(ami,order,coeffs,loga,dloga,ddloga);
        if (!ok)
        {
            WARNING("The computation of the logClassTensor2 is not successful");
        }
    }
}


inline void expClassTensor2(const classTensor2& a, const int order, classTensor2& expa, classTensor4* dexpa, classTensor6* ddexpa)
{
    
    int ndim = a.getDim();
    const classTensor2 I(ndim,1.);
    const classTensor4 I4(ndim,1.,true);
  
    // Choose approximation order 
    if (order == 1) 
    {
        // Order 1 : exp(A) = I + A
        expa = I;
        expa.axpy(a,1);
        // First derivative
        if(dexpa!=NULL)
        {
            (*dexpa) = I4;
        }
        // Second derivative
        if(ddexpa!=NULL)
        {
            ddexpa->setAll(0);
        }
    }
    else if (order > 1)
    {
        vector<double> coeffs(order+1);
        coeffs[0] =1.;
        coeffs[1] =1;
        for (int i=2; i< order+1; i++)
        {
            coeffs[i] = coeffs[i-1]/double(i);
        }
        
    
        bool ok = polynomial(a,order,coeffs,expa,dexpa,ddexpa);
        if (!ok)
        {
            WARNING("The computation of the expClassTensor2 is not successful");
        }
    }
}
};
#endif