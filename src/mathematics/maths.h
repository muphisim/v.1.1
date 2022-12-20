//
//
// File author(s): see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _maths_H_
#define _maths_H_

#include "configuration.h"
#include <petscsnes.h>
#include <petscksp.h>
#include <map>
#include <set>


/*! \brief This class is to manage the large matrices, i.e., KK and M
*/
class classMatrix {

protected:
    int _r, _c, _size;
    string _type;

public:
    classMatrix(int r, int c, string type);
    virtual void preallocate()=0;
    virtual bool isPreallocated() = 0;
    virtual void insertSparsityPattern(int i, int j) = 0;

    // val should be ind1.size()*ind2.size()
    virtual void setValues(vector<int> ind1, vector<int> ind2, vector<double> val) = 0;

    virtual int getRowNum();

    virtual int getColNum();

    virtual void setIJ(int i, int j, double val) = 0;

    virtual double getIJ(int i, int j) = 0;

    virtual void addIJ(int i, int j, double val) = 0;

    virtual void assembly() = 0;

    virtual vector<double> getRow(int row) = 0;

    virtual vector<double> getCol(int col) = 0;

    virtual void getPetscMat(Mat &A) = 0;// All type of matrices should have a getpetsc if you want to use petsc solvers

    virtual void getPetscMatSeq(Mat &A) = 0;// Used when the computation is on sequential when solving for a linear system in DivisionRandomVariablesRandomVariables
    
    virtual void zeroEntries() = 0;

    inline int getNumRows() { return _r; };

    inline int getSize() { return _size; };

    inline string getType() { return _type; };

    //virtual void multFactor(double factor)=0;
    virtual ~classMatrix() {};
};


/*! \brief This class is to manage fourth order tensors
*/
class classCompactMatrix : public classMatrix {

protected:
    vector<double> _data;
    Mat _Petsc;
public:
    classCompactMatrix(int r, int c, string type);
    
    virtual void preallocate(){};
    virtual bool isPreallocated() {return true;};
    virtual void insertSparsityPattern(int i, int j){}

    virtual void setValues(vector<int> ind1, vector<int> ind2, vector<double> val);

    virtual double getIJ(int i, int j);

    virtual void setIJ(int i, int j, double val);

    virtual void addIJ(int i, int j, double val);

    virtual vector<double> getRow(int row);

    virtual vector<double> getCol(int col);

    virtual void assembly();

    virtual void getPetscMat(Mat &A);

    virtual void getPetscMatSeq(Mat &A);

    virtual void zeroEntries();

///virtual void multFactor(double factor);
    ~classCompactMatrix() {
        MatDestroy(&_Petsc);
    };
};

/*! \brief This class is to manage the sparsity of a CRS matrix 
*/
class sparsityPattern
{
  protected:
    std::map<int, std::set<int> > _allocatedPosition; // row as key,
    int _numRows;
#ifdef PARALLEL        
    vector<int> _numRowsRanks;
#endif 

  public:
    sparsityPattern(int r): _numRows(r)
    {

    }
    ~sparsityPattern(){}
    int getNumSettingRows() {return _allocatedPosition.size();};
    void insert(int i, int j);
#ifdef PARALLEL    
    void setNumRowsRanks(const vector<int>& vv) {_numRowsRanks=vv;}
    void getRowPattern(std::vector<int>& nByRow, std::vector<int>& nByOffRow);
#else
   void getRowPattern(std::vector<int>& nByRow);
#endif 
};


/*! \brief This class is to manage fourth order tensors
*/
class classCRSMatrix : public classMatrix {

protected:
    Mat _data;
    vector<int> _indices;
    sparsityPattern _sparPattern;
    bool _isPreallocated;
    
public:
    classCRSMatrix(int r, int c, string type);
    
    virtual void preallocate();
    virtual bool isPreallocated() {return _isPreallocated;};
    virtual void insertSparsityPattern(int i, int j);
    
    virtual void setValues(vector<int> ind1, vector<int> ind2, vector<double> val);

    virtual double getIJ(int i, int j);

    virtual void setIJ(int i, int j, double val);

    virtual void addIJ(int i, int j, double val);

    virtual vector<double> getRow(int row);

    virtual vector<double> getCol(int col);

    virtual void assembly();

    virtual void getPetscMat(Mat &A);

    virtual void getPetscMatSeq(Mat &A);

    virtual void zeroEntries();
//virtual void multFactor(double factor);

    ~classCRSMatrix() {
        MatDestroy(&_data);
    };
};


/*! \brief This class is to manage fourth order tensors
*/
class classLargeMatrix {

protected:
    int size;
    vector<double> values;
    vector<int> rows;
    vector<int> columns;

public:
    classLargeMatrix(int size);

    inline void add(double value, int row, int column) {
        this->values.push_back(value);
        this->rows.push_back(row);
        this->columns.push_back(column);
    };

    inline void printAll() {
        INFO("The matrix is of size %ix%i. All other values are zero", size, size);
        for (int i = 0; i < this->values.size(); i++) {
            INFO("%i%i = %f", rows[i], columns[i], values[i]);
        }
    };

    void getCRS(vector<double> &val, vector<int> &col, vector<int> &rowPtr);

    void getArrays(vector<double> &val, vector<int> &col, vector<int> &row);

    inline void clear() {
        values.clear();
        columns.clear();
        rows.clear();
    };

    ~classLargeMatrix() {

    };
};

/*! \brief This class is to manage firts orde tensors
*/
class classTensor1 {

    protected:
        vector<double> values;
        int ndim;

    public:
        classTensor1(int ndim=0);
        classTensor1(int ndim, const vector<double>& vec);
        classTensor1(const classTensor1& src): values(src.values),ndim(src.ndim){}
        classTensor1& operator=(const classTensor1& src)
        {
            values = src.values;
            ndim = src.ndim;
            return *this;
        }
        int getDim() const {return ndim;};
        vector<double>& getData() {return values;};
        const vector<double>& getData() const {return values;};
        
        void setAll(double s);
        inline void scale(double s)
        {
            for (int i=0; i< values.size(); i++)
            {
                values[i] *= s;
            }
        }
        double operator()(int i) const 
        {
            return values[i];
        }
        double& operator()(int i) 
        {
            return values[i];
        }

        inline double get(int i) const {
            return values[i];
        }

        inline void setValues(int i, double tmp) {
            values[i] = tmp;
        }
        inline void addValues(int i,double tmp) {
            values[i] += tmp;
        }

        ~classTensor1() {

        };
};

/*! \brief This class is to manage second order tensors
*/
class classTensor2 {

    protected:
        vector<double> values;
        int ndim;

    public:
        classTensor2(int ndim=0, double v=0);
        classTensor2(int ndim, const vector<double>& vec);
        classTensor2(const classTensor2& src): values(src.values),ndim(src.ndim){}
        classTensor2& operator=(const classTensor2& src)
        {
            values = src.values;
            ndim = src.ndim;
            return *this;
        }
        int getDim() const {return ndim;};
        bool isnan() const
        {
            for (int i=0; i<values.size(); i++){
                if (values[i] != values[i]) return true;
            }
            return false;
        }
        double trace() const
        {
            double tr= 0;
            for (int i=0; i< ndim; i++)
            {
                tr += (*this)(i,i);
            }
            return tr;
        }
        void dev(classTensor2& dev) const
        {
            double tr = trace();
            dev = (*this);
            for (int i=0; i< ndim; i++)
            {
                dev(i,i) -= (tr/3);
            }
        }
        vector<double>& getData() {return values;};
        const vector<double>& getData() const {return values;};
        void setAll(double s);
        inline void scale(double s)
        {
            for (int i=0; i< values.size(); i++)
            {
                values[i] *= s;
            }
        }
        double operator()(int i, int j) const 
        {
            return values[i*ndim+j];
        }
        double& operator()(int i, int j) 
        {
            return values[i*ndim+j];
        }

        inline double get(int i, int j) const {
            return values[i*ndim+j];
        }

        inline void setValues(int i, int j, double tmp) {
            values[i*ndim+j] = tmp;
        }
        inline void addValues(int i, int j, double tmp) {
            values[i*ndim+j] += tmp;
        }
        inline void axpy(const classTensor2& src, double alpha)
        {
            // this = this+ alpha*src
            for (int i=0; i< ndim; i++)
            {
                for (int j=0; j < ndim; j++)
                {
                    (*this)(i,j) += alpha*src(i,j);
                }
            }
        }
        inline void addDiagonal(double v)
        {
            for (int i=0; i< ndim; i++)
            {
                 (*this)(i,i) += v;
            }
        }
        
        inline void printData(string name) const
        {
            INFO("classTensor2: %s:",name.c_str(),values.size());
            for (int i=0; i< ndim; i++)
            {
                for (int j=0; j < ndim; j++)
                {
                    printf(" %g",(*this)(i,j));
                }
                printf("\n");
            }
            printf("\n");
            
        }
        ~classTensor2() {

        };
};

/*! \brief This class is to manage second order tensors
*/
class classTensor3 {

    protected:
        vector<double> values;
        int ndim;

    public:
        classTensor3(int ndim);
        classTensor3(const classTensor3& src): values(src.values),ndim(src.ndim){}
        classTensor3& operator=(const classTensor3& src)
        {
            values = src.values;
            ndim = src.ndim;
            return *this;
        }
        int getDim() const {return ndim;};
        vector<double>& getData() {return values;};
        const vector<double>& getData() const {return values;};
        void setAll(double s);
        inline void scale(double s)
        {
            for (int i=0; i< values.size(); i++)
            {
                values[i] *= s;
            }
        }
        double operator()(int i, int j, int k) const 
        {
            return values[i*ndim*ndim+j*ndim+k];
        }
        double& operator()(int i, int j, int k) 
        {
            return values[i*ndim*ndim+j*ndim+k];
        }

        inline double get(int i, int j, int k) const {
            return values[i*ndim*ndim+j*ndim+k];
        }

        inline void setValues(int i, int j, int k, double tmp) {
            values[i*ndim*ndim+j*ndim+k] = tmp;
        }
        inline void addValues(int i, int j, int k, double tmp) {
            values[i*ndim*ndim+j*ndim+k] += tmp;
        }

        ~classTensor3() {

        };
};

/*! \brief This class is to manage fourth order tensors
*/
class classTensor4 {

protected:
    vector<double> values;
    int ndim;
    static int ind3d[3][3][3][3];
    static int ind2d[2][2][2][2];

public:
    classTensor4(int ndim);
    classTensor4(int ndim, double v1, bool sym); // identity tensor
    classTensor4(const classTensor4& src): values(src.values),ndim(src.ndim){}
    classTensor4& operator=(const classTensor4& src)
    {
        values = src.values;
        ndim = src.ndim;
        return *this;
    }
    int getDim() const {return ndim;};
    inline void scale(double s)
    {
        for (int i=0; i< values.size(); i++)
        {
            values[i] *= s;
        }
    }
    inline void axpy(const classTensor4& src, double alp)
    {
        for (int i=0; i< values.size(); i++)
        {
            values[i] += src.values[i]*alp;
        }
    }
    inline double get(int i, int j, int k, int l) const {
        if (ndim == 2) {
            return values[ind2d[i][j][k][l]];
        } else {
            return values[ind3d[i][j][k][l]];
        }
    }

    inline void setValues(int i, int j, int k, int l, double tmp) {
        if (ndim == 2) {
            values[ind2d[i][j][k][l]] = tmp;
        } else {
            values[ind3d[i][j][k][l]] = tmp;
        }
    }
    vector<double>& getData() {return values;};
    const vector<double>& getData() const {return values;};
    void setAll(double s);
    double operator()(int i, int j, int k, int l) const 
    {
        if (ndim == 2) {
            return values[ind2d[i][j][k][l]];
        } else {
            return values[ind3d[i][j][k][l]];
        }
    }
    double& operator()(int i, int j, int k, int l) 
    {
        if (ndim == 2) {
            return values[ind2d[i][j][k][l]];
        } else {
            return values[ind3d[i][j][k][l]];
        }
    }

    ~classTensor4() {

    };
};

/*! \brief This class is to manage fourth order tensors
*/
class classTensor6 {

protected:
    int ndim;
    int nstoch;
    vector<double> values;

public:
    classTensor6(int ndim, int nstoch);
    classTensor6(const classTensor6& src): ndim(src.ndim),nstoch(src.nstoch),values(src.values){}
    classTensor6& operator=(const classTensor6& src)
    {
        ndim = src.ndim;
        nstoch = src.nstoch;
        values = src.values;
        return *this;
    }
    inline double get(int n, int m, int i, int j, int k, int l) const {
        return values[i + ndim * j + ndim * ndim * k + ndim * ndim * ndim * l + ndim * ndim * ndim * ndim * m + ndim * ndim * ndim * ndim * nstoch * n];

    }
    inline double operator()(int n, int m, int i, int j, int k, int l) const
    {
        return values[i + ndim * j + ndim * ndim * k + ndim * ndim * ndim * l + ndim * ndim * ndim * ndim * m + ndim * ndim * ndim * ndim * nstoch * n];
    }
    inline double& operator()(int n, int m, int i, int j, int k, int l)
    {
        return values[i + ndim * j + ndim * ndim * k + ndim * ndim * ndim * l + ndim * ndim * ndim * ndim * m + ndim * ndim * ndim * ndim * nstoch * n];
    }

    inline void setValues(int n, int m, int i, int j, int k, int l, double tmp) {
        values[i + ndim * j + ndim * ndim * k + ndim * ndim * ndim * l + ndim * ndim * ndim * ndim * m + ndim * ndim * ndim * ndim * nstoch * n] = tmp;
    }
    
    void setAll(double s)
    {
        fill(values.begin(), values.end(), s);
    }
    int getDim() const {return ndim;};
    ~classTensor6() {

    };
};


class mVector
{
    protected:
        vector<double> values;
    
    public:
        mVector(int n=0): values(n,0){}
        mVector(const mVector& src): values(src.values){}
        mVector& operator = (const mVector& src)
        {
            values = src.values;
            return *this;
        }
        ~mVector(){}
        const vector<double>& getData() const {return values;};
        vector<double>& getData() {return values;};
        void setAll(double s)
        {
            fill(values.begin(),values.end(),s);
        }
        void resize(int n, bool resetValue)
        {
            values.resize(n);
            if (resetValue)
            {
                setAll(0.);
            }
        }
        int size() const {return values.size();};
        double operator()(int i) const {return values[i];};
        double& operator()(int i) {return values[i];};
        void axpy(const mVector& src, double alp)
        {
            for (int i=0; i< values.size(); i++)
            {
                values[i] += src.values[i]*alp;
            }
        }
        void printData(string name) const{
            INFO("vector %s of %ld elements:",name.c_str(),values.size());
            for (int i=0; i< values.size(); i++)
            {
                printf(" %f",values[i]);
            }
            printf("\n");
        }
};
class mMatrix
{
    protected:
        vector<double> values;
        int nrows, ncols;
        
    public:
        mMatrix(int r=0, int c=0): values(r*c,0),nrows(r),ncols(c){}
        mMatrix(const mMatrix& src): values(src.values),nrows(src.nrows),ncols(src.ncols){}
        mMatrix& operator = (const mMatrix& src)
        {
            values = src.values;
            nrows = src.nrows;
            ncols = src.ncols;
            return *this;
        }
        ~mMatrix(){}
        const vector<double>& getData() const {return values;};
        vector<double>& getData() {return values;};
        void mult(const mVector& a, mVector& res) const
        {
            if (ncols != a.size())
            {
                INFO("Matrix-vector multiplication impossible: dimension incompatibility");
                exit(-1);
            }
            res.resize(nrows,true);
            for (int i=0; i< nrows; i++)
            {
                for (int j=0; j< ncols; j++)
                {
                    res(i) += (*this)(i,j)*a(j);
                }
            }
            
        }
        void copy(const mMatrix& src, int rowStart, int colStart)
        {
            int crow = std::min(nrows,rowStart+src.sizeRows());
            int ccol = std::min(ncols,colStart+src.sizeCols());
            for (int i=rowStart; i < crow; i++)
            {
                for (int j= colStart; i< ccol; j++)
                {
                    (*this)(i,j) = src(i-rowStart,j-colStart);
                }
            }
        }
        void setAll(double s)
        {
            fill(values.begin(),values.end(),s);
        }
        void resize(int r, int c, bool resetValue)
        {
            nrows = r;
            ncols = c;
            values.resize(r*c);
            if (resetValue)
                setAll(0);
        }
        int sizeRows() const {return nrows;};
        int sizeCols() const {return ncols;};
        double operator()(int i, int j) const {return values[i*ncols+j];};
        double& operator()(int i, int j) {return values[i*ncols+j];};
        void printData(string name) const
        {
            INFO("matrix %s of %ld elements, rows %d, cols %d:",name.c_str(),values.size(),nrows,ncols);
            for (int i=0; i< nrows; i++)
            {
                for (int j=0; j< ncols; j++)
                {
                    printf(" %f",(*this)(i,j));
                }
                printf("\n");
            }
        }
        void printDiagonalData(string name) const
        {
            INFO("diagonal data : matrix %s of %ld elements, rows %d, cols %d:",name.c_str(),values.size(),nrows,ncols);
            for (int i=0; i< nrows; i++)
            {
                
                printf(" %f",(*this)(i,i));
                
            }
            printf("\n");
        }
};
    
extern void largeMatrixMultVector(classMatrix *&A, const vector<double> &mm, vector<double> &val);

extern void largeMatricesSum(double sign, classMatrix *&A, classMatrix *&B);

extern vector<double> devSTensor3(const vector<double> &A, int ndim);

extern void dyadicProduct(const vector<double> &A, int m, const vector<double> &mm, int l, vector<double> &C);

extern void printVector(const vector<double>& A, string name, bool nonZeroOnly=false);
extern void printVector(const vector<int>& A, string name);

extern void setAll(vector<double> &A, double s);
extern void setAll(vector<vector<double> > &A, double s);
extern void scale(vector<double> &A, double s);
extern void scale(vector<vector<double> > &A, double s);

extern void transpose(const vector<double> &t, int n, vector<double> &val);

extern void transposeStochastic(const vector<double> &A, int n, int nstoch, vector<double> &C);

extern void computeGL(const vector<double> &FCurr, int ndim, vector<double> &E);

extern void InverseTensorStochastic(const vector<double> &A, int ndim, int nstoch,
                                    const vector<double> &C, vector<double> &D);

extern void computeCauchy(const vector<double> &FCurr, int ndim, const vector<double> &PK1, vector<double> &Cauchy);

extern void computeVMS(vector<double> &Cauchy, int ndim, vector<double> &VMS);

extern void computeVolumetricStrain(vector<double> &E, int ndim, vector<double> &volumetricStrain);

extern void computeEquivalentStrain(const vector<double> &E, int ndim, vector<double> &equivalentStrain);

extern void computeGLstoch(const vector<double> &FCurr, int ndim, int nstoch, const vector<double>&C_Stoch, vector<double> &E);

extern void computeCauchyStoch(const vector<double> &FCurr, int ndim, int nstoch, const vector<double> &C_stoch, const vector<double> &PK1, vector<double> &Cauchy);

extern void computeStochasticVMS(vector<double> &Cauchy, int ndim, int nstoch, vector<double> &VMS, const vector<double> &C_stoch);

extern void computeVolumetricStrainStochastic(const vector<double> &E, int ndim, int nstoch, vector<double> &volumetricStrain);

extern void computeEquivalentStrainStochastic(const vector<double> &E, int ndim, int nstoch, const vector<double> &C_stoch, vector<double> &equivalentStrain);

extern double weightedAverage(const vector<double> &a, int size);

extern double norm2(const vector<double> &a, int size);

extern double determinantTensor(const vector<double> &A, int size);

extern void multTensor4Tensor(const classTensor4 *AA, vector<double> &B, vector<double> &C, int ndim);

extern void sumTensorTensor(const vector<double> &t, const vector<double> &m, vector<double> &val, int ndim, int sum);

extern void multTensorVector(const vector<double> &A, int m, int n, const vector<double> &mm, int l, vector<double> &val);

extern void multTensorVectorStochastic(const vector<double> &A, int ndim, int nstoch, const vector<double> &B, const vector<double> &C, vector<double> &D);

extern void transposeGeneralMatrix(const vector<double> &A, int m, int n, vector<double> &C);

extern void multGeneralMatrices(const vector<double> &A, int m, int n, const vector<double> &B, int kk, int ll, vector<double> &C);

extern void multTensorTensor3(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C);

extern void multSTensor3FirstTranspose(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C);

extern void multSTensor3SecondTranspose(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C);

extern void multSTensor3FirstTransposeStochastic(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C,
                                     int nstoch, const vector<double>& multiplier);

extern void multSTensor3SecondTransposeStochastic(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C,
                                     int nstoch, const vector<double>& multiplier);

extern void multTensorTensor3Stochastic(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C,
                                        int nstoch, const vector<double>& multiplier);

extern void multiSTensorScalar(vector<double> &A, int m, int n, double s);

extern void multiSTensorScalarStochastic(vector<double> &A, int m, int n, double s, int nstoch);

extern void multiSTensorRandomVariableStochastic(const vector<double> &A, int m, int n, const vector<double>& s, int nstoch,
                                                 const vector<double>& C, vector<double> &D);

extern void multiSVectorRandomVariableStochastic(const vector<double> &A, int m,const vector<double>& s, int nstoch,
                                                 const vector<double>& C, vector<double> &D);

extern void lumpMass(const vector<double> &A, int numDof, vector<double> &C);

extern void inverse(const vector<double> &A, int nn, vector<double> &C);


extern bool eigen(vector<double>&  F, int ndim, vector<double> &WR, vector<double> &WI, vector<double> &VL,
                  vector<double> &VR);

extern void dyadicProductSecondTensorSecondTensor(const vector<double> &t, int a, int b, const vector<double> &mm, int c, int d,
                                                  classTensor4 *AA);

extern void fourthTensorTranspose(int ndim, const classTensor4 *AA, classTensor4 *AATranspose);

extern void multFourthTensorFourthTensor(int ndim, const classTensor4 *AA, const classTensor4 *BB, classTensor4 *result);

extern void tensorProductSecondTensorSecondTensor(const vector<double> &t, int a, int b, const vector<double> &mm, int c, int d,
                                                  classTensor4 *AA);

extern double doubleContractionSecondTensors(const vector<double> &A, const vector<double> &B, int m);

extern vector<double> DEVSTensor3(const vector<double> &A, int ndim, const vector<double> &C, const vector<double> &Cinv);

extern void sumFourthTensorFourthTensor(int ndim, const classTensor4 *AA, const classTensor4 *BB, classTensor4 *result);

extern void ConvertToFortran(char *fstring, std::size_t fstring_len, const char *cstring);

extern void ConvertToFortran2Darray(double *array, int files, int columns, double *farray);

vector<int> local2global(const vector<vector<int> >& table_poly, int j);

double getTrace(const vector<double> &A, int ndim);

vector<double> symSTensor3(const vector<double> &a, int ndim);

vector<double> vectorOperator(double sign, const vector<double>& a, const vector<double>& b);

extern void DeterminantStochasticTensor(const vector<double> &A, int ndim, int nstoch, vector<double> &D, const vector<double> &C);

int indexofSmallestElement(const vector<double>& array);

extern void getCofactor(const vector<vector<double> > &A, vector<vector<double> > &temp, int p, int q, int n);

extern double determinant(const vector<vector<double> > &A, int n, int N);

extern int factorial(int n);

extern void isInCircle(bool &A, const vector<double>& n0, const vector<double>& center, const vector<double>& coordinates, double radius);

extern double distance (const vector<double>& center, const vector<double>& coordinates);

double QuantileNormal(double p);

double QuantileUniform(double p);

double evaluate_haar(int j, int k, double y);

int binomialCoefficients(int n, int k);

int evaluate_wavelet(double y);

extern void squareRootIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C);

extern void logIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C);

extern void cubicRootIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C);

extern void nthPowerIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C, const double &n);

extern void expIntegral(const vector<double> &A, vector<double> &Result,  const vector<double>& C);

extern void sigmoid(const vector<double> &A, vector<double> &Result, const vector<double>& C);

vector<double> linspace(double a, double b, int num);

extern void build_ThirdOrderTensor_C(int order, int number_stochastic, vector<double> &C, string distribution);

extern void DL_sqrt_1minusx(vector<double> &A, vector<double> C, vector<double> &Result);

extern void InterpolationSqrt(vector<double> &A, vector<double> &Result, vector<double> C);

extern void polyfitsqrt(vector<double> &A, vector<double> C, vector<double> &Result);

extern void matlabfit(double r, vector<double> c, vector<double> C, vector<double> &Result);

extern void piecewisefit(double r, vector<double> mu, vector<double> C, vector<double> &Result);

extern vector<double> NormalToHaar(int order, string distribution);

extern void ConvertToHaar(vector<double> &A, vector<double> &B, int order, const vector<double>& Stochastic_parameters,
                          const vector<double>& Stochastic_numbers, string distribution);

extern void Build_C_Haar(int order, int number_stochastic, vector<double> &C);

extern void Build_C_Mixed(int order, int resolution, int number_func, vector<double> &C, string distribution);

double evaluate_mixed(int i, int order, double gauss);

extern void adjoint(const vector<vector<double> >  &A, vector<vector<double> > &adj, int N);

extern bool inverse(const vector<vector<double> > &A, vector<vector<double> > &inverse, int N);

extern void vectorArithmetic(const vector<double> &a, const vector<double> &b, vector<double> &c, int ndim, bool add);

extern double dotProduct(const vector<double> &a, const vector<double> &b, int ndim);

extern vector<int> get_position(const vector<int>& box);

extern void DivisionRandomVariableRandomVariableStochastic(const vector<double> &A, const vector<double> &B, const vector<double>& C,
                                                           vector<double> &D);

extern void MultiplicationRandomVariableScalar(vector<double> &A, double scalar);

extern void linearCombinationVariableStochastic(vector<double> &A, double a, double b);

extern void sum_nonnegative_integers_k(int n, int k, vector<vector<int>> &double_poly, int &row_index);

extern bool is_scalar(const vector<double> &A);

extern bool balls_right(const vector<int>& box, int n);

extern void
multiRandomVariableRandomVariableStochastic(const vector<double> &A, const vector<double> &B, const vector<double>& C, vector<double> &D);

extern vector<double> product(const vector<double>& A, const vector<double>& B);

extern int find_index_shifted(const vector<int>& box, int n);

void linearSystemSequential(classMatrix *&KK2, const vector<double> &rhs, vector<double> &xx2, bool &flagConver);

int calculationFactorial(int i);

double calculationMomentUniform(int i);

void findCommonInterval(int i, int j, int k, bool &flagCommonInterval, vector<int> &indices, vector<double> &interval, int &max);
#endif
