//
//
// File authors: see Authors.txt
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file maths.cpp
  \brief This file contains all functions related to mathematical operations. Some of the following functions are just functions to prepare the data to call Lapack and Blas libraries
*/

#include "maths.h"
#include "commandLine.h"

classMatrix::classMatrix(int r, int c, string type) {
    _r = r;
    _c = c;
    _size = _r * _c;
    _type = type;
};

int classMatrix::getRowNum() {
    return _r;
};

int classMatrix::getColNum() {
    return _c;
};

classCompactMatrix::classCompactMatrix(int r, int c, string type) : classMatrix(r, c, type) {
    _data.resize(_size, 0.);
    _Petsc = NULL;
};

double classCompactMatrix::getIJ(int i, int j) {
    return _data[i * _c + j];
};

void classCompactMatrix::setIJ(int i, int j, double val) {
    _data[i * _c + j] = val;
};

void classCompactMatrix::addIJ(int i, int j, double val) {
    _data[i * _c + j] += val;
};

void classCompactMatrix::assembly() {
};

vector<double> classCompactMatrix::getRow(int row) {
    vector<double> val;
    val.assign(&_data[row * _c], &_data[row * _r + _c]);
    return val;
};

vector<double> classCompactMatrix::getCol(int col) {
    vector<double> val(_r, 0);
    for (int i = 0; i < _c; i++) {
        val[i] = _data[i * _r + col];
    }
    return val;
};

void classCompactMatrix::setValues(vector<int> ind1, vector<int> ind2, vector<double> val) {
    int ind1S = ind1.size();
    int ind2S = ind2.size();
    for (int i = 0; i < ind1S; i++) {
        int k = ind1[i];
        for (int j = 0; j < ind2S; j++) {
            int l = ind2[j];
            _data[_c * k + l] = val[j + i * ind2S];
            _data[_c * l + k] = val[j + i * ind2S];
        }
    }
};

void classCompactMatrix::getPetscMat(Mat &A) {

    MatCreate(PETSC_COMM_WORLD, &_Petsc);
    MatSetSizes(_Petsc, PETSC_DECIDE, PETSC_DECIDE, _r, _c);
    MatSetFromOptions(_Petsc);
    MatSetType(_Petsc, MATAIJ);
    //
    std::vector<int> nByRow(_r,_c);
    MatSeqAIJSetPreallocation(_Petsc, 0, &nByRow[0]);
    //
    MatSetUp(_Petsc);
    double small = 1E-15;
    for (int i = 0; i < _r; i++) {
        for (int j = 0; j < _c; j++) {
            if (fabs(_data[i * _c + j]) > small) {
                MatSetValues(_Petsc, 1, &i, 1, &j, &_data[i * _c + j], INSERT_VALUES);
            }
        }
    }
    MatAssemblyBegin(_Petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_Petsc, MAT_FINAL_ASSEMBLY);
    A = _Petsc;
};

void classCompactMatrix::getPetscMatSeq(Mat &A) {
    MatCreate(PETSC_COMM_SELF, &_Petsc);//Sequential so the communication is not the same.
    MatSetSizes(_Petsc, PETSC_DECIDE, PETSC_DECIDE, _r, _c);
    MatSetFromOptions(_Petsc);
    MatSetType(_Petsc, MATAIJ);
    //
    std::vector<int> nByRow(_r,_c);
    MatSeqAIJSetPreallocation(_Petsc, 0, &nByRow[0]);
    //
    MatSetUp(_Petsc);
    double small = 1E-15;
    for (int i = 0; i < _r; i++) {
        for (int j = 0; j < _c; j++) {
            if (fabs(_data[i * _c + j]) > small) {
                MatSetValues(_Petsc, 1, &i, 1, &j, &_data[i * _c + j], INSERT_VALUES);
            }
        }
    }
    MatAssemblyBegin(_Petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_Petsc, MAT_FINAL_ASSEMBLY);
    A = _Petsc;
};

void classCompactMatrix::zeroEntries() {
    fill(_data.begin(), _data.end(), 0.);
    MatDestroy(&_Petsc);
};


void sparsityPattern::insert(int i, int j)
{
  std::set<int>& rowi = _allocatedPosition[i];
  int sizeCur = rowi.size();
  rowi.insert(j);
};

#ifdef PARALLEL    
void sparsityPattern::getRowPattern(std::vector<int>& nByRow, std::vector<int>& nByOffRow)
{
    // make a same map for all procs
    int rank = GeneralOptions::commRank;
    int nprcs = GeneralOptions::commSize;
    //
    std::vector<int> localBuffer(2*_numRows,0);
    for (std::map<int, std::set<int> >::const_iterator it = _allocatedPosition.begin(); it != _allocatedPosition.end(); it++)
    {
        int row= it->first;
        const std::set<int>& allCols = it-> second;
        
        int localStart = 0;
        int localEnd = _numRowsRanks[0];
        bool found = false;
        if ((localStart <= row) and (row < localEnd))
        {
            found = true;
        }
        else
        {
            for (int i=0; i< nprcs-1; i++)
            {
                localStart += _numRowsRanks[i];
                localEnd += _numRowsRanks[i+1];
                if ((localStart <= row) and (row < localEnd))
                {
                    found = true;
                    break;
                }
            }
            
        }
        if (!found)
        {
            ERROR("row %d cannot be found in a matrix of %d rows localStart =%d localEnd=%d",row,_numRows,localStart,localEnd);
            exit(-1);
        }
        
        for (std::set<int>::const_iterator itset= allCols.begin(); itset !=allCols.end(); itset++)
        {
            if ((localStart <= *itset) and (*itset < localEnd))
            {
                localBuffer[row] ++;
            }
            else
            {
                localBuffer[row+_numRows] ++;
            }
        };
    };
    

    int Nel=localBuffer.size();
    std::vector<int> buffer(Nel,0);  
    MPI_Allreduce(&localBuffer[0], &buffer[0], Nel, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    
    int localStart = 0;
    int localEnd = _numRowsRanks[0];
    for (int i=0; i< rank; i++)
    {
        localStart += _numRowsRanks[i];
        localEnd += _numRowsRanks[i+1];
    }
    
    int localSize = localEnd - localStart;
    nByRow.resize(localSize,0);
    nByOffRow.resize(localSize,0);
  
    for (int j=0; j< localSize; j++)
    {
        nByRow[j] = buffer[j+localStart];
        nByOffRow[j] = buffer[j+localStart+_numRows];
    } 
};
#else
void sparsityPattern::getRowPattern(std::vector<int>& nByRow)
{
  nByRow.resize(_numRows);
  for (int i=0; i< _numRows; i++)
  {
    std::map<int, std::set<int> >::iterator itF = _allocatedPosition.find(i);
    if (itF ==_allocatedPosition.end())
    {
      nByRow[i] = 0;
    }
    else
    {
      nByRow[i] = itF->second.size();
    }
  }
};
#endif 


classCRSMatrix::classCRSMatrix(int r, int c, string type) : classMatrix(r, c, type), _sparPattern(r), _isPreallocated(false){
    MatCreate(PETSC_COMM_WORLD, &_data);
    MatSetSizes(_data, PETSC_DECIDE, PETSC_DECIDE, _r, _c);
    MatSetFromOptions(_data);
    MatSetType(_data, MATAIJ);
    MatSetUp(_data);

    _indices.resize(_r, 0);
    for (int k = 0; k < _r; k++) {
        _indices[k] = k;
    }
};

void classCRSMatrix::preallocate()
{
    if (_isPreallocated) return;
    _isPreallocated = true;
    
    int localVal = _sparPattern.getNumSettingRows();
    int globalVal = localVal;
#ifdef PARALLEL
   MPI_Allreduce(&localVal, &globalVal, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
#endif // PARALLEL
    if (globalVal == 0) return;    
    INFO("start pre-allocating matrix");
#ifdef PARALLEL 
    PetscInt localSize;
    MatGetLocalSize(_data,&localSize,PETSC_NULL);
    int rank = GeneralOptions::commRank;
    int nprcs = GeneralOptions::commSize;
    
    vector<int> allLocalSizes(nprcs,0);
    MPI_Allgather(&localSize, 1, MPI_INT,&allLocalSizes[0], 1, MPI_INT, PETSC_COMM_WORLD);
 
    _sparPattern.setNumRowsRanks(allLocalSizes);
    std::vector<int> nByRow, nByRowOffDiag;
    _sparPattern.getRowPattern(nByRow,nByRowOffDiag);
    MatMPIAIJSetPreallocation(_data, 0, &nByRow[0], 0, &nByRowOffDiag[0]);
#else
    std::vector<int> nByRow;
    _sparPattern.getRowPattern(nByRow);
    MatSeqAIJSetPreallocation(_data, 0, &nByRow[0]);
#endif //  
  INFO("done pre-allocating matrix");
}

double classCRSMatrix::getIJ(int i, int j) {
    double val = 0;
    MatGetValues(_data, 1, &i, 1, &j, &val);
    return val;
};

void classCRSMatrix::setIJ(int i, int j, double val) {
    preallocate();
    MatSetValues(_data, 1, &i, 1, &j, &val, INSERT_VALUES);
};

void classCRSMatrix::addIJ(int i, int j, double val) {
    preallocate();
    MatSetValues(_data, 1, &i, 1, &j, &val, ADD_VALUES);
};

void classCRSMatrix::assembly() {
    MatAssemblyBegin(_data, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_data, MAT_FINAL_ASSEMBLY);
};

vector<double> classCRSMatrix::getRow(int row) {
    vector<double> val(_r, 0);
    MatGetValues(_data, 1, &row, _r, &_indices[0], &val[0]);
    return val;
};

vector<double> classCRSMatrix::getCol(int col) {
    vector<double> val(_c, 0);
    MatGetValues(_data, _c, &_indices[0], 1, &col, &val[0]);
    return val;
};

void classCRSMatrix::insertSparsityPattern(int i, int j)
{
  _sparPattern.insert(i,j);
}

void classCRSMatrix::setValues(vector<int> ind1, vector<int> ind2, vector<double> val) {
    preallocate();
    int ind1S = ind1.size();
    int ind2S = ind2.size();
    MatSetValues(_data, ind1S, &ind1[0], ind2S, &ind2[0], &val[0], INSERT_VALUES); // Columns
};

void classCRSMatrix::getPetscMat(Mat &A) {
    A = _data;
};

void classCRSMatrix::getPetscMatSeq(Mat &A) {
    A = _data;
};

void classCRSMatrix::zeroEntries() {
    preallocate();
    MatZeroEntries(_data);
}
// void classCRSMatrix::multFactor(double factor){
//   cout<<"to timplement"<<endl;
//   exit(0);
// };


/*! \brief The matrices should be square A=sign*B + A
  @param[in] sign +/-
  @param[inout] A First tensor
  @param[in] B Second tensor
 */
extern void largeMatricesSum(double sign, classMatrix *&A, classMatrix *&B) {
    int Asize = A->getSize();
    int Bsize = B->getSize();
    if (Asize != Bsize) {
        ERROR("The sizes of the matrices should be the same");
        exit(0);
    }
    string AType = A->getType();
    string BType = B->getType();
    int n = A->getNumRows();

    if (AType == BType && BType == "CRS") {
        Mat AA;
        A->getPetscMat(AA);
        Mat BB;
        B->getPetscMat(BB);
        MatAXPY(AA, sign, BB, DIFFERENT_NONZERO_PATTERN);
        A->assembly();
    } else if (AType == BType && BType == "COMPACT") {
        // Petsc style
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A->setIJ(i, j, A->getIJ(i, j) + sign * B->getIJ(i, j));
            }
        }
        A->assembly();
    } else {
        ERROR("All matrices should be of the same type or the format is not known");
        exit(0);
    }

}

/*! \brief Multiplication of a matrix by a vector. ndim x ndim XXX ndim
  @param[in] A Second order tensor
  @param[in] m Number of rows
  @param[in] n Number of columns
  @param[in] mm Vector
  @param[in] l Size of the vector
  @param[out] val Result
 */
extern void largeMatrixMultVector(classMatrix *&A, const vector<double> &mm, vector<double> &val) {
    int Arow = A->getNumRows();
    int Bsize = mm.size();
    if (Arow != Bsize) {
        ERROR("The number of columns of the matrix and the vector size should be the same");
        exit(0);
    }
    Mat AA;
    A->getPetscMat(AA);
    PetscInt size = mm.size();
    Vec a, b;
    VecCreateSeqWithArray(PETSC_COMM_SELF,size,size,&mm[0],&a);
    VecCreateSeqWithArray(PETSC_COMM_SELF,size,size,&val[0],&b);
    MatMult(AA,a,b);
    VecDestroy(&a);
    VecDestroy(&b);
    /*for (int i = 0; i < Bsize; i++) {
        vector<double> row = A->getRow(i);
        for (int j = 0; j < Bsize; j++) {
            val[i] += row[j] * mm[j];
        }
    }*/
}


classLargeMatrix::classLargeMatrix(int size) {
    this->size = size;
}

void classLargeMatrix::getArrays(vector<double> &val, vector<int> &col, vector<int> &row) {
    val = values;
    col = columns;
    row = rows;
};

void classLargeMatrix::getCRS(vector<double> &val, vector<int> &col, vector<int> &rowPtr) {
    val = values;
    col = columns;
    for (int i = 0; i < this->size; i++) {
        int pos = distance(rows.begin(), find(rows.begin(), rows.end(), i));
        if (pos == rows.size()) {
            ERROR("The matrix is singular. One row is full of zeros");
            exit(0);
        } else {
            rowPtr.push_back(pos);
        }
    }
    rowPtr.push_back(val.size() + 1); // It should be like that
};

classTensor1::classTensor1(int ndim_): ndim(ndim_),values(ndim_, 0){
}

classTensor1::classTensor1(int ndim_, const vector<double>& vec): ndim(ndim_),values(vec)
{
    
}


void classTensor1::setAll(double s)
{
    fill(values.begin(),values.end(), s);
}



classTensor2::classTensor2(int ndim_, double v): ndim(ndim_),values(ndim_ * ndim_, 0)
{
    for (int i=0; i < ndim_; i++)
    {
        values[ndim_*i+i] = v;
    }
}

classTensor2::classTensor2(int ndim_, const vector<double>& vec): ndim(ndim_), values(&vec[0],&vec[0]+ndim_*ndim_)
{
    
}

void classTensor2::setAll(double s)
{
    fill(values.begin(),values.end(), s);
}

classTensor3::classTensor3(int ndim_): ndim(ndim_),values(ndim_* ndim_*ndim_, 0){
}

void classTensor3::setAll(double s)
{
    fill(values.begin(),values.end(), s);
}


int classTensor4::ind3d[3][3][3][3] = {{{{0,  1,  2},  {3,  4,  5},  {6,  7,  8}},  {{9,  10, 11}, {12, 13, 14}, {15, 16, 17}}, {{18, 19, 20}, {21, 22, 23}, {24, 25, 26}}},
                             {{{27, 28, 29}, {30, 31, 32}, {33, 34, 35}}, {{36, 37, 38}, {39, 40, 41}, {42, 43, 44}}, {{45, 46, 47}, {48, 49, 50}, {51, 52, 53}}},
                             {{{54, 55, 56}, {57, 58, 59}, {60, 61, 62}}, {{63, 64, 65}, {66, 67, 68}, {69, 70, 71}}, {{72, 73, 74}, {75, 76, 77}, {78, 79, 80}}}};

int classTensor4::ind2d[2][2][2][2] = {{{{0, 1}, {2,  3}},  {{4,  5},  {6,  7}}},
                         {{{8, 9}, {10, 11}}, {{12, 13}, {14, 15}}}};

classTensor4::classTensor4(int ndim_):ndim(ndim_),values(ndim_ * ndim_ * ndim_ * ndim_, 0) 
{
};

classTensor4::classTensor4(int ndim_, double v1, bool sym): ndim(ndim_),values(ndim_ * ndim_ * ndim_ * ndim_, 0)
{
    for (int i=0; i<ndim_; i++)
    {
        for (int j=0;  j< ndim_; j++)
        {
            for (int k=0; k<ndim_; k++)
            {
                for (int l=0;  l< ndim_; l++)
                {
                    if (sym)
                    {
                        if ((i == k) && (j == l))
                        {
                            if (ndim_ == 2)
                            {
                                values[ind2d[i][j][k][l]] += 0.5*v1;
                            } 
                            else 
                            {
                                values[ind3d[i][j][k][l]]+= 0.5*v1;
                            }
                        }
                        if ((i == l) && (j == k))
                        {
                            if (ndim_ == 2)
                            {
                                values[ind2d[i][j][k][l]] += 0.5*v1;
                            } 
                            else 
                            {
                                values[ind3d[i][j][k][l]]+= 0.5*v1;
                            }
                        }
                    }
                    else
                    {
                        if ((i == k) && (j == l))
                        {
                            if (ndim_ == 2)
                            {
                                values[ind2d[i][j][k][l]] += v1;
                            } 
                            else 
                            {
                                values[ind3d[i][j][k][l]] += v1;
                            }
                        }
                    }
                }
            }
        }
    }
    
}

void classTensor4::setAll(double s)
{
    fill(values.begin(),values.end(), s);
}

classTensor6::classTensor6(int _ndim, int _nstoch) : ndim (_ndim), nstoch(_nstoch),
 values(ndim * ndim * ndim * ndim * nstoch * nstoch, 0){
}

extern int factorial(int n) {
    int factorial = 1;

    if (n == 1) {
        return factorial;
    } else {
        for (int i = 1; i <= n; ++i) {
            factorial *= i;
        }
        return factorial;
    }
}

void printVector(const vector<double>& A, string name, bool nonZeroOnly)
{
    if (nonZeroOnly)
    {
        INFO("Non-zero values in vector %s of %ld elements:",name.c_str(),A.size());
        for (int i=0; i< A.size(); i++)
        {
            if (fabs(A[i]) > 0)
                printf("%d %g\n",i,A[i]);
        }
        printf("\n");
    }
    else
    {
        INFO("vector %s of %ld elements:",name.c_str(),A.size());
        for (int i=0; i< A.size(); i++)
        {
            printf(" %g",A[i]);
        }
        printf("\n");
        
    }
}
void printVector(const vector<int>& A, string name)
{
    INFO("vector %s of %ld elements:",name.c_str(),A.size());
    for (int i=0; i< A.size(); i++)
    {
        printf(" %d",A[i]);
    }
    printf("\n");
}

/*! \brief set all values in a vector*/

extern void setAll(vector<double> &A, double s)
{
    fill(A.begin(),A.end(), s);
};

extern void setAll(vector<vector<double> > &A, double s)
{
    for (int i=0; i< A.size(); i++)
    {
        setAll(A[i],s);
    }
};


extern void scale(vector<double> &A, double s)
{
    for (int i=0; i< A.size(); i++)
    {
        A[i] *= s;
    }
};

extern void scale(vector<vector<double> > &A, double s)
{
    for (int i=0; i< A.size(); i++)
    {
        scale(A[i],s);
    }
};

/*! \brief Build the third order tensor <psi_i*psi_j,psi_k>/<psi_k,psi_k> for PC (continuous)
  @param[in] order of the PC expansion
  @param[in] number_stochastic
  @param[inout] C The third order tensor
 */
extern void build_ThirdOrderTensor_C(int order, int number_stochastic, vector<double> &C, string distribution) {
    vector<double> esperance(3 * order + 1, 0);
    vector<vector<double>> poly_local(order + 1, vector<double>(3 * order + 1));
    poly_local[0][0] = 1;
    double expectation;
    if(distribution.compare(0, 8, "GAUSSIAN") == 0){
        INFO("Distributions of random variables are gaussians");
       for (int i = 0; i < esperance.size() ; i++) {
            if (i % 2 == 0) {
                expectation = calculationFactorial(i);
                esperance[i] = expectation;
            }
       }
    }

    if(distribution.compare(0, 7, "UNIFORM") == 0){
        INFO("Distributions of random variables are uniforms");
        for (int i = 0; i < esperance.size() ; i++) {
            expectation = calculationMomentUniform(i);
            esperance[i] = expectation;
        }
    }

    for (int i = 1; i <= order; i++) {
        vector<double> rhs(i + 1, 0);
        vector<double> coeff(i + 1, 0);
        rhs[i] = 1;
        vector<vector<double>> KK((i + 1), vector<double>(i + 1));
        vector<vector<double>> KKinv((i + 1), vector<double>(i + 1));
        for (int j = 0; j <= i; j++) {
            for (int k = 0; k <= i; k++) {
                if (j == i) {
                    if (k == i) {
                        KK[j][k] = 1;
                    } else {
                        KK[j][k] = 0;
                    }
                } else {
                    KK[j][k] = esperance[j + k];
                }
            }
        }
        inverse(KK, KKinv, i + 1);
        for (int l = 0; l <= i; l++) {
            coeff[l] = KKinv[l][i];
        }
        std::copy(coeff.begin(), coeff.begin() + i + 1, poly_local[i].begin());
    }
    int card = factorial(number_stochastic + order) / (factorial(number_stochastic) * factorial(order));

    vector<vector<vector<double>>> global_poly(card, vector<vector<double>>(number_stochastic,
                                                                               vector<double>(3 * order + 1)));

    for (int i = 0; i < number_stochastic; i++) {
        global_poly[0][i][0] = 1;
    }
    vector<vector<int>> table_poly(card, vector<int>(number_stochastic));

    int row_index = 1;
    for (int j = 1; j <= order; j++) {
        sum_nonnegative_integers_k(number_stochastic, j, table_poly, row_index);
    }

    for (int i = 0; i < card; i++) {
        std::reverse(table_poly[i].begin(), table_poly[i].end());
    }
    for (int k = 1; k < card; k++) {
        for (int l = 0; l < number_stochastic; l++) {
            std::copy(poly_local[table_poly[k][l]].begin(), poly_local[table_poly[k][l]].end(),
                      global_poly[k][l].begin());
        }
    }

    vector<double> product_poly;
    vector<double> product_poly_k;
    double temp = 1;
    double espe = 0;
    for (int i = 0; i < card; i++) {
        for (int j = 0; j < card; j++) {
            for (int k = 0; k < card; k++) {
                for (int l = 0; l < number_stochastic; l++) {
                    product_poly = product(global_poly[i][l], global_poly[j][l]);
                    product_poly = product(product_poly, global_poly[k][l]);
                    product_poly_k = product(global_poly[k][l], global_poly[k][l]);
                    espe = dotProduct(product_poly, esperance, 3 * order + 1) /
                           dotProduct(product_poly_k, esperance, 3 * order + 1);
                    temp = temp * espe;
                }
                C[i + j * card + card * card * k] = temp;
                temp = 1;
            }
        }

    }


}

/*! \brief computes the Haar coordinates of a random gaussian
  @param[in] order order of the PC expansion
  @param[out] Projection of the random variable on the Haar basis
 */
extern vector<double> NormalToHaar(int order, string distribution) {
    int num = 100000; // Default
    double d_eps = 1. / (num - 1);
    vector<double> cdf = linspace(0.001, 0.999, num);
    vector<double> projection(pow(2, order + 1), 0);
    for (int j = 0; j <= order; j++) {
        for (int k = 0; k <= pow(2, j) - 1; k++) {
            int index = pow(2, j) + k;
            for (int l = 0; l < num; l++) {
                if(distribution.compare(0, 8, "GAUSSIAN") == 0){
                projection[index] += QuantileNormal(cdf[l]) * evaluate_haar(j, k, cdf[l]) * d_eps;
                }
                if(distribution.compare(0, 7, "UNIFORM") == 0){
                projection[index] += QuantileUniform(cdf[l]) * evaluate_haar(j, k, cdf[l]) * d_eps;
                }
            }
        }
    }
    return projection;
}

/*! \brief Build the third order tensor <psi_i*psi_j,psi_k>/<psi_k,psi_k> for PC (continuous)
  @param[in] order order of the PC expansion
  @param[in] Stochastic_parameters list of the stochastic variables (integers)
  @param[in] Stochastic_numbers value of the variables (integers)
  @param[inout] A,B Random variables to convert from Chaos to Haar
 */
extern void ConvertToHaar(vector<double> &A, vector<double> &B, int order, const vector<double>& Stochastic_parameters,
                          const vector<double>& Stochastic_numbers, string distribution) {
    vector<double> projection = NormalToHaar(order, distribution);
    int total = 0;
    for (int i = 0; i < Stochastic_numbers.size(); i++) {
        total += Stochastic_numbers[i];
    }

    int card = binomialCoefficients(pow(2, order + 1) + total - 1, total);

    vector<vector<int>> table_poly(card, vector<int>(total));

    int row_index = 1;
    for (int j = 1; j <= pow(2, order + 1) - 1; j++) {
        sum_nonnegative_integers_k(total, j, table_poly, row_index);
    }

    vector<vector<int>> localtoglobal(total, vector<int>(pow(2, order + 1) - 1));
    vector<int> intermediate;
    for (int j = 0; j < total; j++) {
        intermediate = local2global(table_poly, j);
        for (int k = 0; k < intermediate.size(); k++) {
            localtoglobal[j][k] = intermediate[k];
        }
    }


    int count = 1;
    for (int j = 0; j < Stochastic_parameters.size(); j++) {
        if (Stochastic_parameters[j] == 1) {
            for (int k = 0; k < Stochastic_numbers[j]; k++) {
                double temp = A[count + k];
                A[count + k] = 0;
                for (int l = 0; l < localtoglobal[count - 1].size(); l++) {
                    A[localtoglobal[count - 1][l]] = projection[l + 1] * temp;
                }
                count++;
            }
        }
        if (Stochastic_parameters[j] == 2) {
            for (int k = 0; k < Stochastic_numbers[j]; k++) {
                double temp = B[count + k];
                B[count + k] = 0;
                for (int l = 0; l < localtoglobal[count - 1].size(); l++) {
                    B[localtoglobal[count - 1][l]] = projection[l + 1] * temp;
                }
                count++;
            }
        }
    }


}

/*! \brief Build the third order tensor <psi_i*psi_j,psi_k>/<psi_k,psi_k> for Haar (piecewise-continuous)
  @param[in] order order of the Haar expansion (number of levels)
  @param[in] Stochastic_numbers Total number of random variables (integers)
  @param[inout] C Third order tensor
 */

extern void Build_C_Haar(int order, int number_stochastic, vector<double> &C) {
    if (number_stochastic > 1) {
        ERROR("MuPhiSim does not allow simulation with more than one random variable with Haar Polynomials");
        exit(-1);
    }
    int num = 1000000; // Default
    vector<double> cdf = linspace(0.001, 0.999, num);

    int card = binomialCoefficients(pow(2, order + 1) + number_stochastic - 1, number_stochastic);

    vector<vector<int>> table_poly(card, vector<int>(number_stochastic));

    int row_index = 1;
    for (int j = 1; j <= pow(2, order + 1) - 1; j++) {
        sum_nonnegative_integers_k(number_stochastic, j, table_poly, row_index);
    }

    vector<vector<double>> polynomes(card, vector<double>(num));
    double temp = 1;
    for (int i = 0; i < card; i++) {
        for (int l = 0; l < num; l++) {
            for (int j = 0; j < table_poly[0].size(); j++) {
                int index = table_poly[i][j];
                if (index == 0) {
                    temp = temp;
                } else {
                    int index_1 = floor(log(index) / log(2));
                    int index_2 = index - pow(2, index_1);
                    temp *= evaluate_haar(index_1, index_2, cdf[l]);
                }
            }
            polynomes[i][l] = temp;
            temp = 1;
        }
    }
    bool flagCommonInterval;
    vector<int> indices(6, 0);
    vector<double> interval(2, 0);
    for (int i = 0; i < card; i++) {
        for (int j = 0; j < card; j++) {
            for (int k = 0; k < card; k++) {
                int max = 0;
                flagCommonInterval = false;
                findCommonInterval(i, j, k, flagCommonInterval, indices, interval, max);
                double val_1 = 0;
                double val_2 = 0;
                if (max == 0) {
                    C[i * card * card + j * card + k]=1;
                } else {
                    if (flagCommonInterval == true) {
                        int max_j = floor(log(max) / log(2));
                        int max_k = max - pow(2,max_j);
                        val_1 = pow(2, -max_j - 1);
                        val_2 = pow(2, -max_j - 1);

                        if (i == 0) {
                            val_1 = val_1;
                            val_2 = val_2;
                        } else {
                            val_1 *= evaluate_haar(indices[0], indices[1], pow(2, -max_j) * max_k);
                            val_2 *= evaluate_haar(indices[0], indices[1], pow(2, -max_j) * (max_k + 0.5));
                        }

                        if (j == 0) {
                            val_1 = val_1;
                            val_2 = val_2;
                        } else {
                            val_1 *= evaluate_haar(indices[2], indices[3], pow(2, -max_j) * max_k);
                            val_2 *= evaluate_haar(indices[2], indices[3], pow(2, -max_j) * (max_k + 0.5));


                        }

                        if (k == 0) {
                            val_1 = val_1;
                            val_2 = val_2;
                        } else {
                            val_1 *= evaluate_haar(indices[4], indices[5], pow(2, -max_j) * max_k);
                            val_2 *= evaluate_haar(indices[4], indices[5], pow(2, -max_j) * (max_k + 0.5));


                        }
                        C[i*card*card + j*card +k] = val_1 + val_2;
                    }
                }
            }
        }
    }
}

/*! \brief Build the third order tensor <psi_i*psi_j,psi_k>/<psi_k,psi_k> for Haar (piecewise-continuous)
@param[in] order order of the polynomial expansion
@param[in] resolution Depth of the approximation
@param[inout] C Third order tensor
*/
extern void Build_C_Mixed(int order, int resolution, int number_stochastic, vector<double> &C, string distribution) {


    vector<double> weight(8,0);
    weight[0] = 0.3626837833783620;
    weight[1] = 0.3626837833783620;
    weight[2] = 0.3137066458778873;
    weight[3] = 0.3137066458778873;
    weight[4] = 0.2223810344533745;
    weight[5] = 0.2223810344533745;
    weight[6] = 0.1012285362903763;
    weight[7] = 0.1012285362903763;

    vector<double>gauss(8,0);

    gauss[0] = -0.1834346424956498;
    gauss[1] = 0.1834346424956498;
    gauss[2] = -0.5255324099163290;
    gauss[3] = 0.5255324099163290;
    gauss[4] = -0.7966664774136267;
    gauss[5] = 0.7966664774136267;
    gauss[6] = -0.9602898564975363;
    gauss[7] = 0.9602898564975363;

    for(int i=0; i<8;i++){
        gauss[i] = (gauss[i]+1)/2;
    }
    int num = 100000; // Default
    double d_eps = 1. / num;
    vector<double> cdf = linspace(0.001, 0.999, num);

    int card = factorial(number_stochastic + (order+1)*(2*pow(2,resolution)) - 1 ) / (factorial(number_stochastic) * factorial((order+1)*(2*pow(2,resolution)) -1));

    int total_order = (order+1)*(2*pow(2,resolution));
    vector<vector<int>> table_poly(card, vector<int>(number_stochastic));

    int row_index = 1;
    for (int j = 1; j <= total_order -1; j++) {
        sum_nonnegative_integers_k(number_stochastic, j, table_poly, row_index);
    }

    for (int i = 0; i < card; i++) {
        std::reverse(table_poly[i].begin(), table_poly[i].end());
    }

    vector<double> stock(number_stochastic,0);


    for(int i=0; i<card; i++){
        for(int j=0; j<card; j++){
            for(int k=0; k<card; k++){
                for (int l = 0; l < number_stochastic; l++) {
                    double val =0;
                    for(int m=0; m<num;m++){
                    val+= d_eps*evaluate_mixed(table_poly[i][l], order, cdf[m])*evaluate_mixed(table_poly[j][l], order, cdf[m])*evaluate_mixed(table_poly[k][l], order, cdf[m]);
                    }
                    stock[l] = val;
                }
                double  temp = 1;
                for (int l = 0; l < number_stochastic; l++){
                    temp = temp*stock[l];
                }
            C[i*card*card + j*card + k] =temp;

            }
        }
    }

}
/*! \brief Multiplication between two tensors. Just used for square tensors, e.g., ndim x ndim XXX ndim x ndim
  @param[in] A First tensor
  @param[in] m Number of rows first tensor
  @param[in] n Number of columns first tensor
  @param[in] B Second tensor
  @param[in] l Number of rows second tensor
  @param[out] C Result
 */
extern double doubleContractionSecondTensors(const vector<double> &A, const vector<double> &B, int m) {
    double temp = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            temp += A[i * m + j] * B[i * m + j];
        }
    }
    return temp;
}

/*! \brief Dyadic product between two vectors
  @param[in] ndim Dimension of the domain
  @param[in] AA Fourth order tensor
  @param[in] BB Fourth order tensor
  @param[in] result Fourth order tensor
 */
extern void sumFourthTensorFourthTensor(int ndim, const classTensor4 *AA, const classTensor4 *BB, classTensor4 *result) {
    double temp = 0;
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                for (int l = 0; l < ndim; l++) {
                    temp = AA->get(i, j, k, l) + BB->get(i, j, k, l);
                    result->setValues(i, j, k, l, temp);
                }
            }
        }
    }

};


/*! \brief Dyadic product between two vectors
  @param[in] ndim Dimension of the domain
  @param[in] AA Fourth order tensor
  @param[in] BB Fourth order tensor
  @param[in] result Fourth order tensor
 */
extern void multFourthTensorFourthTensor(int ndim, const classTensor4 *AA, const classTensor4 *BB, classTensor4 *result) {
    double temp = 0;
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                for (int l = 0; l < ndim; l++) {
                    temp = 0;
                    for (int p = 0; p < ndim; p++) {
                        for (int q = 0; q < ndim; q++) {
                            temp += AA->get(i, j, p, q) * BB->get(p, q, k, l);
                            // cout<<"temp "<< temp<<endl;
                        }
                    }
                    result->setValues(i, j, k, l, temp);
                }
            }
        }
    }

};


/*! \brief Dyadic product between two vectors
  @param[in] ndim Dimension of the domain
  @param[in] AA Fourth order tensor
  @param[in] AATranspose Fourth order tensor
 */
extern void fourthTensorTranspose(int ndim, const classTensor4 *AA, classTensor4 *AATranspose) {
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                for (int l = 0; l < ndim; l++) {
                    AATranspose->setValues(i, j, k, l, AA->get(k, l, i, j));
                }
            }
        }
    }
}


/*! \brief Dyadic product between two vectors
  @param[in] t Second order tensor
  @param[in] a number of rows
  @param[in] b number of columns
  @param[in] mm Second order tensor
  @param[in] c number of rows of mm
  @param[in] d number of columns of mm
  @param[in] a number of rows
  @param[out] AA Fourth order tensor
 */
extern void tensorProductSecondTensorSecondTensor(const vector<double> &t, int a, int b, const vector<double> &mm, int c, int d,
                                                  classTensor4 *AA) {

    for (int i = 0; i < a; i++) {
        for (int j = 0; j < b; j++) {
            for (int k = 0; k < c; k++) {
                for (int l = 0; l < d; l++) {
                    AA->setValues(i, j, k, l, t[i * b + k] * mm[j * d + l]);
                }
            }
        }
    }
}


/*! \brief Dyadic product between two vectors
  @param[in] t First second order tensor
  @param[in] a Number of rows
  @param[in] b Number of Columns
  @param[in] mm Second order tensor
  @param[in] c Number of rows of mm
  @param[in] d Number of Columns of mm
  @param[out] AA Fourth order tensor
 */
extern void dyadicProductSecondTensorSecondTensor(const vector<double> &t, int a, int b, const vector<double> &mm, int c, int d,
                                                  classTensor4 *AA) {

    for (int i = 0; i < a; i++) {
        for (int j = 0; j < b; j++) {
            for (int k = 0; k < c; k++) {
                for (int l = 0; l < d; l++) {
                    AA->setValues(i, j, k, l, t[i * b + j] * mm[k * d + l]);
                }
            }
        }
    }
}


/*! \brief This function calculates the deviatoric part of a tensor
  @param[in] A Tensor to calculate its deviatoric part
  @param[in] ndim Dimension of the tensor
  @param[in] C Cauchy-Green tensor
  @param[in] Cinv Cauchy-Green inverse tensor
  @return The deviatoric part of the tensor A
 */
extern vector<double> DEVSTensor3(const vector<double> &A, int ndim, const vector<double> &C, const vector<double> &Cinv) {

    double contract = doubleContractionSecondTensors(A, C, ndim);
    vector<double> DEVA(ndim * ndim, 0.);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            DEVA[i * ndim + j] = A[i * ndim + j] - 1. / double(ndim) * contract * Cinv[i * ndim + j];
        }
    }
    return DEVA;
}

/*! \brief This function calculates the deviatoric part of a tensor
  @param[in] A Tensor to calculate its deviatoric part
  @param[in] ndim Dimension of the tensor
  @return The deviatoric part of the tensor A
 */
extern vector<double> devSTensor3(const vector<double> &A, int ndim) {
    vector<double> I(ndim * ndim, 0.);
    double trace = 0.;//a(0,0)+a(1,1)+a(2,2);
    for (int i = 0; i < ndim; i++) {
        I[i * ndim + i] = 1; // defining the kronecker delta
        trace += A[i * ndim + i];
    }
    vector<double> devA(ndim * ndim, 0.);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            devA[i * ndim + j] = A[i * ndim + j] - 1. / double(ndim) * trace * I[i * ndim + j];
        }
    }
    return devA;
}

double getTrace(const vector<double> &A, int ndim) {
    double trace = 0;
    for (int i = 0; i < ndim; i++) {
        trace += A[i * ndim + i];
    }
    return trace;
}

/*! \brief Dyadic product between two vectors
  @param[in] t First vector
  @param[in] m Size of the first vector
  @param[in] mm Second vector
  @param[in] l Size of the second vector
  @param[in] C Tensor with the result
 */
extern void dyadicProduct(const vector<double> &t, int m, const vector<double> &mm, int l, vector<double> &C) {
    if (m != l) {
        cout << "ERROR: Wrong dyadic multiplication" << endl;
    }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < l; j++)
            C[i * m + j] = t[i] * mm[j];
}

/*! \brief Transpose of a tensor. It should be a square tensor.
  @param[in] A Tensor to transpose
  @param[in] n Size of the tensor
  @param[in] C Tensor with the result
 */
extern void transpose(const vector<double> &A, int n, vector<double> &C) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = A[j * n + i];
        }
    }
}

/*! \brief Transpose of a tensor. It should be a square tensor.
  @param[in] A Tensor to transpose
  @param[in] n Size of the tensor
  @param[in] C Tensor with the result
 */
extern void transposeStochastic(const vector<double> &A, int n, int nstoch, vector<double> &C) {
    setAll(C,0);
    for (int k = 0; k < nstoch; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                C[k * n * n + i * n + j] = A[k * n * n + j * n + i];
            }
        }
    }
}

/*! \brief Compute GreenLagrangeTensor.
  @param[in] FCurr Deformation Gradient
  @param[in] ndim size of FCurr, must be square
  @param[in&out] E GreenLagrange tensor
 */
extern void computeGL(const vector<double> &FCurr, int ndim, vector<double> &E) {

    vector<double> rightCauchy(ndim *ndim, 0);
    vector<double> delta(ndim *ndim, 0);

    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }
    multSTensor3FirstTranspose(FCurr, ndim, ndim, FCurr, ndim, rightCauchy);
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            E[i * ndim + j] = 0.5 * (rightCauchy[i * ndim + j] - delta[i * ndim + j]);
        }
    }
}

/*! \brief Compute Cauchy Tensor.
  @param[in] FCurr Deformation Gradient
  @param[in] ndim size of FCurr, must be square
  @param[in] PK1 First Piola Kirchhoff tensor
  @param[in&out] Cauchy Cauchy Tensor
 */
extern void computeCauchy(const vector<double> &FCurr, int ndim, const vector<double> &PK1, vector<double> &Cauchy) {
    double detJinv = 1./determinantTensor(FCurr, ndim);
    multSTensor3SecondTranspose(PK1, ndim, ndim, FCurr, ndim, Cauchy);
    multiSTensorScalar(Cauchy,ndim,ndim,detJinv);
}

/*! \brief Compute VMS.
  @param[in] Cauchy Cauchy Tensor
  @param[in] ndim size of Cauchy
  @param[in&out] VMS Von Mises Stress
 */
extern void computeVMS(vector<double> &Cauchy, int ndim, vector<double> &VMS) {
    if (ndim == 3) {
        VMS[0] = sqrt(0.5 * ((Cauchy[0] - Cauchy[4]) * (Cauchy[0] - Cauchy[4]) +
                          (Cauchy[4] - Cauchy[8]) * (Cauchy[4] - Cauchy[8]) +
                          (Cauchy[8] - Cauchy[0]) * (Cauchy[8] - Cauchy[0]) +
                          6 * (Cauchy[1] * Cauchy[1] + Cauchy[2] * Cauchy[2] + Cauchy[5] * Cauchy[5])));
    } else {
        VMS[0] = sqrt(
                Cauchy[0] * Cauchy[0] + Cauchy[3] * Cauchy[3] - Cauchy[0] * Cauchy[3] + 3 * Cauchy[1] * Cauchy[1]);
    }
}

/*! \brief Compute Volumetric Strain.
  @param[in] E GreenLagrange tensor
  @param[in] ndim size of FCurr, must be square
  @param[in&out] volumetricStrain Volumetric Strain
 */
extern void computeVolumetricStrain(vector<double> &E, int ndim, vector<double> &volumetricStrain) {
    setAll(volumetricStrain,0);  /*! Just to make sure*/
    for(int i=0; i<ndim; i++){
        volumetricStrain[0] += E[i+i*ndim];
    }

}

/*! \brief Compute Equivalent Strain
  @param[in] E GreenLagrange tensor
  @param[in] ndim size of FCurr, must be square
  @param[in&out] equivalentStrain Equivalent Strain
 */
extern void computeEquivalentStrain(const vector<double> &E, int ndim, vector<double> &equivalentStrain) {

    if (ndim == 3) {
        equivalentStrain[0] = 0.6666*sqrt(1.5 *(E[0]*E[0]+E[4]*E[4]+E[8]*E[8]) + 3*(E[1]*E[1] + E[2]*E[2] + E[5]*E[5]));
    } else {
        equivalentStrain[0] = 0.6666*sqrt(1.5*(E[0]*E[0]+E[3]*E[3]) + 3*(E[1]*E[1]));
    }

}

/*! \brief Compute Stochastic GreenLagrangeTensor.
  @param[in] FCurr Stochastic Deformation Gradient
  @param[in] ndim Physical dimension
  @param[in] nstoch Stochastic dimension
  @param[in] Cstoch Stochastic Tensor for Algebra computation
  @param[in&out] E Stochastic GreenLagrange tensor
 */
extern void computeGLstoch(const vector<double> &FCurr, int ndim, int nstoch, const vector<double> &C_Stoch, vector<double> &E) {
    vector<double> delta(ndim *ndim, 0);
    for (int i=0; i< ndim; i++)
    {
        delta[i * ndim + i]  = 1;
    }

    vector<double>rightCauchy(ndim*ndim*nstoch,0);

    multSTensor3FirstTransposeStochastic(FCurr, ndim, ndim, FCurr, ndim, rightCauchy, nstoch, C_Stoch);


    for (int m = 0; m < nstoch; m++) {
        for (int i = 0; i < ndim; i++) {
            for (int j = 0; j < ndim; j++) {
                if (m == 0) {
                    E[m * ndim * ndim + i * ndim + j] =
                            0.5 * (rightCauchy[m * ndim * ndim + i * ndim + j] - delta[i * ndim + j]);
                } else {
                    E[m * ndim * ndim + i * ndim + j] = 0.5 * (rightCauchy[m * ndim * ndim + i * ndim + j]);
                }
            }
        }
    }

}

/*! \brief Compute Stochastic Cauchy Tensor.
 @param[in] FCurr Stochastic Deformation Gradient
  @param[in] ndim Physical dimension
  @param[in] nstoch Stochastic dimension
  @param[in] Cstoch Stochastic Tensor for Algebra computation
  @param[in] PK1 Stochastic PK1 Tensor
  @param[in&out] Cauchy Stochastic Cauchy tensor
 */
extern void computeCauchyStoch(const vector<double> &FCurr, int ndim, int nstoch, const vector<double> &C_stoch, const vector<double> &PK1, vector<double> &Cauchy) {
    vector<double> detJ(nstoch,0);
    vector<double> detJinv(nstoch,0);
    vector<double> ones(nstoch,0);
    vector<double> J_Cauchy(nstoch*ndim*ndim,0);
    ones[0]=1;

    DeterminantStochasticTensor(FCurr, ndim, nstoch, detJ, C_stoch);

    DivisionRandomVariableRandomVariableStochastic(ones, detJ, C_stoch, detJinv);
    multSTensor3SecondTransposeStochastic(PK1, ndim, ndim, FCurr, ndim, J_Cauchy, nstoch, C_stoch);
    multiSTensorRandomVariableStochastic(J_Cauchy, ndim, ndim, detJinv, nstoch, C_stoch,
                                     Cauchy);
}

/*! \brief Compute Stochastic VMS.
  @param[in] Cauchy Stochastic Cauchy tensor
  @param[in] ndim Physical dimension
  @param[in] nstoch Stochastic dimension
  @param[in] Cstoch Stochastic Tensor for Algebra computation
  @param[in&out] VMS Stochastic VMS
 */
extern void computeStochasticVMS(vector<double> &Cauchy, int ndim, int nstoch, vector<double> &VMS, const vector<double> &C_stoch) {
    vector<double> sig_xx(nstoch, 0);
    vector<double> sig_yy(nstoch, 0);
    vector<double> sig_zz(nstoch, 0);
    vector<double> sig_xy(nstoch, 0);
    vector<double> sig_yz(nstoch, 0);
    vector<double> sig_xz(nstoch, 0);
    vector<double> sig_xx_square(nstoch, 0);
    vector<double> sig_yy_square(nstoch, 0);
    vector<double> sig_zz_square(nstoch, 0);
    vector<double> sig_xy_square(nstoch, 0);
    vector<double> sig_yz_square(nstoch, 0);
    vector<double> sig_xz_square(nstoch, 0);
    vector<double> A(nstoch, 0);
    vector<double> B(nstoch, 0);
    vector<double> C(nstoch, 0);
    vector<double> A_square(nstoch, 0);
    vector<double> B_square(nstoch, 0);
    vector<double> C_square(nstoch, 0);
    vector<double> VMS_Square(nstoch, 0);
    if (ndim == 3) {
        for (int i = 0; i < nstoch; i++) {
            sig_xx[i] = Cauchy[i * ndim * ndim + 0];
            sig_yy[i] = Cauchy[i * ndim * ndim + 4];
            sig_zz[i] = Cauchy[i * ndim * ndim + 8];
            sig_xy[i] = Cauchy[i * ndim * ndim + 1];
            sig_yz[i] = Cauchy[i * ndim * ndim + 2];
            sig_xz[i] = Cauchy[i * ndim * ndim + 5];
            A[i] = sig_xx[i] - sig_yy[i];
            B[i] = sig_yy[i] - sig_zz[i];
            C[i] = sig_zz[i] - sig_xx[i];
        }

        multiRandomVariableRandomVariableStochastic(A, A, C_stoch, A_square);
        multiRandomVariableRandomVariableStochastic(B, B, C_stoch, B_square);
        multiRandomVariableRandomVariableStochastic(C, C, C_stoch, C_square);

        multiRandomVariableRandomVariableStochastic(sig_xy, sig_xy, C_stoch, sig_xy_square);
        multiRandomVariableRandomVariableStochastic(sig_yz, sig_yz, C_stoch, sig_yz_square);
        multiRandomVariableRandomVariableStochastic(sig_xz, sig_xz, C_stoch, sig_xz_square);

        for (int i = 0; i < nstoch; i++) {
            VMS_Square[i] = 0.5 * (A[i] + B[i] + C[i] + 6 * (sig_xy[i] + sig_xz[i] + sig_yz[i]));
        }

        squareRootIntegral(VMS_Square, VMS, C_stoch);


    } else {
        for (int i = 0; i < nstoch; i++) {
            sig_xx[i] = Cauchy[i * ndim * ndim + 0];
            sig_yy[i] = Cauchy[i * ndim * ndim + 3];
            sig_xy[i] = Cauchy[i * ndim * ndim + 1];
            A[i] = sig_xx[i] - sig_yy[i];
            B[i] = sig_yy[i] - sig_zz[i];
            C[i] = sig_zz[i] - sig_xx[i];
        }
        multiRandomVariableRandomVariableStochastic(sig_xx, sig_xx, C_stoch, sig_xx_square);
        multiRandomVariableRandomVariableStochastic(sig_yy, sig_yy, C_stoch, sig_yy_square);
        multiRandomVariableRandomVariableStochastic(sig_xx, sig_yy, C_stoch, A);
        multiRandomVariableRandomVariableStochastic(sig_xy, sig_xy, C_stoch, sig_xy_square);


        for (int i = 0; i < nstoch; i++) {
            VMS_Square[i] = sig_xx_square[i] + sig_yy_square[i] + 3 * sig_xy_square[i] - A[i];
        }

        squareRootIntegral(VMS_Square, VMS, C_stoch);
    }
}

/*! \brief Compute Stochastic Volumetric Strain.
  @param[in] E Stochastic GreenLagrangeTensor
  @param[in] ndim Physical dimension
  @param[in] nstoch Stochastic dimension
  @param[in] Cstoch Stochastic Tensor for Algebra computation
  @param[in&out] volumetricStrain Stochastic VolumetricStrain
 */
    extern void computeVolumetricStrainStochastic(const vector<double> &E, int ndim, int nstoch, vector<double> &volumetricStrain) {
        setAll(volumetricStrain,0);  /*! Just to make sure*/
        for(int i=0; i<nstoch; i++) {
            for (int j = 0; j < ndim; j++) {
                volumetricStrain[i] += E[i*ndim*ndim + j + j * ndim];
            }
        }
    }

/*! \brief Compute Stochastic Equivalent Strain.
  @param[in] E Stochastic GreenLagrangeTensor
  @param[in] ndim Physical dimension
  @param[in] nstoch Stochastic dimension
  @param[in&out] EquivalentStrain Stochastic Equivalent Strain
 */
    extern void computeEquivalentStrainStochastic(const vector<double> &E, int ndim, int nstoch, const vector<double> &C_stoch, vector<double> &equivalentStrain) {
        setAll(equivalentStrain,0);
        vector<double> eps_xx(nstoch, 0);
        vector<double> eps_yy(nstoch, 0);
        vector<double> eps_zz(nstoch, 0);
        vector<double> eps_xy(nstoch, 0);
        vector<double> eps_yz(nstoch, 0);
        vector<double> eps_xz(nstoch, 0);
        vector<double> eps_xx_square(nstoch, 0);
        vector<double> eps_yy_square(nstoch, 0);
        vector<double> eps_zz_square(nstoch, 0);
        vector<double> eps_xy_square(nstoch, 0);
        vector<double> eps_yz_square(nstoch, 0);
        vector<double> eps_xz_square(nstoch, 0);
        vector<double> A(nstoch, 0);
        vector<double> B(nstoch, 0);
        vector<double> C(nstoch, 0);
        vector<double> A_square(nstoch, 0);
        vector<double> B_square(nstoch, 0);
        vector<double> C_square(nstoch, 0);
        vector<double> eps_Eq_Square(nstoch, 0);
        if (ndim == 3) {
            for (int i = 0; i < nstoch; i++) {
                eps_xx[i] = E[i * ndim * ndim + 0];
                eps_yy[i] = E[i * ndim * ndim + 4];
                eps_zz[i] = E[i * ndim * ndim + 8];
                eps_xy[i] = E[i * ndim * ndim + 1];
                eps_yz[i] = E[i * ndim * ndim + 2];
                eps_xz[i] = E[i * ndim * ndim + 5];
                A[i] = eps_xx[i] - eps_yy[i];
                B[i] = eps_yy[i] - eps_zz[i];
                C[i] = eps_zz[i] - eps_xx[i];
            }

            multiRandomVariableRandomVariableStochastic(A, A, C_stoch, A_square);
            multiRandomVariableRandomVariableStochastic(B, B, C_stoch, B_square);
            multiRandomVariableRandomVariableStochastic(C, C, C_stoch, C_square);

            multiRandomVariableRandomVariableStochastic(eps_xy, eps_xy, C_stoch, eps_xy_square);
            multiRandomVariableRandomVariableStochastic(eps_yz, eps_yz, C_stoch, eps_yz_square);
            multiRandomVariableRandomVariableStochastic(eps_xz, eps_xz, C_stoch, eps_xz_square);

            for (int i = 0; i < nstoch; i++) {
                eps_Eq_Square[i] = 0.5 * (A[i] + B[i] + C[i] + 6 * (eps_xy[i] + eps_xz[i] + eps_yz[i]));
            }

            squareRootIntegral(eps_Eq_Square, equivalentStrain, C_stoch);


        } else {
            for (int i = 0; i < nstoch; i++) {
                eps_xx[i] = E[i * ndim * ndim + 0];
                eps_yy[i] = E[i * ndim * ndim + 3];
                eps_xy[i] = E[i * ndim * ndim + 1];
                A[i] = eps_xx[i] - eps_yy[i];
                B[i] = eps_yy[i] - eps_zz[i];
                C[i] = eps_zz[i] - eps_xx[i];
            }
            multiRandomVariableRandomVariableStochastic(eps_xx, eps_xx, C_stoch, eps_xx_square);
            multiRandomVariableRandomVariableStochastic(eps_yy, eps_yy, C_stoch, eps_yy_square);
            multiRandomVariableRandomVariableStochastic(eps_xx, eps_yy, C_stoch, A);
            multiRandomVariableRandomVariableStochastic(eps_xy, eps_xy, C_stoch, eps_xy_square);


            for (int i = 0; i < nstoch; i++) {
                eps_Eq_Square[i] = eps_xx_square[i] + eps_yy_square[i] + 3 * eps_xy_square[i] - A[i];
            }

            squareRootIntegral(eps_Eq_Square, equivalentStrain, C_stoch);
        }
    }
/*! \brief WeightedAverage of a vector
  @param[in] a Vector
  @param[in] size Size of the vector
  @return The weighted average
 */
extern double weightedAverage(const vector<double> &a, int size) {
    double temp = 0;

    for (int i = 0; i < size; i++) {
        temp += a[i];
    }
    temp = temp / size;

    return temp;
}

/*! \brief Euclidean norm of a vector
  @param[in] a Vector
  @param[in] size Size of the vector
  @return norm
 */
extern double norm2(const vector<double> &a, int size) {
    double temp = 0;
    for (int i = 0; i < size; i++) {
        temp += a[i] * a[i];
    }
    temp = sqrt(temp);
    return temp;
}

/*! \brief Determinant of a tensor. ndim x ndim
  @param[in] A Tensor
  @param[in] size Size of the tensor size x size
  @return determinant
 */
extern double determinantTensor(const vector<double> &A, int size) {
    double det = 0;

    if (size == 1) {
        det = A[0];
    } else if (size == 2) {
        det = A[0] * A[3] - A[2] * A[1];
    } else {
        det = A[0] * (A[4] * A[8] - A[7] * A[5]) -
              A[1] * (A[3] * A[8] - A[6] * A[5]) +
              A[2] * (A[3] * A[7] - A[6] * A[4]);
    }
    return det;
}

/*! \brief Multiplication of a fourth order tensor by a second order tensor (2D and 3D) ndim x ndim x ndim xndim  XXX ndim x ndim
  @param[in] AA Fourth order tensor
  @param[in] B Second order tensor
  @param[out] C  Result
  @param[in] ndim Dimension of the domain
 */
extern void multTensor4Tensor(const classTensor4 *AA, const vector<double> &B, vector<double> &C, int ndim) {
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < ndim; k++) {
                for (int l = 0; l < ndim; l++) {
                    C[i * ndim + j] += AA->get(i, j, k, l) * B[k * ndim + l];
                }
            }
        }
    }
}

/*! \brief Sum two tensors (2D and 3D) ndim x ndim
  @param[in] A Second order tensor
  @param[in] B Second order tensor
  @param[out] C Result
  @param[in] ndim Dimension of the domain
  @param[in] sum If it is equal to 1 or -1
 */
extern void sumTensorTensor(const vector<double> &A, const vector<double> &B, vector<double> &C, int ndim, int sum) {
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            C[i * ndim + j] = A[i * ndim + j] + sum * B[i * ndim + j];
        }
    }
}

/*! \brief Multiplication of a matrix by a vector. ndim x ndim XXX ndim
  @param[in] A Second order tensor
  @param[in] m Number of rows
  @param[in] n Number of columns
  @param[in] mm Vector
  @param[in] l Size of the vector
  @param[out] val Result
 */
extern void multTensorVector(const vector<double> &A, int m, int n, const vector<double> &mm, int l, vector<double> &val) {
    if (n != l) {
        cout << "ERROR: Wrong matrix-vector multiplication" << endl;
        exit(-1);
    }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            val[i] += A[i * m + j] * mm[j];
}

/*! \brief Multiplication of a matrix by a vector. ndim x ndim XXX ndim
  @param[in] A Second order tensor
  @param[in] m Number of rows
  @param[in] n Number of columns
  @param[in] mm Vector
  @param[in] l Size of the vector
  @param[out] val Result
 */
extern void multTensorVectorStochastic(const vector<double> &A, int ndim, int nstoch, const vector<double> &B, const vector<double> &C, vector<double> &D) {

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            for (int k = 0; k < nstoch; k++) {
                for (int l = 0; l < nstoch; l++) {
                    for (int m = 0; m < nstoch; m++) {
                        D[i + ndim * m] += A[k * ndim * ndim + i * ndim + j] * B[l * ndim + j] * C[k + l *nstoch + m *nstoch*nstoch];
                    }
                }
            }
        }
    }
}

/*! \brief transpose a general matrix. m x n -> n x m
  @param[in] A First Matrix
  @param[in] m Number of rows first matrix
  @param[in] n Number of columns first matrix
  @param[out] C transposed matrix
 */
extern void transposeGeneralMatrix(const vector<double> &A, int m, int n, vector<double> &C) {

    fill(C.begin(), C.end(), 0);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[j * m + i] = A[i * n + j];
        }
    }
}

/*! \brief Multiplication between two matrices. m x n XXX kk x ll
  @param[in] A First Matrix
  @param[in] m Number of rows first matrix
  @param[in] n Number of columns first matrix
  @param[in] B Second Matrix
  @param[in] kk Number of rows second matrix
  @param[in] ll Number of columns second matrix
  @param[out] C Result
 */
extern void multGeneralMatrices(const vector<double> &A, int m, int n, const vector<double> &B, int kk, int ll, vector<double> &C) {
    if (n != kk) {
        cout << "ERROR: Wrong general matrix-matrix multiplication" << endl;
        exit(-1);
    }
    fill(C.begin(), C.end(), 0);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < ll; j++) {
            for (int k = 0; k < kk; k++) {
                C[i * ll + j] += A[i * kk + k] * B[k * ll + j];
            }
        }
    }
}

/*! \brief Multiplication between two tensors. Just used for square tensors, e.g., ndim x ndim XXX ndim x ndim
  @param[in] A First tensor
  @param[in] m Number of rows first tensor
  @param[in] n Number of columns first tensor
  @param[in] B Second tensor
  @param[in] l Number of rows second tensor
  @param[out] C Result
 */
extern void multTensorTensor3(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C) {
    if (n != l) {
        cout << "ERROR: Wrong matrix multiplication" << endl;
        exit(-1);
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < m; k++) {
                C[i * m + j] += A[i * m + k] * B[k * m + j];
            }
        }
    }
}

/*! \brief Multiplication between two 'stochastic' tensors (third order tensor)
  @param[in] A First tensor
  @param[in] m Number of rows first tensor
  @param[in] n Number of columns first tensor
  @param[in] B Second tensor
  @param[in] l Number of rows second tensor
  @param[out] C Result
  @param[in] nstoch number of stochastic components
  @param[in] multiplier Third order tensor that computes the product between two random variables
 */
extern void
multTensorTensor3Stochastic(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C, int nstoch,
                            const vector<double>& multiplier) {
    if (n != l) {
        cout << "ERROR: Wrong matrix multiplication" << endl;
        exit(-1);
    }
    setAll(C,0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int p = 0; p < nstoch; p++) {
                for (int u = 0; u < nstoch; u++) {
                    for (int v = 0; v < nstoch; v++) {
                        for (int k = 0; k < m; k++) {
                            C[p * m * m + i * m + j] +=
                                    multiplier[p * nstoch * nstoch + v * nstoch + u] * A[u * m * m + i * m + k] *
                                    B[v * m * m + k * m + j];
                        }
                    }
                }
            }
        }
    }
}

/*! \brief Multiplication between the transpose of the first tensor and the second one. Just used for square tensors, e.g., ndim x ndim XXX ndim x ndim
  @param[in] A First tensor
  @param[in] m Number of rows first tensor
  @param[in] n Number of columns first tensor
  @param[in] B Second tensor
  @param[in] l Number of rows second tensor
  @param[out] C Result
 */
extern void multSTensor3FirstTranspose(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C) {
    if (n != l) {
        cout << "ERROR: Wrong matrix multiplication" << endl;
        exit(-1);
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < l; k++) {
                C[i * m + j] += A[k * l + i] * B[k * m + j];
            }
        }
    }
}

/*! \brief Multiplication between the transpose of the first tensor and the second one for stochastic quantities. Just used for square tensors, e.g., ndim x ndim XXX ndim x ndim
  @param[in] A First tensor
  @param[in] m Number of rows first tensor
  @param[in] n Number of columns first tensor
  @param[in] B Second tensor
  @param[in] l Number of rows second tensor
  @param[out] C Result
  @param[in] nstoch number of stochastic components
  @param[in] multiplier Third order tensor that computes the product between two random variables
 */
extern void multSTensor3FirstTransposeStochastic(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C,
                                     int nstoch, const vector<double>& multiplier) {
    if (n != l) {
        cout << "ERROR: Wrong matrix multiplication" << endl;
        exit(-1);
    }
    setAll(C,0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int p = 0; p < nstoch; p++) {
                for (int k = 0; k < l; k++) {
                    for (int u = 0; u < nstoch; u++) {
                        for (int v = 0; v < nstoch; v++) {
                            C[p * m * m + i * m + j] +=
                                    multiplier[p * nstoch * nstoch + v * nstoch + u] * A[u * m * m + k * l + i] *
                                    B[v * m * m + k * m + j];
                        }
                    }
                }
            }
        }
    }
}

/*! \brief Multiplication between first tensor and the second transposed of the second one for stochastic quantities. Just used for square tensors, e.g., ndim x ndim XXX ndim x ndim
  @param[in] A First tensor
  @param[in] m Number of rows first tensor
  @param[in] n Number of columns first tensor
  @param[in] B Second tensor
  @param[in] l Number of rows second tensor
  @param[out] C Result
  @param[in] nstoch number of stochastic components
  @param[in] multiplier Third order tensor that computes the product between two random variables
 */
extern void multSTensor3SecondTransposeStochastic(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C,
                                     int nstoch, const vector<double>& multiplier) {
    if (n != l) {
        cout << "ERROR: Wrong matrix multiplication" << endl;
        exit(-1);
    }
    setAll(C,0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int p = 0; p < nstoch; p++) {
                for (int k = 0; k < l; k++) {
                    for (int u = 0; u < nstoch; u++) {
                        for (int v = 0; v < nstoch; v++) {
                            C[p * m * m + i * m + j] +=
                                    multiplier[p * nstoch * nstoch + v * nstoch + u] * A[u * m * m + i * l + k] *
                                    B[v * m * m + j * m + k];
                        }
                    }
                }
            }
        }
    }
}

/*! \brief Exponential of random variable
  @param[in] A Input probabilist ditribution
  @param[inout] Result
  @param[in] C Third order tensor that computes the product between two random variables
 */
extern void expIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C){
    int size = A.size();
    vector<double> initial_point(size,0);
    initial_point[0] = A[0];
    int s = 1;
    double ds = 0.0001;
    int num_step = s/ds ;
    vector<double> u_old(size,0);
    vector<double> u_new(size,0);
    u_old[0] = exp(A[0]);
    u_new[0] = exp(A[0]);
    double sum = 0;

    for (int l=0; l<num_step; l++) {
        for (int k = 0; k < size; k++) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    sum = sum + ds* (A[j]-initial_point[j]) * C[i + size * j + size*size*k]*u_new[i];
                }
            }
            u_new[k] = u_new[k] + sum;
            sum = 0;
        }
    }
    Result = u_new;

};

/*! \brief Sigmoid of random variable (we have x and "return" 1/(1+exp(x))
  @param[in] A Input probabilist ditribution
  @param[inout] Result
  @param[in] C Third order tensor that computes the product between two random variables
 */
extern void sigmoid(const vector<double> &A, vector<double> &Result, const vector<double>& C){
    int size = A.size();
    vector<double> ones(size,0);
    ones[0] = 1;
    vector<double> exp_A(size,0);
    expIntegral(A, exp_A, C);
    linearCombinationVariableStochastic(exp_A, 1, 1);
    DivisionRandomVariableRandomVariableStochastic(ones, exp_A, C, Result);
};

/*! \brief Square root of a probabilistic distribution (careful for distributions having large non-zero probabilities of giving negative results)
  @param[in] A Input probabilist ditribution
  @param[inout] Result
  @param[in] C Third order tensor that computes the product between two random variables
 */
extern void squareRootIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C){
    int size = A.size();
    vector<double> initial_point(size,0);
    double tol = 1e-2;
    initial_point[0] = A[0];
    if(abs(A[0]) < tol){
        initial_point[0] = initial_point[0] + 0.1;
    }
    vector<double> ones(size,0);
    ones[0] = 1;
    vector<double> inversion(size,0);
    int s = 1;
    double ds = 0.01;
    int num_step = s/ds ;
    vector<double> u_old(size,0);
    vector<double> u_new(size,0);
    if(abs(A[0]) < tol){
        u_old[0] = sqrt(A[0]+0.1);
        u_new[0] = sqrt(A[0]+0.1);
        }else{
            u_old[0] = sqrt(A[0]);
            u_new[0] = sqrt(A[0]);
    }

    double sum = 0;

    for (int l=0; l<num_step; l++) {
        DivisionRandomVariableRandomVariableStochastic(ones, u_new, C,
                                                       inversion);
        for (int k = 0; k < size; k++) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                sum = sum + ds* (A[j]-initial_point[j]) * C[i + size * j + size*size*k]*0.5*inversion[i];
                }
            }
            u_new[k] = u_new[k] + sum;
            sum = 0;
        }
    }
Result = u_new;

};

/*! \brief Ln of a probabilistic distribution (careful for distributions having large non-zero probabilities of giving negative results)
  @param[in] A Input probabilist ditribution
  @param[inout] Result
  @param[in] C Third order tensor that computes the product between two random variables
 */
extern void logIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C){
    int size = A.size();
    vector<double> initial_point(size,0);
    vector<double> u_new(size,0);
    initial_point[0] = 1;
    vector<double> ones(size,0);
    ones[0] = 1;
    vector<double> inversion(size,0);
    vector<double> A_curr(size,0);
    int s = 1;
    double ds = 0.01;
    int num_step = s/ds ;
    double s_curr = 0;
    double sum = 0;

    for (int l=0; l<num_step; l++) {
    s_curr = ds *l;
        for (int a=0; a < size; a++){
            A_curr[a] = initial_point[a]*(1-s_curr) + A[a]*s_curr;
        }
        DivisionRandomVariableRandomVariableStochastic(ones, A_curr, C,
                                                       inversion);
        for (int k = 0; k < size; k++) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                sum = sum + ds* (A[j]-initial_point[j]) * C[i + size * j + size*size*k]*inversion[i];
                }
            }
            u_new[k] = u_new[k] + sum;
            sum = 0;
        }
    }
Result = u_new;

};

/*! \brief Square root of a probabilistic distribution (careful for distributions having large non-zero probabilities of giving negative results)
  @param[in] A Input probabilist ditribution
  @param[inout] Result
  @param[in] C Third order tensor that computes the product between two random variables
 */
extern void InterpolationSqrt(vector<double> &A, vector<double> &Result, vector<double> C){

    int size = A.size();
    vector<double> A_square(size,0);
    vector<double> A_cube(size,0);
    vector<double> A_fourth(size,0);
    vector<double> A_fifth(size,0);

    multiRandomVariableRandomVariableStochastic(A, A, C, A_square);
    multiRandomVariableRandomVariableStochastic(A_square, A, C, A_cube);
    multiRandomVariableRandomVariableStochastic(A_cube, A, C, A_fourth);
    multiRandomVariableRandomVariableStochastic(A_fourth, A, C, A_fifth);


    for(int i=0; i<A.size();i++){
        Result[i] =  162.648*A[i] -114.9309*A_square[i] - 92.0879*A_cube[i] + 141.467*A_fourth[i] - 44.9639*A_fifth[i];
    }
    Result[0] = Result[0]-52;

};


/*! \brief Cubic root of a probabilistic distribution (careful for same reasons as square root)
  @param[in] A Input probabilist ditribution
  @param[inout] Result
  @param[in] C Third order tensor that computes the product between two random variables
 */
extern void cubicRootIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C){
    int size = A.size();
    vector<double> initial_point(size,0);
    initial_point[0] = A[0];
    vector<double> ones(size,0);
    ones[0] = 1;
    vector<double> inversion(size,0);
    int s = 1;
    double ds = 0.001;
    int num_step = s/ds ;
    vector<double> u_old(size,0);
    vector<double> u_new(size,0);
    u_old[0] = cbrt(A[0]);
    u_new[0] = cbrt(A[0]);
    double sum = 0;
    double s_curr = 0;
    vector<double> Acurr(size,0);
    vector<double> g(size,0);
    for (int l=0; l<num_step; l++) {
        s_curr = ds *l;
        for (int a=0; a < size; a++){
            Acurr[a] = initial_point[a]*(1-s_curr) + A[a]*s_curr;
        }
        DivisionRandomVariableRandomVariableStochastic(ones, Acurr, C,
                                                       inversion);
        multiRandomVariableRandomVariableStochastic(inversion, u_new, C, g);
        for (int k = 0; k < size; k++) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    sum = sum + ds* (A[j]-initial_point[j]) * C[i + size * j + size*size*k]*0.33*g[i];
                }
            }
            u_new[k] = u_new[k] + sum;
            sum = 0;
        }
        for (int a=0; a < size; a++){
            g[a] = 0;
        }
    }
    Result = u_new;

};

/*! \brief Cnth power of a probabilistic distribution (careful for same reasons as square root)
  @param[in] A Input probabilist ditribution
  @param[inout] Result
  @param[in] C Third order tensor that computes the product between two random variables
 */
extern void nthPowerIntegral(const vector<double> &A, vector<double> &Result, const vector<double>& C, const double &n){
    int size = A.size();
    vector<double> initial_point(size,0);
    initial_point[0] = 1.0;
    vector<double> ones(size,0);
    ones[0] = 1.0;
    vector<double> inversion(size,0);
    int s = 1.0;
    double ds = 0.05;
    int num_step = s/ds ;
    vector<double> u_old(size,0);
    vector<double> u_new(size,0);
    u_old[0] = 1;
    u_new[0] = 1;
    double sum = 0;
    double s_curr = 0;
    vector<double> Acurr(size,0);
    vector<double> g(size,0);
    for (int l=0; l<num_step; l++) {
        s_curr = ds *(l+1);
        for (int a=0; a < size; a++){
            Acurr[a] = initial_point[a]*(1-s_curr) + A[a]*s_curr;
        }
        DivisionRandomVariableRandomVariableStochastic(ones, Acurr, C,
                                                       inversion);
        multiRandomVariableRandomVariableStochastic(inversion, u_new, C, g);
        for (int k = 0; k < size; k++) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    sum = sum + ds* (A[j]-initial_point[j]) * C[i + size * j + size*size*k]*n*g[i];
                }
            }
            u_new[k] = u_new[k] + sum;
            sum = 0;
        }
        for (int a=0; a < size; a++){
            g[a] = 0;
        }
    }
    Result = u_new;

};

/*! \brief Compute sqrt(1-x) for probabilistic distributions
  @param[in] A Probabilistic distribution x
  @param[in] C Third order tensor that computes product
  @param[in] Result The result of the operation
 */
extern void DL_sqrt_1minusx(vector<double> &A, vector<double> C, vector<double> &Result){
    int size = A.size();
    vector<double> A_square(size,0);
    vector<double> A_cube(size,0);
    vector<double> A_fourth(size,0);
    vector<double> A_fifth(size,0);
    vector<double> A_sixth(size,0);
    vector<double> A_seventh(size,0);
    vector<double> A_eighth(size,0);
    vector<double> A_ninth(size,0);
    vector<double> A_tenth(size,0);

    multiRandomVariableRandomVariableStochastic(A, A, C, A_square);
    multiRandomVariableRandomVariableStochastic(A_square, A, C, A_cube);
    multiRandomVariableRandomVariableStochastic(A_cube, A, C, A_fourth);
    multiRandomVariableRandomVariableStochastic(A_fourth, A, C, A_fifth);
    multiRandomVariableRandomVariableStochastic(A_fifth, A, C, A_sixth);
    multiRandomVariableRandomVariableStochastic(A_sixth, A, C, A_seventh);
    multiRandomVariableRandomVariableStochastic(A_seventh, A, C, A_eighth);
    multiRandomVariableRandomVariableStochastic(A_eighth, A, C, A_ninth);
    multiRandomVariableRandomVariableStochastic(A_ninth, A, C, A_tenth);

    for(int i=0; i<A.size();i++){
        Result[i] =  -0.5*A[i] - (1/8)*A_square[i] - 1/16*A_cube[i] - 0*5/128*A_fourth[i] - 0*(21/1024)*A_fifth[i];
    }
    Result[0] = Result[0]+1;
};

/*! \brief Compute sqrt(x) for x between 0 and 0.5
  @param[in] A Probabilistic distribution x
  @param[in] C Third order tensor that computes product
  @param[in] Result The result of the operation
 */
extern void polyfitsqrt(vector<double> &A, vector<double> C, vector<double> &Result){
    int size = A.size();
    vector<double> A_square(size,0);
    vector<double> A_cube(size,0);
    vector<double> A_fourth(size,0);
    vector<double> A_fifth(size,0);
    vector<double> A_sixth(size,0);
    //vector<double> A_seventh(size,0);
    //vector<double> A_eighth(size,0);
    //vector<double> A_ninth(size,0);
    //vector<double> A_tenth(size,0);

    multiRandomVariableRandomVariableStochastic(A, A, C, A_square);
    multiRandomVariableRandomVariableStochastic(A_square, A, C, A_cube);
    multiRandomVariableRandomVariableStochastic(A_cube, A, C, A_fourth);
    multiRandomVariableRandomVariableStochastic(A_fourth, A, C, A_fifth);
    multiRandomVariableRandomVariableStochastic(A_fifth, A, C, A_sixth);
    //multiRandomVariableRandomVariableStochastic(A_sixth, A, C, A_seventh);
    //multiRandomVariableRandomVariableStochastic(A_seventh, A, C, A_eighth);
    //multiRandomVariableRandomVariableStochastic(A_eighth, A, C, A_ninth);
    //multiRandomVariableRandomVariableStochastic(A_ninth, A, C, A_tenth);

    for(int i=0; i<A.size();i++){
        Result[i] =  4.73*A[i] - 21.79*A_square[i] +63.97*A_cube[i] - 100.88*A_fourth[i] + 79.59*A_fifth[i] - 24.63*A_sixth[i];
    }
    //Result[0] = Result[0]+0.0581;
};

/*! \brief Compute sqrt(x) for x between 0 and 0.5
  @param[in] A Probabilistic distribution x
  @param[in] C Third order tensor that computes product
  @param[in] Result The result of the operation
 */
extern void matlabfit(double r, vector<double> mu, vector<double> C, vector<double> &Result){
    int size = mu.size();
    vector<double> mu_square(size,0);
    vector<double> mu_cube(size,0);
    vector<double> mu_fourth(size,0);
    vector<double> mu_fifth(size,0);
    vector<double> mu_sixth(size,0);


    multiRandomVariableRandomVariableStochastic(mu, mu, C, mu_square);
    multiRandomVariableRandomVariableStochastic(mu_square, mu, C, mu_cube);
    multiRandomVariableRandomVariableStochastic(mu_cube, mu, C, mu_fourth);
    multiRandomVariableRandomVariableStochastic(mu_fourth, mu, C, mu_fifth);

    double r_square = pow(r,2);
    double r_cube = pow(r,3);
    double r_fourth = pow(r,4);
    double r_fifth = pow(r,5);

    for(int i=0; i<mu.size();i++) {
        //Result[i] = -89.77 * mu[i] + 11.79 * r * mu[i] + 301 * mu_square[i] + 3.676 * r_square * mu[i] -
                   // 46.75 * r * mu_square[i] - 441.2 * mu_cube[i] + 0.1383 * r_cube * mu[i] -
                    //4.532 * r_square * mu_square[i] + 51.16 * r * mu_cube[i] + 233.3 * mu_fourth[i];
         Result[i] =  -308.4*mu[i] + 53.14 * r * mu[i] + 1455*mu_square[i] + 14.61*r_square * mu[i] - 274.3*r *mu_square[i] - 3400*mu_cube[i] + 2.35 *r_cube *mu[i] - 42.88*r_square*mu_square[i] +577.3 * r *mu_cube[i] + 3906 * mu_fourth[i] - 1.245*r_fourth*mu[i] + 7.515*r_cube *mu_square[i] +4.315 * r_square*mu_cube[i] - 385.1*r*mu_fourth[i] - 1769 * mu_fifth[i];
    }

    Result[0] = Result[0] + 26.49 - 3.847*r - 0.5*r_square - 1.262 * r_cube + 0.299*r_fourth + 0.00638 * r_fifth ;

    //Result[0] = Result[0] + 10.24 - 0.4629*r - 1.104*r_square +0.3184 * r_cube - 0.089*r_fourth;
    //Result[0] = Result[0]+0.0581;

};
/*! \brief Compute sqrt(x) for x between 0 and 0.5
  @param[in] A Probabilistic distribution x
  @param[in] C Third order tensor that computes product
  @param[in] Result The result of the operation
 */
extern void piecewisefit(double r, vector<double> c, vector<double> C, vector<double> &Result){
    int size = c.size();
    vector<double> c_square(size,0);
    vector<double> c_cube(size,0);
    vector<double> c_fourth(size,0);
    vector<double> c_fifth(size,0);
    vector<double> c_sixth(size,0);

    vector<double> Approx_left(size,0);
    vector<double> Approx_right(size,0);
    double ponderation = 0;
    double opposite = 0;
    ponderation = 1./(1+exp(-10*(r-2)));
    opposite = 1-ponderation;

    multiRandomVariableRandomVariableStochastic(c, c, C, c_square);
    multiRandomVariableRandomVariableStochastic(c_square, c, C, c_cube);
    multiRandomVariableRandomVariableStochastic(c_cube, c, C, c_fourth);
    multiRandomVariableRandomVariableStochastic(c_fourth, c, C, c_fifth);

    double r_square = pow(r,2);
    double r_cube = pow(r,3);
    double r_fourth = pow(r,4);
    double r_fifth = pow(r,5);

    for(int i=0; i<c.size();i++) {
        Approx_left[i] =  8.47239335*c[i]  -33.58843728 * r * c[i] -3.42862098*c_square[i]  -2.78289426*r_square * c[i] + 28.74281167*r *c_square[i]  -1.66435262*c_cube[i] + 2.95018441 *r_cube *c[i] -0.74586295*r_square*c_square[i] -10.29536899 * r *c_cube[i] + 1.4132674 * c_fourth[i] +0.6575345*r_fourth*c[i] -1.42889958*r_cube *c_square[i] +0.77001881* r_square*c_cube[i] +1.2698798*r*c_fourth[i] -0.25225718 * c_fifth[i];
        Approx_right[i] =  48.55190949 * c[i]  + 12.2104203*r * c[i]  -63.71634681*c_square[i] -2.56950708*r_square * c[i] -5.3740038*r *c_square[i]  + 37.82952312*c_cube[i] + 7.39139734*r_cube *c[i] -14.52935508*r_square*c_square[i] +15.32140001*r *c_cube[i] -15.64975996*c_fourth[i] -2.78251198*r_fourth*c[i] +6.08457483*r_cube *c_square[i] -6.15529065*r_square*c_cube[i] + 2.40404198*r*c_fourth[i] + 0.97871166*c_fifth[i];
    }

    Approx_left[0] = Approx_left[0] -4.11752294 + 14.09682494*r + 3.14625107*r_square -1.30161716*r_cube - -0.65432379*r_fourth -0.12694093*r_fifth;
    Approx_right[0] = Approx_right[0] + 388.80279454 -845.64060137*r + 699.80960989*r_square -290.09981251*r_cube + 58.78778992*r_fourth -4.4852782 *r_fifth;

    for(int i=0; i<c.size();i++) {
        Result[i] = Approx_left[i] *(1-ponderation) + Approx_right[i]*ponderation;
    }


};
/*! \brief Multiplication between the transpose of the second tensor and the first one. Just used for square tensors, e.g., ndim x ndim XXX ndim x ndim
  @param[in] A First tensor
  @param[in] m Number of rows first tensor
  @param[in] n Number of columns first tensor
  @param[in] B Second tensor
  @param[in] l Number of rows second tensor
  @param[out] C Result
 */
extern void multSTensor3SecondTranspose(const vector<double> &A, int m, int n, const vector<double> &B, int l, vector<double> &C) {
    if (n != l) {
        cout << "ERROR: Wrong matrix multiplication" << endl;
        exit(-1);
    }
    setAll(C,0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < l; k++) {
                C[i * m + j] += A[i * m + k] * B[j * m + k];
            }
        }
    }
}

/*! \brief multiplication of a matrix with a scalar
  @param[inout] A Matrix
  @param[in] m Number of rows
  @param[in] n Number of columns
  @param[in] s Scalar number for multiplication
 */
extern void multiSTensorScalar(vector<double> &A, int m, int n, double s) {

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            A[i * m + j] = A[i * m + j] * s;
        }
    }

}


/*! \brief multiplication of a random matrix with a scalar
  @param[inout] A Matrix
  @param[in] m Number of rows
  @param[in] n Number of columns
  @param[in] s Scalar number for multiplication
  @param[in] nstoch Stochastic dimension
 */
extern void multiSTensorScalarStochastic(vector<double> &A, int m, int n, double s, int nstoch) {

    for (int p = 0; p < nstoch; p++) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[p * m * m + i * m + j] = A[p * m * m + i * m + j] * s;
            }
        }
    }
}

/*! \brief multiplication of a third (random) order tensor with a random variable
  @param[inout] A random Matrix
  @param[in] m Number of rows
  @param[in] n Number of columns
  @param[in] s random variable
  @param[in] nstoch Stochastic dimension
  @param[in] C Third order tensor that computes product.
  @param[out] D Result
 */
extern void
multiSTensorRandomVariableStochastic(const vector<double> &A, int m, int n, const vector<double>& s, int nstoch, const vector<double>& C,
                                     vector<double> &D) {
    setAll(D,0);
    for (int p = 0; p < nstoch; p++) {
        for (int k = 0; k < nstoch; k++) {
            for (int l = 0; l < nstoch; l++) {
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        D[p * m * m + i * m + j] +=
                                A[k * m * m + i * m + j] * s[l] * C[p * nstoch * nstoch + k * nstoch + l];
                    }
                }
            }
        }
    }
}

/*! \brief Determinant of a random tensor
  @param[in] A Random Matrix
  @param[in] ndim spatial dimension
  @param[in] nstoch stochastic dimension
  @param[inout] D Stochastic Jacobian
 */

extern void DeterminantStochasticTensor(const vector<double> &A, int ndim, int nstoch,
                                     vector<double> &D, const vector<double> &C) {
setAll(D,0);
    if(ndim ==2){
        vector<double> a(nstoch,0);
        vector<double> b(nstoch,0);
        vector<double> c(nstoch,0);
        vector<double> d(nstoch,0);
        // The determinant is ad-bc
        vector<double> bc(nstoch,0);
        vector<double> ad(nstoch,0);
        for(int i=0; i<nstoch; i++){
            a[i] = A[i*ndim*ndim + 0 * ndim + 0];
            b[i] = A[i*ndim*ndim + 0 * ndim + 1];
            c[i] = A[i*ndim*ndim + 1 * ndim + 0];
            d[i] = A[i*ndim*ndim + 1 * ndim + 1];
        }
        multiRandomVariableRandomVariableStochastic(a, d, C, ad);
        multiRandomVariableRandomVariableStochastic(b, c, C, bc);
        for(int i=0; i<nstoch; i++){
            D[i] = ad[i] - bc[i];
        }

    }

    if(ndim ==3){
        vector<double> a(nstoch,0);
        vector<double> b(nstoch,0);
        vector<double> c(nstoch,0);
        vector<double> d(nstoch,0);
        vector<double> e(nstoch,0);
        vector<double> f(nstoch,0);
        vector<double> g(nstoch,0);
        vector<double> h(nstoch,0);
        vector<double> I(nstoch,0);
        // The determinant is a
        vector<double> ei(nstoch,0);
        vector<double> fh(nstoch,0);
        vector<double> bi(nstoch,0);
        vector<double> ch(nstoch,0);
        vector<double> bf(nstoch,0);
        vector<double> ce(nstoch,0);
        // Another set of vector to create
        vector<double> ei_minus_fh(nstoch,0);
        vector<double> bi_minus_ch(nstoch,0);
        vector<double> bf_minus_ce(nstoch,0);
        //
        vector<double> a_times_ei_minus_fh(nstoch,0);
        vector<double> d_times_bi_minus_ch(nstoch,0);
        vector<double> g_times_bf_minus_ce(nstoch,0);
        for(int i=0; i<nstoch; i++){
            a[i] = A[i*ndim*ndim + 0 * ndim + 0];
            b[i] = A[i*ndim*ndim + 0 * ndim + 1];
            c[i] = A[i*ndim*ndim + 0 * ndim + 2];
            d[i] = A[i*ndim*ndim + 1 * ndim + 0];
            e[i] = A[i*ndim*ndim + 1 * ndim + 1];
            f[i] = A[i*ndim*ndim + 1 * ndim + 2];
            g[i] = A[i*ndim*ndim + 2 * ndim + 0];
            h[i] = A[i*ndim*ndim + 2 * ndim + 1];
            I[i] = A[i*ndim*ndim + 2 * ndim + 2];
        }
        multiRandomVariableRandomVariableStochastic(e, I, C, ei);
        multiRandomVariableRandomVariableStochastic(f, h, C, fh);
        multiRandomVariableRandomVariableStochastic(b, I, C, bi);
        multiRandomVariableRandomVariableStochastic(c, h, C, ch);
        multiRandomVariableRandomVariableStochastic(b, f, C, bf);
        multiRandomVariableRandomVariableStochastic(c, e, C, ce);
        for(int i=0; i<nstoch; i++){
            ei_minus_fh[i] = ei[i] - fh[i];
            bi_minus_ch[i] = bi[i] - ch[i];
            bf_minus_ce[i] = bf[i] - ce[i];
        }
        multiRandomVariableRandomVariableStochastic(a, ei_minus_fh, C, a_times_ei_minus_fh);
        multiRandomVariableRandomVariableStochastic(d, bi_minus_ch, C, d_times_bi_minus_ch);
        multiRandomVariableRandomVariableStochastic(g, bf_minus_ce, C, g_times_bf_minus_ce);
        for(int i=0; i<nstoch; i++){
            D[i] = a_times_ei_minus_fh[i] - d_times_bi_minus_ch[i] + g_times_bf_minus_ce[i] ;
        }
    }



}

/*! \brief multiplication of a third order tensor with a random variable
  @param[inout] A Matrix
  @param[in] m Number of rows
  @param[in] n Number of columns
  @param[in] s random variable
  @param[in] nstoch Stochastic dimension
  @param[in] C Third order tensor that computes product.
  @param[out] D Result
 */
extern void
multiSVectorRandomVariableStochastic(const vector<double> &A, int m, const vector<double>& s, int nstoch, const vector<double>& C,
                                     vector<double> &D) {

    for (int p = 0; p < nstoch; p++) {
        for (int k = 0; k < nstoch; k++) {
            for (int l = 0; l < nstoch; l++) {
                for (int i = 0; i < m; i++) {
                        D[p * m + i ] +=
                                A[k * m + i ] * s[l] * C[p * nstoch * nstoch + k * nstoch + l];

                }
            }
        }
    }
}

/*! \brief multiplication of two random variables
  @param[inout] A random variable
  @param[in] B random variable
  @param[in] C Third order tensor that computes product.
  @param[out] D Result
 */

extern void
multiRandomVariableRandomVariableStochastic(const vector<double> &A, const vector<double> &B, const vector<double>& C, vector<double> &D) {
// 	multiply A by B and store the result in D
    int stochastic_dim = A.size();
    for(int i=0; i<stochastic_dim; i++){
        D[i] = 0; // In case D is not empty
    }
    for (int k = 0; k < stochastic_dim; k++) {
        for (int i = 0; i < stochastic_dim; i++) {
            for (int j = 0; j < stochastic_dim; j++) {
                D[k] += C[i + j * stochastic_dim + k * stochastic_dim * stochastic_dim] * A[i] * B[j];
            }
        }
    }
}

/*! \brief Division of two random variables
  @param[inout] A random variable
  @param[in] B random variable
  @param[in] C Third order tensor that computes product.
  @param[out] D Result
 */

extern void DivisionRandomVariableRandomVariableStochastic(const vector<double> &A, const vector<double> &B, const vector<double>& C,
                                                           vector<double> &D) {
// Compute A/B and store the result in D
    int stochastic_dim = A.size();
// build the matrix sum_i BiCijk
    classMatrix *KKp;
    KKp = new classCompactMatrix(stochastic_dim, stochastic_dim, "COMPACT");
    double val = 0;
    for (int j = 0; j < stochastic_dim; j++) {
        for (int k = 0; k < stochastic_dim; k++) {
            for (int i = 0; i < stochastic_dim; i++) {
                val = val + B[i] * C[i + stochastic_dim * k + stochastic_dim * stochastic_dim * j];
            }
            KKp->addIJ(j, k, val);
            val = 0;
        }
    }
    KKp->assembly();
    bool flagConver = true;
    linearSystemSequential(KKp, A, D, flagConver);
    delete KKp;

}

/*! \brief Division of two random variables
  @param[inout] A random variable
  @param[in] B random variable
  @param[in] C Third order tensor that computes product.
  @param[out] D Result
 */

extern void InverseTensorStochastic(const vector<double> &A, int ndim, int nstoch,
                                                           const vector<double> &C, vector<double> &D) {
// Compute A/B and store the result in D
    setAll(D,0);
    int total_size = ndim*ndim*nstoch;
    vector<double> Identity(total_size,0);
    vector<double> Delta(ndim*ndim,0);
    for(int i=0; i<ndim; i++){
        Identity[i*ndim + i] = 1;
    }
    for(int i=0; i<ndim; i++){
        Delta[i*ndim + i] = 1;
    }
// build the matrix sum_i BiCijk
    classMatrix *KKp;
    KKp = new classCompactMatrix(total_size, total_size, "COMPACT");
    for (int i = 0; i < nstoch; i++) {
        for (int k = 0; k < nstoch; k++) {
            for (int l = 0; l < nstoch; l++) {
                for (int m = 0; m < ndim; m++) {
                    for (int n = 0; n < ndim; n++) {
                        for (int o = 0; o < ndim; o++) {
                            for (int p = 0; p < ndim; p++) {
                                KKp->addIJ(p + o * ndim + k * ndim * ndim, n + m * ndim + l * ndim * ndim,
                                           Delta[p + ndim * n] * A[i * ndim * ndim + o * ndim + m] *
                                           C[i + nstoch * l + nstoch * nstoch * k]);
                            }
                        }
                    }
                }
            }
        }
    }
    KKp->assembly();
    bool flagConver = true;
    linearSystemSequential(KKp, Identity, D, flagConver);
    delete KKp;
}

/*! \brief Multiplication of a random variable by a scalar plus another scalar (a*A +b)
  @param[inout] A random variable
  @param[in] a scalar
  @param[in] b scalar
 */

extern void linearCombinationVariableStochastic(vector<double> &A, double a, double b) {
    // cmpute A*a + b
    MultiplicationRandomVariableScalar(A, a);
    A[0] = A[0] + b;


}

/*! \brief Multiplication of a random variable by a scalar (a*A)
  @param[inout] A random variable
  @param[in] scalar scalar
  @param[in] b scalar
 */
extern void MultiplicationRandomVariableScalar(vector<double> &A, double scalar) {
    int stochastic_dim = A.size();
    for (int i = 0; i < stochastic_dim; i++) {
        A[i] = A[i] * scalar;
    }


}


/*! \brief This function a linear system of equations Ax=B in sequential. For more information about what this function does, please go to PETSC documentation
  @param[in] KK2 Matrix on the left hand side
  @param[in] rhs Right hand side vector
  @param[out] xx2 Array solution
  @param[out] flagConver Flag to check the convergence
*/
void linearSystemSequential(classMatrix *&KK2, const vector<double> &rhs, vector<double> &xx2, bool &flagConver) {

    Vec x, b;      /* approx solution, RHS, exact solution */
    Mat A;            /* linear system matrix */
    KSP ksp;         /* linear solver context */
    PC pc;           /* preconditioner context */
    // PetscReal      norm;  /* norm of solution error */
    PetscInt n = 10, its = 1000;
    PetscBool nonzeroguess;// = PETSC_FALSE;
    nonzeroguess = PETSC_FALSE;
    n = rhs.size();
    PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);
    PetscOptionsGetBool(NULL, NULL, "-nonzero_guess", &nonzeroguess, NULL);
    VecCreateSeq(PETSC_COMM_SELF, n, &x);
    PetscObjectSetName((PetscObject) x, "Solution");
    VecSetFromOptions(x);
    VecCreateSeq(PETSC_COMM_SELF, n, &b);
    PetscObjectSetName((PetscObject) b, "Righ hand side");
    VecSetFromOptions(b);

    int ind[n];
    for (int i = 0; i < n; i++) {
        ind[i] = i;
    }

    VecSetValues(b, n, ind, &rhs[0], INSERT_VALUES);
    //PetscInt start,end;
    //VecGetOwnershipRange(b,&start,&end);
    PetscScalar *array, *array3;
    VecGetArray(b, &array);
    KK2->getPetscMatSeq(A);
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    //   VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    //PCSetType(pc,PCJACOBI);
    KSPSetType(ksp, KSPGMRES);
    //KSPSetType(ksp,KSPCG);
    KSPSetTolerances(ksp, 1E-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    //VecView(x,PETSC_VIEWER_STDOUT_WORLD);
    KSPSolve(ksp, b, x);
    //VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    VecGetArray(x, &array3);
    PetscMemcpy(&xx2[0], array3, n * sizeof(double));
    KSPGetIterationNumber(ksp, &its);
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    if (reason < 0) {
        if (reason == -11) {
            PetscPrintf(PETSC_COMM_WORLD, "PETSC has failed. New attempt with Jacobi preconditioner\n");

            VecDestroy(&x);
            KSPDestroy(&ksp);

            VecCreate(PETSC_COMM_WORLD, &x);
            PetscObjectSetName((PetscObject) x, "Solution");
            VecSetSizes(x, PETSC_DECIDE, n);
            VecSetFromOptions(x);

            KSPCreate(PETSC_COMM_WORLD, &ksp);
            KSPSetOperators(ksp, A, A);
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCJACOBI);
            KSPSolve(ksp, b, x);
            VecGetArray(x, &array3);
            PetscMemcpy(&xx2[0], array3, n * sizeof(double));
            KSPGetIterationNumber(ksp, &its);
            KSPGetConvergedReason(ksp, &reason);
            if (reason < 0) {
                if (reason == -11) {
                    PetscPrintf(PETSC_COMM_WORLD, "PETSC has failed. New attempt without any preconditioner\n");

                    VecDestroy(&x);
                    KSPDestroy(&ksp);

                    VecCreate(PETSC_COMM_WORLD, &x);
                    PetscObjectSetName((PetscObject) x, "Solution");
                    VecSetSizes(x, PETSC_DECIDE, n);
                    VecSetFromOptions(x);

                    KSPCreate(PETSC_COMM_WORLD, &ksp);
                    KSPSetOperators(ksp, A, A);
                    KSPGetPC(ksp, &pc);
                    PCSetType(pc, PCNONE);
                    KSPSolve(ksp, b, x);
                    VecGetArray(x, &array3);
                    PetscMemcpy(&xx2[0], array3, n * sizeof(double));
                    KSPGetIterationNumber(ksp, &its);
                    KSPGetConvergedReason(ksp, &reason);
                    if (reason < 0) {
                        PetscPrintf(PETSC_COMM_WORLD, "PETSC has failed\n");
                        PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
                        ERROR("The linear system of equations did not converge");
                        flagConver = false;
                    } else {
                        PetscPrintf(PETSC_COMM_WORLD, "PETSC succeeded\n");
                    }
                } else {
                    PetscPrintf(PETSC_COMM_WORLD, "PETSC has failed\n");
                    PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
                    ERROR("The linear system of equations did not converge");
                    flagConver = false;
                }
            } else {
                PetscPrintf(PETSC_COMM_WORLD, "Jacobi preconditioner succeeded\n");
            }
        } else {
            PetscPrintf(PETSC_COMM_WORLD, "PETSC has failed\n");
            PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
            ERROR("The linear system of equations did not converge");
            flagConver = false;
        }
    } else {
        //PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D and iterations %D\n", reason, its);
    }
    VecDestroy(&x);
    VecDestroy(&b);
    KSPDestroy(&ksp);
}

/*! \brief Multiplication of a random variable by a scalar plus another scalar (a*A)
  @param[in] A random variable
  @param[out] scalar a boolean telling whether the variable is a scalar or not
 */
extern bool is_scalar(const vector<double> &A) {
    bool scalar = false;
    int stochastic_dim = A.size();
    int i = 1;
    while (A[i] == 0) {
        i = i + 1;
    }
    if (i == stochastic_dim - 1) {
        scalar = true;
    }
    return scalar;
}

/*! \brief Lump mass matrix. All values of a row are summed up and put at the diagonal of the matrix (here represented as a vector c)
  @param[in] A Mass matrix
  @param[in] numDof Total number of degrees of fredom
  @param[out] c Vector that is the lumped diagonal of the mass matrix
 */
extern void lumpMass(const vector<double> &A, int numDof, vector<double> &c) {
    for (int i = 0; i < numDof; i++) {
        double temp = 0;
        for (int j = 0; j < numDof; j++) {
            temp += A[i * numDof + j];
        }
        c[i] = temp;
    }
}

/*! \brief Calculate the inverse of a tensor in 2D and 3D. ndim x ndim
  @param[in] A Tensor
  @param[in] nn Size of the tensor
  @param[out] Ainv Result
 */
extern void inverse(const vector<double> &A, int nn, vector<double> &Ainv) {
    double det = 0;
    //To find cofactor and determinant
    det = determinantTensor(A, nn);
    //cout<<"The determinant of the matrix is "<<det<<endl;
    // if(det==0){
    //   det=1E-12;
    // }

    double idet = 1.0 / det;

    if (nn == 3) {
        Ainv[0] = idet * (A[nn + 1] * A[nn * 2 + 2] - A[nn * 1 + 2] * A[nn * 2 + 1]);
        Ainv[nn] = -idet * (A[nn * 1 + 0] * A[nn * 2 + 2] - A[nn * 1 + 2] * A[nn * 2 + 0]);
        Ainv[nn * 2] = idet * (A[nn * 1 + 0] * A[nn * 2 + 1] - A[nn * 1 + 1] * A[nn * 2 + 0]);

        Ainv[1] = -idet * (A[nn * 0 + 1] * A[nn * 2 + 2] - A[nn * 0 + 2] * A[nn * 2 + 1]);
        Ainv[nn + 1] = idet * (A[nn * 0 + 0] * A[nn * 2 + 2] - A[nn * 0 + 2] * A[nn * 2 + 0]);
        Ainv[2 * nn + 1] = -idet * (A[nn * 0 + 0] * A[nn * 2 + 1] - A[nn * 0 + 1] * A[nn * 2 + 0]);

        Ainv[2] = idet * (A[nn * 0 + 1] * A[nn * 1 + 2] - A[nn * 0 + 2] * A[nn * 1 + 1]);
        Ainv[nn + 2] = -idet * (A[nn * 0 + 0] * A[nn * 1 + 2] - A[nn * 0 + 2] * A[nn * 1 + 0]);
        Ainv[2 * nn + 2] = idet * (A[nn * 0 + 0] * A[nn * 1 + 1] - A[nn * 0 + 1] * A[nn * 1 + 0]);
    } else if (nn == 2) {
        Ainv[0] = idet * A[nn * 1 + 1];
        Ainv[nn] = -idet * A[nn * 1];

        Ainv[1] = -idet * A[1];
        Ainv[nn + 1] = idet * A[0];

    } else {
        cout << "ERROR: Not implemented" << endl;
        exit(-1);
    }

}

/*! \brief Auxiliar function for the calculation of the square root of a general matrix
 */
static void set_entry(double *A, int i, int j, double val, int N) {
    A[j * N + i] = val;
}

/*! \brief Auxiliar function for the calculation of the square root of a general matrix
 */
static double get_entry(const double *A, int i, int j, int N) {
    return A[j * N + i];
}
/*! \brief External function to for the calculation of the square root of a general matrix
 */
extern "C" void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
                        double *A, int *LDAp, double *VLp, double *VUp,
                        int *ILp, int *IUp, double *ABSTOLp, int *Mp,
                        double *W, double *Z, int *LDZp, int *ISUPPZ,
                        double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
                        int *INFOp);
/*! \brief External function to for the calculation of the square root of a general matrix
 */
extern "C" double dlamch_(char *CMACHp);

/*! \brief Auxiliar function for the calculation of the square root of a general matrix
 */
static int dsyevr(char JOBZ, char RANGE, char UPLO, int N,
                  double *A, int LDA, double VL, double VU,
                  int IL, int IU, double ABSTOL, int *M,
                  double *W, double *Z, int LDZ, int *ISUPPZ,
                  double *WORK, int LWORK, int *IWORK, int LIWORK) {

    int INFO;
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU,
            &IL, &IU, &ABSTOL, M, W, Z, &LDZ, ISUPPZ,
            WORK, &LWORK, IWORK, &LIWORK, &INFO);
    return INFO;
}

/*! \brief External function to for the calculation of the square root of a general matrix
 */
static double dlamch(char CMACH) {
    return dlamch_(&CMACH);
}


/*! \brief DGEEV computes for an N-by-N real nonsymmetric matrix A, the eigenvalues and, optionally, the left and/or right eigenvectors.. This function is external, and is part of LAPACK library.
The right eigenvector v(j) of A satisfies A * v(j) = lambda(j) * v(j) where lambda(j) is its eigenvalue.
The left eigenvector u(j) of A satisfies u(j)**H * A = lambda(j) * u(j)**H where u(j)**H denotes the conjugate transpose of u(j).
The computed eigenvectors are normalized to have Euclidean norm equal to 1 and largest component real.

@param[in] jobvl CHARACTER*1 = 'N': left eigenvectors of A are not computed; 'V': left eigenvectors of A are computed
@param[in] jobvr CHARACTER*1 = 'N': right eigenvectors of A are not computed; 'V': right eigenvectors of A are computed
@param[in] n INTEGER = The order of the matrix A. N >= 0
@param[inout] a DOUBLE PRECISION array, dimension (LDA,N) = On entry, the N-by-N matrix A. On exit, A has been overwritten.
@param[in]  lda INTEGER = The leading dimension of the array A.  LDA >= max(1,N)
@param[out] wr DOUBLE PRECISION array, dimension (N)
@param[out] wi DOUBLE PRECISION array, dimension (N)
WR and WI contain the real and imaginary parts, respectively, of the computed eigenvalues.  Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.
@param[out] vl DOUBLE PRECISION array, dimension (LDVL,N)
If JOBVL = 'V', the left eigenvectors u(j) are stored one after another in the columns of VL, in the same order as their eigenvalues.
If JOBVL = 'N', VL is not referenced. If the j-th eigenvalue is real, then u(j) = VL(:,j), the j-th column of VL.
If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and u(j+1) = VL(:,j) - i*VL(:,j+1)
@param[in] ldvl INTEGER. The leading dimension of the array VL.  LDVL >= 1; if JOBVL = 'V', LDVL >= N.
@param[out] vr DOUBLE PRECISION array, dimension (LDVR,N)
If JOBVR = 'V', the right eigenvectors v(j) are stored one after another in the columns of VR, in the same order as their eigenvalues.
If JOBVR = 'N', VR is not referenced.
If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1).
@param[in] ldvr INTEGER
The leading dimension of the array VR.  LDVR >= 1; if
JOBVR = 'V', LDVR >= N.
@param[out] work DOUBLE PRECISION array, dimension (MAX(1,LWORK))
On exit, if info = 0, WORK(1) returns the optimal LWORK.
@param[in] lwork   (input) INTEGER
The dimension of the array WORK.  LWORK >= max(1,3*N), and if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good  performance, LWORK must generally be larger.
If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
@param[out] info    (output) INTEGER = 0:  successful exit, < 0:  if INFO = -i, the i-th argument had an illegal value, > 0:  if INFO = i, the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed; elements i+1:N of WR and WI contain eigenvalues which have converged.
*/
extern "C" void
dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr,
       int *ldvr, double *work, int *lwork, int *info);


static void swap(double *a, int inca, double *b, int incb, int n) {
    double tmp;
    for (int i = 0; i < n; i++, a += inca, b += incb) {
        tmp = (*a);
        (*a) = (*b);
        (*b) = tmp;
    }
}

/*! \brief This function sorts the eigenvalues/vectors in ascending order according to their real part. Warning: this will screw up the ordering if we have complex eigenvalues. It is coupled to dgeev_ function
  @param[in] n The order of the matrix A. N >= 0
  @param[inout] wr DOUBLE PRECISION array, dimension (N)
  @param[inout] wi DOUBLE PRECISION array, dimension (N)
  WR and WI contain the real and imaginary parts, respectively, of the computed eigenvalues.  Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.
  @param[out] VL DOUBLE PRECISION array, dimension (LDVL,N) the left eigenvectors v(j) are stored one after another in the columns of VR, in the same order as their eigenvalues.
  @param[out] VR DOUBLE PRECISION array, dimension (LDVR,N) the right eigenvectors v(j) are stored one after another in the columns of VR, in the same order as their eigenvalues.
 */
static void eigenSort(int n, double *wr, double *wi, double *VL, double *VR) {
    // Sort
    for (int i = 0; i < n - 1; i++) {
        int k = i;
        double ek = wr[i];
        // search for something to swap
        for (int j = i + 1; j < n; j++) {
            const double ej = wr[j];
            if (ej < ek) {
                k = j;
                ek = ej;
            }
        }
        if (k != i) {
            swap(&wr[i], 1, &wr[k], 1, 1);
            swap(&wi[i], 1, &wi[k], 1, 1);
            swap(&VL[n * i], 1, &VL[n * k], 1, n);
            swap(&VR[n * i], 1, &VR[n * k], 1, n);
        }
    }
};

/*! \brief This function prepares the data to call the dgeev_ function and calculate the eigen values and eigen vectors
  @param[in] F Tensor to calculate its eigen values
  @param[in] ndim Dimension of the tensor
  @param[out] WR DOUBLE PRECISION array, dimension (N)
  @param[out] WI DOUBLE PRECISION array, dimension (N)
  WR and WI contain the real and imaginary parts, respectively, of the computed eigenvalues.  Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.
  @param[out] VL DOUBLE PRECISION array, dimension (LDVL,N) the left eigenvectors v(j) are stored one after another in the columns of VR, in the same order as their eigenvalues.
  @param[out] VR DOUBLE PRECISION array, dimension (LDVR,N) the right eigenvectors v(j) are stored one after another in the columns of VR, in the same order as their eigenvalues.
  @return Successfull calculation
 */
extern bool
eigen(vector<double>& F, int ndim, vector<double> &WR, vector<double> &WI, vector<double> &VL, vector<double> &VR) {

    int N = ndim, info;
    int lwork = 10 * N;
    double *work = new double[lwork];

    dgeev_("V", "V", &N, &F[0], &N, &WR[0], &WI[0], &VL[0], &N, &VR[0], &N, work, &lwork, &info);

    delete[] work;

    if (info > 0) {
        ERROR("QR Algorithm failed to compute all the eigenvalues %i %i", info, info);
        exit(0);
    } else if (info < 0) {
        ERROR("Wrong %d-th argument in eig", info);
        exit(0);
    }
    eigenSort(N, &WR[0], &WI[0], &VL[0], &VR[0]);

    return true;
};


extern void ConvertToFortran(char *fstring, std::size_t fstring_len, const char *cstring) {
    std::size_t inlen = strlen(cstring);
    std::size_t cpylen = std::min(inlen, fstring_len);
    std::copy(cstring, cstring + cpylen, fstring);
    std::fill(fstring + cpylen, fstring + fstring_len, ' ');
};

extern void ConvertToFortran2Darray(double *array, int files, int columns, double *farray) {
    for (int i = 0; i < files; i++) {
        for (int j = 0; j < columns; j++) {
            farray[j * files + i] = array[i * columns + j];
        }
    }
};


vector<double> symSTensor3(const vector<double> &a, int ndim) {
    vector<double> syms_a(ndim * ndim, 0.);
    vector<double> trans_a(ndim * ndim, 0.);
    transpose(a, ndim, trans_a);
    sumTensorTensor(a, trans_a, syms_a, ndim, 1.);
    multiSTensorScalar(syms_a, ndim, ndim, 0.5);
    return syms_a;
}

vector<double> vectorOperator(double sign, const vector<double>& a, const vector<double>& b) {
    vector<double> c(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        c[i] = a[i] + sign * b[i];
    }
    return c;
}

int indexofSmallestElement(const vector<double>& array) {
    int index = 0;

    for (int i = 1; i < array.size(); i++) {
        if (array[i] < array[index]) {
            index = i;
        }
    }

    return index;
}

void getCofactor(const vector<vector<double> > &A, vector<vector<double> > &temp, int p, int q, int n) {
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q) {
                temp[i][j++] = A[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
double determinant(const vector<vector<double> > &A, int n, int N) {
    double D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];

    vector<vector<double>> temp; // To store cofactors
    vector<double> inner;
    inner.resize(N, 0);
    for (int i = 0; i < N; i++) {
        temp.push_back(inner);
    }

    int sign = 1;  // To store sign multiplier

    // Iterate for each element of first row
    for (int f = 0; f < n; f++) {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1, N);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(const vector<vector<double> > &A, vector<vector<double>> &adj, int N) {
    if (N == 1) {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1;
    vector<vector<double>> temp;
    vector<double> inner;
    inner.resize(N, 0);
    for (int i = 0; i < N; i++) {
        temp.push_back(inner);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, N - 1, N));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(const vector<vector<double> > &A, vector<vector<double> > &inverse, int N) {
    // Find determinant of A[][]
    double det = determinant(A, N, N);
    if (det == 0) {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    vector<vector<double>> adj;
    vector<double> inner;
    inner.resize(N, 0);
    for (int i = 0; i < N; i++) {
        adj.push_back(inner);
    }
    adjoint(A, adj, N);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            inverse[i][j] = adj[i][j] / float(det);

    return true;
}

double dotProduct(const vector<double> &a, const vector<double> &b, int ndim) {
    double c = 0;
    for (int i = 0; i < ndim; i++) {
        c += a[i] * b[i];
    }
    return c;
}

void vectorArithmetic(const vector<double> &a, const vector<double> &b, vector<double> &c, int ndim, bool add) {
    int f;
    if (add) {
        f = 1;
    } else {
        f = -1;
    }
    for (int i = 0; i < ndim; i++) {
        c[i] = a[i] + f * b[i];
    }
}

/*! \brief compute all vectors of non negative integers of size n whose sum is equal to k
  @param[in] n
  @param[in] k
  @param[in] rowindex
  @param[inout] double_poly
 */

void sum_nonnegative_integers_k(int n, int k, vector<vector<int>> &double_poly, int &row_index) {
    if (n == 1) {
        double_poly[row_index][0] = k;
        row_index = row_index + 1;
    } else {

// fill n + k - 1 boxes with n - 1 balls
        vector<int> box(n + k - 1, 0);
// initialise the position of the balls
        for (int i = 0; i < n - 1; i++) {
            box[i] = 1;
        }
        vector<int> test(n, 0);
        bool all_right = balls_right(box, n - 1);

        vector<int> pos = get_position(box);

        for (int i = 0; i < n; i++) {
            if (i == 0) {
                test[0] = pos[0];
            } else if (i == n - 1) {
                test[n - 1] = n + k - 2 - pos[n - 2];
            } else {
                test[i] = pos[i] - pos[i - 1] - 1;
            }
        }
        std::copy(test.begin(), test.begin() + test.size(), double_poly[row_index].begin());

        row_index = row_index + 1;
        int previous_shift_ball = find_index_shifted(box, n - 1);
        while (!all_right) {
            vector<int> test(n, 0);
            int new_shift_ball = find_index_shifted(box, n - 1);
            int incr = 1;
            if (new_shift_ball == previous_shift_ball) {
                box[pos[new_shift_ball]] = 0;
                box[pos[new_shift_ball] + 1] = 1;
            } else {
                box[pos[new_shift_ball]] = 0;
                box[pos[new_shift_ball] + 1] = 1;
                for (int i = new_shift_ball + 1; i < n - 1; i++) {
                    box[pos[i]] = 0;
                    box[pos[new_shift_ball] + 1 + incr] = 1;
                    incr = incr + 1;
                }
            }
            pos = get_position(box);
            for (int i = 0; i < n; i++) {
                if (i == 0) {
                    test[0] = pos[0];
                } else if (i == n - 1) {
                    test[n - 1] = n + k - 2 - pos[n - 2];
                } else {
                    test[i] = pos[i] - pos[i - 1] - 1;
                }
            }

            std::copy(test.begin(), test.end(), double_poly[row_index].begin());
            row_index = row_index + 1;
            previous_shift_ball = new_shift_ball;
            all_right = balls_right(box, n - 1);

        }
    }
}

/*! \brief return binomial coefficient
  @param[in] n
  @param[in] k
  @param[out] value of the binomial coefficient
 */

int binomialCoefficients(int n, int k) {
    if (k == 0 || k == n) {
        return 1;
    } else {
        return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);
    }
}

/*! \brief tell if one ball of the box is at the far right or not
  @param[in] box
  @param[in] n
  @param[out] return a boolean telling wether the ball is at the right of the box or not
 */

bool balls_right(const vector<int>& box, int n) {
    bool flag_allright = false;
    int c = 0;

    for (int j = 1; j <= n; j++) {
        if (box[box.size() - j] == 1) {
            c = c + 1;
        }
    }
    if (c == n) {
        flag_allright = true;
    }
    return flag_allright;
}

/*! \brief return the index of the vector whose component is equal to 1
  @param[in] vector containing values of 0 and 1
  @param[inout] index of component equals to 1
 */
vector<int> get_position(const vector<int>& box) {
    vector<int> pos;

    for (int i = 0; i < box.size(); i++) {
        if (box[i] == 1) {
            pos.push_back(i);
        }
    }
return pos;
}

/*! \brief map of index of the univariate polynomials to the global one
  @param[in] table poly coefficient of the polynomials
  @param[in] j index of the random variable
  @param[out] return mapping
 */

vector<int> local2global(const vector<vector<int> >& table_poly, int j) {

    int size_1 = table_poly.size();
    int size_2 = table_poly[0].size();
    vector<int> localtoglobal;
    for (int k = 1; k < size_1; k++) {
        bool flag_cross = true;
        for (int l = 0; l < size_2; l++) {
            if (!(l == j)) {
                if (!(table_poly[k][l] == 0)) {
                    flag_cross = false;
                }
            }
        }
        if (flag_cross) {
            localtoglobal.push_back(k);
        }
    }
    return localtoglobal;
}

/*! \brief analog of the matlab linspace function
  @param[in] a start point
  @param[in] b end point
  @param[in] num number of points desired
  @param[out] return result
 */
vector<double> linspace(double a, double b, int num) {

    vector<double> result(num, 0);

    for (int i = 0; i < num; i++) {
        result[i] = a + i * (b - a) / (num - 1);
    }

    return result;

}

int find_index_shifted(const vector<int>& box, int n) {
    bool can_shift_right = false;
    int shift_index = 0;
    if (box[box.size() - 1] == 0) {
        shift_index = n - 1;
        can_shift_right = true;
    }
    int s = 1;
    while (!can_shift_right) {
        if (box[box.size() - s ] == 0) {
            shift_index = n  - s;
            can_shift_right = true;
        }
        s = s + 1;
    }
    return shift_index;
}

/*! \brief Creates a mapping between the CDF (comprised between 0 and 1) and the value of random gaussian variable
  @param[in] p value of the gaussian variable
  @param[out] return result
 */
double QuantileNormal(double p) {
// 	https://books.google.co.uk/books?hl=fr&lr=&id=RSC002LFbY0C&oi=fnd&pg=PA1&ots=Ghrwe2WiiS&sig=GhIwST5vX8XR3uIqL3SyEkxrPOc&redir_esc=y#v=onepage&q&f=false page 933, rational approxmiation of quantile func
    double quantile;
    double a_0 = 2.30753;
    double a_1 = 0.27061;
    double b_0 = 0.99229;
    double b_1 = 0.04481;
    double t;
    if (p > 0.5) {
        t = sqrt(-2 * log(1 - p));
        quantile = t - (a_0 + a_1 * t) / (1 + b_0 * t + b_1 * t * t);
    } else {
        t = sqrt(-2 * log(p));
        quantile = -t + (a_0 + a_1 * t) / (1 + b_0 * t + b_1 * t * t);
    }
    return quantile;
}

/*! \brief Creates a mapping between the CDF (comprised between 0 and 1) and the value of random uniform variable
  @param[in] p value of the gaussian variable
  @param[out] return result
 */
double QuantileUniform(double p) {
    double quantile;
    quantile = 2*p -1;
    return quantile;
}

extern void isInCircle(bool &A, const vector<double>& n0, const vector<double>& center, const vector<double>& coordinates, double radius){
    int ndim = center.size();
    vector<double> difference(ndim,0);
    vector<double> projection(ndim,0);
    vector<double> newcenter(ndim,0);
    double dist = 0;
    for(int i=0; i<ndim; i++){
        difference[i] = center[i] - coordinates [i];
        projection[i] = difference[i] * n0[i];
        newcenter[i] = center[i] - projection[i];
    }

    dist = distance(newcenter, coordinates);

    if(dist<radius){
        A = true;
    }



}

double distance (const vector<double>& center, const vector<double>& coordinates) {
    int ndim = center.size();
    double distance = 0;
    double sumsquare = 0;

    for (int i=0;i<ndim;i++){
        sumsquare += pow((center[i]-coordinates[i]),2);
    }
    distance = sqrt(sumsquare);
return distance;
}

/*! \brief Evaluate 2^(j/2)*Haar(2^j*p-k)
  @param[in] j
  @param[in] k
  @param[in] p quantile value comprised between 0 and 1
  @param[out] return result
 */
double evaluate_haar(int j, int k, double p) {
    int value;
    value = evaluate_wavelet(pow(2, j) * p - k);
    return value * pow(2, j / 2.0);
}

/*! \brief Evaluate Haar(y)
  @param[in] p quantile value comprised between 0 and 1
  @param[out] return result
 */
int evaluate_wavelet(double y) {
    int value;
    if ((y >= 0) && (y < 0.5)) {
        value = 1;
    } else {
        if ((y >= 0.5) && (y < 1)) {
            value = -1;
        } else {
            value = 0;
        }
    }
    return value;
}


/*! \brief Evaluate mixed(y)
  @param[in] index of the function Psi
  @param[in]

 */

double evaluate_mixed(int i, int order, double gauss) {
    double result=0;
    if (i < order + 1) {
        switch (i) {
            case 0:
                result = 1;
                break;
            case 1:
                result = sqrt(2)*sqrt(1.5) * (2*gauss-1) ;
                break;
            case 2:
                result = sqrt(2)*sqrt(2.5) * 0.5*(3*pow(2,2*gauss-1)-1);
                break;
        }
    } else {
        i = i - order - 1;
        int remainder = i % (order + 1);
        i = i - remainder;
        i = i / (order + 1) +1;
        int j = floor(log(i) / log(2));
        int k = i - pow(2, j);

        gauss = pow(2, j) * gauss - k;

        if ((gauss < 0) || (gauss > 1)) {
            result = 0;
        } else {
                if(order==1) {
                    switch (remainder) {
                        case 0:
                            if (gauss < 0.5) {
                                result = sqrt(2) * pow(2, (float) j / 2) * pow(-1, 1 + remainder + order) * sqrt(1.5) *
                                         (-1 + -2 * (2 * gauss - 1));
                            } else {
                                result = sqrt(2) * pow(2, (float) j / 2) * sqrt(1.5) * (-1 + 2 * (2 * gauss - 1));
                            }
                            break;
                        case 1:
                            if (gauss < 0.5) {
                                result = sqrt(2) * pow(2, (float) j / 2) * pow(-1, 1 + remainder + order) * sqrt(0.5) *
                                         (-2 - 3 * (2 * gauss - 1));
                            } else {
                                result = sqrt(2) * pow(2, (float) j / 2) * sqrt(0.5) * (-2 + 3 * (2 * gauss - 1));
                            }
                    }
                }
                if(order==2) {
                    switch (remainder) {
                    case 0:
                        if (gauss < 0.5) {
                            result = sqrt(2) * pow(2, (float) j / 2) * pow(-1, 1 + remainder + order) * 1 / 3 *
                                     sqrt(0.5) *
                                     (1 + 24 * (2 * gauss - 1) + 30 * pow(2, 2 * gauss - 1));
                        } else {
                            result = sqrt(2) * pow(2, (float) j / 2) * 1 / 3 * sqrt(0.5) *
                                     (1 - 24 * (2 * gauss - 1) + 30 * pow(2, 2 * gauss - 1));
                        }
                        break;
                    case 1:
                        if (gauss < 0.5) {
                            result = sqrt(2) * pow(2, (float) j / 2) * pow(-1, 1 + remainder + order) * 1 / 2 *
                                     sqrt(1.5) *
                                     (3 + 16 * (2 * gauss - 1) + 15 * pow(2, 2 * gauss - 1));
                        } else {
                            result = sqrt(2) * pow(2, (float) j / 2) * 1 / 2 * sqrt(1.5) *
                                     (3 - 16 * (2 * gauss - 1) + 15 * pow(2, 2 * gauss - 1));
                        }
                        break;
                    case 2:
                        if (gauss < 0.5) {
                            result = sqrt(2) * pow(2, (float) j / 2) * pow(-1, 1 + remainder + order) * 1 / 3 *
                                     sqrt(2.5) *
                                     (4 + 15 * (2 * gauss - 1) + 12 * pow(2, 2 * gauss - 1));
                        } else {
                            result = sqrt(2) * pow(2, (float) j / 2) * 1 / 3 * sqrt(2.5) *
                                     (4 - 15 * (2 * gauss - 1) + 12 * pow(2, 2 * gauss - 1));

                        }
                    }
                }
        }
    }
return result;
}


/*! \brief Compute the result of a product between two polynomials A and B
  @param[in] A vector with the coefficients of the first polynomial
  @param[in] B vector with the coefficients of the second polynomial
  @param[out] return result
 */
vector<double> product(const vector<double>& A, const vector<double>& B) {
    int dim = A.size();
    vector<double> result(dim, 0);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
            result[i] += A[j] * B[i - j];
        }
    }
    return result;
}

int calculationFactorial(int i) {
    int half_i = i / 2;
    int result;
    double resultDouble = 1;
    if (i == 0) {
        result = 1;
    } else {
        for (int j = 0; j < half_i; j++) {
            resultDouble *= (float)(i - j) / (float)2;
        }
    }
    result = (int) resultDouble;

    return result;


}

double calculationMomentUniform(int i) {
    double resultDouble = 0;
    if(i==0){
        resultDouble = 1;
    } else {
        resultDouble = (1 - pow(-1,i+1))/(2*(i+1));
    }

    return resultDouble;
}

void findCommonInterval(int i, int j, int k, bool &flagCommonInterval, vector<int> &indices, vector<double> &interval, int &max){
    max = i;
    vector<double> range(6,0);
    if(max<j){
        max = j;
    }
    if (max<k){
        max = k;
    }
    int temp_j;
    int temp_k;
    if(i==0){
        indices[0] = 0;
        indices[1] =0;
        range[0] = 0;
        range[1] =1;
    }
    else {
        indices[0] = floor(log(i) / log(2));
        indices[1] = i - pow(2,indices[0]);
        range[0] = indices[1] * pow(2, -indices[0]);
        range[1] = (indices[1]+1) * pow(2, -indices[0]);
    }
    if(j==0){
        indices[2] = 0;
        indices[3] =0;
        range[2] = 0;
        range[3] =1;
    }
    else {
        indices[2] = floor(log(j) / log(2));
        indices[3] = j - pow(2,indices[2]);
        range[2] = indices[3] * pow(2, -indices[2]);
        range[3] = (indices[3]+1) * pow(2, -indices[2]);
    }
    if(k==0){
        indices[4] = 0;
        indices[5] =0;
        range[4] = 0;
        range[5] =1;
    }
    else {
        indices[4] = floor(log(k) / log(2));
        indices[5] = k - pow(2,indices[4]);
        range[4] = indices[5] * pow(2, -indices[4]);
        range[5] = (indices[5]+1) * pow(2, -indices[4]);
    }
    if (max==0){
        interval[0] =0;
        interval[1] =1;
    }else {
        temp_j = floor(log(max) / log(2));
        temp_k = max - pow(2,temp_j);

        interval[0] = temp_k * pow(2, -temp_j);
        interval[1] = (temp_k + 1) * pow(2, -temp_j);
    }
    if((range[0] <= interval[0])&&(range[2] <= interval[0])&&(range[4] <= interval[0])&&(range[1] >= interval[1])&&(range[3] >= interval[1])&&(range[5] >= interval[1])){
        flagCommonInterval = true;
    }
}
