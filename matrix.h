#ifndef __MATRIX_H__
#define __MATRIX_H__
#define __MATRIX__VERSION__ 1079
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<omp.h>
//#include<mpich/mpi.h> Resever

/*
Ver1.0.79

Matrix.h is a header which could help in calculate matrix.

General Inverse would be added in next version.

Debug environment: Windows Linux Subsystem with g++

*/

/*
The function is COLUMN FIRST, i.e. matrix

    0   1
    2   3

would be stored as {0,2,1,3}
*/

class Matrix
{
    /*
    Data: The size of matrix(storing in row and cl), where elements would be stored(matadd) and whether matrix 
    could be used(matstat)

    As PRIVATE member, you should always read or change these member by PUBLIC member function, unless the object
    has been constant.
    */
    private:
    int row;
    int cl;
    double* matadd;
    bool matstat;

    /*
    double inner(double*,double*,int) is a function that you should never try to use(that's the reason it is a 
    PRIVATE member function)
    */

    double inner(double* vec1,double* vec2,int len);

    /*
    Matrix scalemult(double) is a out-of-date function which only exists for several supports, its function has
    been replaced by overload of operator*, see the last part.
    */

    Matrix scalemult(double scale);

    /*
    void cp(double*,double*,int) is a function that you should never try to use(that's the reason it is a 
    PRIVATE member function)
    */

    void cp(double* source,double* target,int len);

    public:

    /*
    Default value when you creat an object: the matrix does not exsit (cause it is a 0x0 matrix), nowhere to store
    elements(the address is NULL) and the matrix could not be used.
    */

    Matrix():matstat(false),row(0),cl(0),matadd(NULL)
    {}

    /*
    Following member functions are initial functions
    */

    /*
    void setmat(int,int):function for initial, set row and column, allocate memory, set the state of matrix could be
    used. All the other function should be used after this function.

    What should be input should be row(1st parameter) and column(2nd parameter).

    Initial matrix is a zero matrix.

    To free the memory after use, see following function. Should not try to free memory by yourself. 
    */

    void setmat(int ro,int clo);

    void setmat(int dim);

    /*
    void writematman():Data input function, get elements from standard input stream. UI has been designed.
    */

    void writematman();

    /*
    void writematbyad(double*):Copy elements from a target address. Please check you have arrange the elements in
    right consequence.
    */

    void writematbyad(double* targetadd);

    /*
    Following member functions are information functions.
    */

    /*
    void getinfo(int*,int*,double**): Get information by address, a handful example to use the function:
    ...
    int r,c;
    double* add;
    M.getinfo(&r,&c,&add)   

    You will get default value if matrix is not initialized.
    */

    void getinfo(int* rowadd,int* cladd,double** getmatadd);

    /*
    void display(): Print the matrix on the standard output, UI HAS NOT BEEN DESIGNED, you need to add information
    before or after the function by yourself.
    */

    void display();
    /*
    int indtr(int,int),int* indtr(int,int*):index transfer function
    */
    int indtr(int rn,int cn);

    int* indtr(int ind,int* store);

    /*
    Following functions are free function.
    */

    /*
    void deletepre(): free the memory, reset the matrix to default.
    
    No member function should appear after deletpre.
    */

    void deletepre();

    /*
    Following parts in fact are overload or new defined operator -- =,transposition,+,-,*(scale multiple and 
    matrix multiple).

    A 1x1 matrix is equivalent to a scale when using operator*, which only means it would not cause error.
    */
    
    Matrix& operator =(const Matrix& Ma);

    double operator()(const int rn,const int cn) const;

    double& operator()(const int rn,const int cn);

    double operator()(int* cgt) const;

    double& operator()(int* cgt);

    double operator[](const int ind) const;
    
    double& operator[](const int ind);

    Matrix matT();

    Matrix operator+(const Matrix& Ma);

    friend Matrix operator+(const Matrix& Ma, const Matrix& Mb);

    Matrix operator-(const Matrix& Ma);

    friend Matrix operator-(const Matrix& Ma, const Matrix& Mb);

    Matrix operator-();

    friend Matrix operator *(double scale,const Matrix& Ma);

    friend Matrix operator *(const Matrix& Ma,double scale);

    Matrix operator *(const Matrix& Mb);

    friend Matrix operator*(const Matrix& Ma,const Matrix& Mb);

    bool iszero();

    bool issq();

    bool isunit();

    bool issym();

    bool isasym();

    bool operator==(const Matrix& Ma);

    bool operator!=(const Matrix& Ma);

    friend Matrix abs(const Matrix& Ma);

    friend Matrix getmax(Matrix& Ma);

    friend Matrix getmin(Matrix& Ma);

    friend Matrix diag(double* targets,int len);

    friend Matrix diag(const Matrix& Ma);

    friend Matrix Unit(int dim);

    friend Matrix pow(const Matrix& Ma,const int n);

    friend Matrix pow(const double a,const Matrix& Ma,int ord);

    friend Matrix Lie(const Matrix& MA,const Matrix& MB);

    friend Matrix Dij(int ind,int jnd, int row,int cl,double ele);

    friend Matrix Dij(int ind,int jnd, int dim,double ele);

    friend Matrix Dij(int* tgind, int row,int cl,double ele);

    friend Matrix Dij(int* tgind, int dim,double ele);

    friend Matrix Dij(int ind,int jnd, int* tgdim,double ele);

    friend Matrix Dij(int* tgind, int* tgdim,double ele);

    friend Matrix Pij(int ind,int jnd,int dim);

    friend Matrix Pij(int* tg,int dim);

    friend Matrix Randmat(int rn,int cn);

    friend Matrix Randmat(int dim);
};

Matrix abs(const Matrix& Ma);
Matrix getmax(Matrix& Ma);
Matrix getmin(Matrix& Ma);
Matrix diag(double* targets,int len);
Matrix diag(const Matrix& Ma);
Matrix Unit(int dim);
Matrix pow(const Matrix& Ma,const int n);
Matrix pow(const double a,const Matrix& Ma,int ord);
Matrix Lie(const Matrix& MA,const Matrix& MB);
Matrix Dij(int ind,int jnd, int row,int cl,double ele);
Matrix Dij(int ind,int jnd, int dim,double ele);
Matrix Dij(int* tgind, int row,int cl,double ele);
Matrix Dij(int* tgind, int dim,double ele);
Matrix Dij(int ind,int jnd, int* tgdim,double ele);
Matrix Dij(int* tgind, int* tgdim,double ele);
Matrix Pij(int ind,int jnd,int dim);
Matrix Pij(int* tg,int dim);
Matrix Randmat(int rn,int cn);
Matrix Randmat(int dim);
#endif