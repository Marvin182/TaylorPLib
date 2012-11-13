#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <typeinfo>
#include "CustomException.h"
// #include "esccolors.h"
// #include "ioutil.h"
// #include "TPolyn.h"

namespace LibMatrix {

	class __declspec(dllexport) Matrix
	{
	private:
			int _rows,											// The number of rows
				_cols,											// The number of columns
				_dimT;											// The dimension of the Taylor polynomials 
    		double **_data;										// The pointer to the allocated memory

    		void allocateMemory(bool initialize);
	    	void deallocateMemory();
			void copyFrom(const Matrix &m);

	public:
			//
			// Constructors, destructor
			//
			Matrix::Matrix();														// Default constructor
			Matrix::Matrix(int rows, int cols, bool initialize = true);				// Regular constructor
			Matrix::Matrix(int rows, int cols, int dimT, bool initialize = true);	// Regular constructor
			Matrix::Matrix(const Matrix &m);										// Copy constructor

			Matrix::Matrix(double* values, int rows = 2, int cols = 2);				// easy test construtor		

			// Matrix Matrix::redim(int rows, int cols);							// Reset the current object dimensions
			// Matrix Matrix::redim(int rows, int cols, int dimT);					// Reset the current object dimensions
			Matrix::~Matrix();													// Destructor

			//
			// Accessing properties
			//
			double** data() { return _data; }						// Returns the ptr to the alloc. memory
			int nrows() const { return _rows; }							// Returns the number of rows
			int ncols() const { return _cols; }							// Returns the number of columns
			int dimT() const { return _dimT; }							// Returns the dimension of the type T
			// void setnrows(int rows ) { _rows = rows; }			// Sets the number of rows
			// void setncols(int cols ) { _cols = cols; }			// Sets the number of columns
			// void setdimT(int dimT ) { _dimT = dimT; }			// Sets the dimension of T
			// void setdim(int rows, int cols ) {
			// 	_rows = rows; _cols = cols; }						// Sets numbers of rows and columns
			// void setdim(int rows, int cols, int dimT ) {
			// 	_rows = rows; _cols = cols; _dimT = dimT; }			// Sets numbers of rows and columns and the dimension of the Taylor polynomials

			//
			// Overloaded operators
			//
			double & Matrix::operator()(int i, int j);				// Element
			Matrix Matrix::operator=(const Matrix &m);				// Assignment operators
			bool Matrix::operator==(const Matrix &m);				// Comparison operators
			bool Matrix::operator!=(const Matrix &m);
			Matrix Matrix::operator+(const Matrix &m);				// Arithmetic operators
			Matrix Matrix::operator+=(const Matrix &m);
			Matrix Matrix::operator-(const Matrix &m);
			Matrix Matrix::operator-=(const Matrix &m);
			Matrix Matrix::operator-();
			Matrix Matrix::operator*(double alpha);
			Matrix Matrix::operator*=(double alpha);
			Matrix Matrix::operator*(const Matrix &m);				// simple matrix multiplication

			//
			// Especial matrix multiplications
			//
			void Matrix::mmCaABbC(double alpha, double beta, const Matrix &A, const Matrix &B);
			// void Matrix::bmmCaABbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
			// void Matrix::mmCasABbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
			// void Matrix::mmCaAsBbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
			// void Matrix::mmCaAUTBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
			void Matrix::mmCaAATbC(double alpha, double beta, const Matrix&A);
			void Matrix::mmCaATAbC(double alpha, double beta, const Matrix&A);
			// void Matrix::mmCaATBbC(double alpha, double beta, const Matrix&A, const Matrix&B);
			// void Matrix::mmCaATBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
			// void Matrix::mmCaABTbC(double alpha, double beta, const Matrix&A, const Matrix&B);
			// void Matrix::mmCaABTbC(int r, bool up, double alpha, double beta, const Matrix&A, const Matrix&B);
			// void Matrix::bmmCaABTbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
			// void Matrix::mmCaIBbC(double alpha, double beta, const Matrix&B);
			// void Matrix::mmCaIBbC(double alpha, double beta, int *piv, bool rows, const Matrix&B);
			// void Matrix::mmCaAIbC(double alpha, double beta, const Matrix&A);
			// void Matrix::mmCaAIbC(double alpha, double beta, Matrix &A, int *piv, bool rows);

			//
			// Solving systems of equations
			//
			// int Matrix::utsolve(Matrix &B);					// U X = B, back-substitution
			// int Matrix::utsolve(Matrix &B, Matrix &X, int *piv);// U X = B, back-substitution, pivoting
			// int Matrix::utsolve(T *b);						// U x = b, back-substitution
			// int Matrix::utxsolve(Matrix &B);					// X U = B, back-substitution
			// int Matrix::gsolve(Matrix &B);					// A X = B, Gaussian elimination
			// int Matrix::gsolve(T *b);							// A x = b, Gaussian elimination

			//
			// Other functions
			//
			// double Matrix::fnorm();								// Frobenius norm
			// int Matrix::tcfnorm(int nrTC, double *fn);		// vector of Frobenius norms
			// int Matrix::colnorm(double *c);
			// int Matrix::colnormdown(int pos, int *piv, double *c);
			// int Matrix::cpermutem(int *piv);
			// int Matrix::cpermutem(int *piv, bool trans);
			// int Matrix::rpermutem(int *piv);
			// Matrix Matrix::transpm();
			// int Matrix::shift(Matrix &M);
			// bool Matrix::isId();
			// bool Matrix::isId(double eps);
			// bool Matrix::isId(int m1, int m2, int n1, int n2, double eps);
			// bool Matrix::isZero();
			// bool Matrix::isZero(double eps);
			// int Matrix::set2Id();
			// int Matrix::set2Id(int m, int n);
			// int Matrix::set2Id(int m1, int m2, int n1, int n2);
			// int Matrix::set2zero();
			// int Matrix::set2zero(int m, int n);
			// int Matrix::set2zero(int m1, int m2, int n1, int n2);
			// int Matrix::set2val(double v);
			// int Matrix::set2val(int i, int j, double val);
			// int Matrix::trinvm(Matrix &Am);
			// int Matrix::trinvm(int r, Matrix &Am);
			// bool Matrix::mcompare(Matrix &B);
			// bool Matrix::mcompare(Matrix &B, double eps);
			// bool Matrix::mcompare(Matrix &B, int r, double eps);

			//
			// Print out functions
			//
			void Matrix::printm();
			void Matrix::printm(char *str);
			// void Matrix::fprintm(FILE * fn, char *str);
			// void Matrix::printm(char *str, const char *const color);
			// void Matrix::printdm(char *str);
			// void Matrix::fprintdm(FILE * fn, char *str);
			// void Matrix::printdm(char *str, const char *const color);
			// void Matrix::printdm(char *str, const char *const color, double eps);
			// void Matrix::fprintdm(FILE * fn, char *str, double eps);
			// void Matrix::printm(int *piv, char *str, const char *const color);
			// void Matrix::printdm(int *piv, char *str, const char *const color, double eps);
			// void Matrix::fprintm(FILE * fn, int *piv, char *str);
			// void Matrix::fprintdm(FILE * fn, int *piv, char *str, double eps);
			// void Matrix::printtpm(char *str, const char *const color, double eps);
			// void Matrix::fprinttpm(FILE * fn, char *str, double eps);
			// void Matrix::printtpm(int m1, int m2, int n1, int n2,
			// 		char *str, const char *const color, double eps);
			// void Matrix::fprinttpm(FILE * fn, int m1, int m2, int n1, int n2,
			// 		char *str, double eps);
			// void Matrix::printtpm(int *piv, char *str, const char *const color, double eps);
			// void Matrix::fprinttpm(FILE * fn, int *piv, char *str, double eps);
	};
}

#endif /* _MATRIX_H_ */
