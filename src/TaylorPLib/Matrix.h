#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <typeinfo>
#include <stdarg.h> // only used by the test constructor
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

    		void allocateDataMemory(bool initialize)
    			{ Matrix::allocateMemory(_data, _rows, _cols, initialize); }
    		static void allocateMemory(double **&data, int rows, int cols, bool initialize);
	    	void deallocateDataMemory()
	    		{ Matrix::deallocateMemory(_data, _rows, _cols); }
	    	static void deallocateMemory(double **&data, int rows, int cols);
			void copyFrom(const Matrix &m);

	public:
			//
			// Constructors, destructor
			//
			Matrix::Matrix();														// Default constructor
			Matrix::Matrix(int rows, int cols, bool initialize = true);				// Regular constructor
			Matrix::Matrix(int rows, int cols, int dimT, bool initialize = true);	// Regular constructor
			Matrix::Matrix(const Matrix &m);										// Copy constructor

			Matrix::Matrix(int rows, int cols, double *values);						// Test construtor		

			Matrix::~Matrix();														// Destructor

			//
			// Accessing properties
			//
			double** data() { return _data; }			// Returns the ptr to the alloc. memory
			int nrows() const { return _rows; }			// Returns the number of rows
			int ncols() const { return _cols; }			// Returns the number of columns
			int dimT() const { return _dimT; }			// Returns the dimension of the type T
			double get(int row, int col) const;			// Returns a single element from the matrix

			//
			// Overloaded operators
			//
			double & Matrix::operator()(int row, int col);			// Element
			Matrix Matrix::operator=(const Matrix &m);				// Assignment operators
			bool Matrix::operator==(const Matrix &m) const;			// Comparison operators
			bool Matrix::operator!=(const Matrix &m) const;
			Matrix Matrix::operator+(const Matrix &m) const;		// Arithmetic operators
			Matrix Matrix::operator+=(const Matrix &m);
			Matrix Matrix::operator-(const Matrix &m) const;
			Matrix Matrix::operator-=(const Matrix &m);
			Matrix Matrix::operator-();
			Matrix Matrix::operator*(double alpha) const;
			Matrix Matrix::operator*=(double alpha);
			Matrix Matrix::operator*(const Matrix &m) const;		// simple matrix multiplication

			//
			// Especial matrix multiplications
			//
			void Matrix::mmCaABbC(double alpha, double beta, const Matrix &A, const Matrix &B);
			void Matrix::bmmCaABbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
			void Matrix::mmCasABbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
			void Matrix::mmCaAsBbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
			void Matrix::mmCaAUTBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
			void Matrix::mmCaAATbC(double alpha, double beta, const Matrix&A);
			void Matrix::mmCaATAbC(double alpha, double beta, const Matrix&A);
			void Matrix::mmCaATBbC(double alpha, double beta, const Matrix&A, const Matrix&B);
			void Matrix::mmCaATBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
			void Matrix::mmCaABTbC(double alpha, double beta, const Matrix&A, const Matrix&B);
			void Matrix::mmCaABTbC(int r, bool up, double alpha, double beta, const Matrix&A, const Matrix&B);
			void Matrix::bmmCaABTbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
			void Matrix::mmCaIBbC(double alpha, double beta, const Matrix&B);
			void Matrix::mmCaIBbC(double alpha, double beta, int *piv, bool rows, const Matrix&B);
			void Matrix::mmCaAIbC(double alpha, double beta, const Matrix&A);
			void Matrix::mmCaAIbC(double alpha, double beta, const Matrix &A, int *piv, bool rows);

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
			void Matrix::cpermutem(int *piv, bool trans = false);
			void Matrix::rpermutem(int *piv);
			void Matrix::transpose();
			Matrix Matrix::asTranspose() const;
			void Matrix::shift();
			bool Matrix::isId() const;
			// bool Matrix::isId(double eps);
			// bool Matrix::isId(int m1, int m2, int n1, int n2, double eps);
			bool Matrix::isZero() const;
			// bool Matrix::isZero(double eps);
			void Matrix::set2Id();
			// int Matrix::set2Id(int m, int n);
			// int Matrix::set2Id(int m1, int m2, int n1, int n2);
			void Matrix::set2zero();
			// int Matrix::set2zero(int m, int n);
			// int Matrix::set2zero(int m1, int m2, int n1, int n2);
			void Matrix::set2val(double v);
			// int Matrix::set2val(int i, int j, double val);
			// int Matrix::trinvm(Matrix &Am);
			// int Matrix::trinvm(int r, Matrix &Am);
			// bool Matrix::mcompare(Matrix &B);
			// bool Matrix::mcompare(Matrix &B, double eps);
			// bool Matrix::mcompare(Matrix &B, int r, double eps);

			//
			// Print out functions
			//
			void Matrix::print();
			void Matrix::print(const char *name);
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
};

#endif /* _MATRIX_H_ */
