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

    		void allocateMemory(bool initialize)
    			{ Matrix::allocateMemory(_data, _rows, _cols, initialize); }
    		static void allocateMemory(double **&data, int rows, int cols, bool initialize);
	    	void deallocateMemory()
	    		{ Matrix::deallocateMemory(_data, _rows, _cols); }
	    	static void deallocateMemory(double **&data, int rows, int cols);
			void copyFrom(const Matrix &m);

	public:
			//
			// Constructors, destructor
			//
			Matrix();														// Default constructor
			Matrix(int rows, int cols, bool initialize = true);				// Regular constructor
			Matrix(int rows, int cols, int dimT, bool initialize = true);	// Regular constructor
			Matrix(const Matrix &m);										// Copy constructor

			Matrix(int rows, int cols, double *values);						// Test construtor		

			~Matrix();														// Destructor

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
			double & operator()(int row, int col);			// Element
			Matrix operator=(const Matrix &m);				// Assignment operators
			bool operator==(const Matrix &m) const;			// Comparison operators
			bool operator!=(const Matrix &m) const;
			Matrix operator+(const Matrix &m) const;		// Arithmetic operators
			Matrix operator+=(const Matrix &m);
			Matrix operator-(const Matrix &m) const;
			Matrix operator-=(const Matrix &m);
			Matrix operator-();
			Matrix operator*(double alpha) const;
			Matrix operator*=(double alpha);
			Matrix operator*(const Matrix &m) const;		// simple matrix multiplication

			//
			// Especial matrix multiplications
			//
			void mmCaABbC(double alpha, double beta, const Matrix &A, const Matrix &B);
			void bmmCaABbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
			void mmCasABbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
			void mmCaAsBbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
			void mmCaAUTBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
			void mmCaAATbC(double alpha, double beta, const Matrix&A);
			void mmCaATAbC(double alpha, double beta, const Matrix&A);
			void mmCaATBbC(double alpha, double beta, const Matrix&A, const Matrix&B);
			void mmCaATBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
			void mmCaABTbC(double alpha, double beta, const Matrix&A, const Matrix&B);
			void mmCaABTbC(int r, bool up, double alpha, double beta, const Matrix&A, const Matrix&B);
			void bmmCaABTbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
			void mmCaIBbC(double alpha, double beta, const Matrix&B);
			void mmCaIBbC(double alpha, double beta, int *piv, bool rows, const Matrix&B);
			void mmCaAIbC(double alpha, double beta, const Matrix&A);
			void mmCaAIbC(double alpha, double beta, const Matrix &A, int *piv, bool rows);

			//
			// Solving systems of equations
			//
			// int utsolve(Matrix &B);					// U X = B, back-substitution
			// int utsolve(Matrix &B, Matrix &X, int *piv);// U X = B, back-substitution, pivoting
			// int utsolve(T *b);						// U x = b, back-substitution
			// int utxsolve(Matrix &B);					// X U = B, back-substitution
			// int gsolve(Matrix &B);					// A X = B, Gaussian elimination
			// int gsolve(T *b);							// A x = b, Gaussian elimination

			//
			// Other functions
			//
			// double fnorm();								// Frobenius norm
			// int tcfnorm(int nrTC, double *fn);		// vector of Frobenius norms
			// int colnorm(double *c);
			// int colnormdown(int pos, int *piv, double *c);
			void cpermutem(int *piv, bool trans = false);
			void rpermutem(int *piv);
			void transpose();
			Matrix asTranspose() const;
			void shift();
			bool isId() const;
			// bool isId(double eps);
			// bool isId(int m1, int m2, int n1, int n2, double eps);
			bool isZero() const;
			// bool isZero(double eps);
			void set2Id();
			void set2Id(int top, int bottom, int left, int right);
			void set2IdFromIndices(int firstRow, int lastRow, int firstCol, int lastCol);
			void set2Zero();
			void set2Zero(int top, int bottom, int left, int right);
			void set2ZeroFromIndices(int firstRow, int lastRow, int firstCol, int lastCol);
			void set2Val(double v);
			void set2Val(int top, int bottom, int left, int right, double v);
			void set2ValFromIndices(int firstRow, int lastRow, int firstCol, int lastCol, double v);
			// int trinvm(Matrix &Am);
			// int trinvm(int r, Matrix &Am);
			// bool mcompare(Matrix &B);
			// bool mcompare(Matrix &B, double eps);
			// bool mcompare(Matrix &B, int r, double eps);

			//
			// Print out functions
			//
			void print();
			void print(const char *name);
			// void fprintm(FILE * fn, char *str);
			// void printm(char *str, const char *const color);
			// void printdm(char *str);
			// void fprintdm(FILE * fn, char *str);
			// void printdm(char *str, const char *const color);
			// void printdm(char *str, const char *const color, double eps);
			// void fprintdm(FILE * fn, char *str, double eps);
			// void printm(int *piv, char *str, const char *const color);
			// void printdm(int *piv, char *str, const char *const color, double eps);
			// void fprintm(FILE * fn, int *piv, char *str);
			// void fprintdm(FILE * fn, int *piv, char *str, double eps);
			// void printtpm(char *str, const char *const color, double eps);
			// void fprinttpm(FILE * fn, char *str, double eps);
			// void printtpm(int m1, int m2, int n1, int n2,
			// 		char *str, const char *const color, double eps);
			// void fprinttpm(FILE * fn, int m1, int m2, int n1, int n2,
			// 		char *str, double eps);
			// void printtpm(int *piv, char *str, const char *const color, double eps);
			// void fprinttpm(FILE * fn, int *piv, char *str, double eps);
	};
};

#endif /* _MATRIX_H_ */
