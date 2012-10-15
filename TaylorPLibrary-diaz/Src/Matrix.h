/*!
 * \file Matrix.h
 * \brief Header file for the class \file Matrix.cpp.
 * \author D. Monett D�z
 * \date November, 2006
 * 
 * This is the header of the class for matrix arithmetic.
 * 
 * Included files:
 * 	\file IDException.h and <stdexcept> for handling exceptions.
 * 
 */
#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdexcept>
#include <stdio.h>
#include <math.h>
#include <typeinfo>
#include "IDException.h"
#include "esccolors.h"
#include "ioutil.h"
#include "adolc/adolc.h"									// ADOL-C functionalities
#include "TPolyn.h"											// Taylor arithmetic

template <typename T> class Matrix
{
	private:
		int _rows,											// The number of rows
			_cols,											// The number of columns
			_dimT;											// The dimension of the type T
    	T **_data;											// The pointer to the allocated memory
		int _maxr,											// The max. number of rows
			_maxc;											// The max. number of columns

	public:
		//
		// Constructors, destructor
		//
		Matrix::Matrix();									// Default constructor
		Matrix::Matrix( int r, int c );						// Regular constructor
		Matrix::Matrix( int r, int c, int dimT );			// Regular constructor
		Matrix::Matrix( const Matrix &m );					// Copy constructor
		Matrix Matrix::redim( int r, int c );				// Reset the current object dimensions
		Matrix Matrix::redim( int r, int c, int dimT );		// Reset the current object dimensions
		Matrix::~Matrix();									// Destructor
		//
		// Accessing properties
		//
		T** data() { return _data; }						// Returns the ptr to the alloc. memory
		int nrows() { return _rows; }						// Returns the number of rows
		int ncols() { return _cols; }						// Returns the number of columns
		int dimT() { return _dimT; }						// Returns the dimension of the type T
		int maxr() { return _maxr; }						// Returns the max. number of rows
		int maxc() { return _maxc; }						// Returns the max. number of columns
		void setnrows( int v ) { _rows = v; }				// Sets the number of rows
		void setncols( int v ) { _cols = v; }				// Sets the number of columns
		void setdimT( int v ) { _dimT = v; }				// Sets the dimension of T
		void setdim( int v1, int v2 ) {
			_rows = v1; _cols = v2; }						// Sets both the nr. rows and columns
		void setdim( int v1, int v2, int v3 ) {
			_rows = v1; _cols = v2; _dimT = v3; }			// Sets the dimensions
		void setmaxdim( int v1, int v2 ) {
			_maxr = v1; _maxc = v2; }						// Sets both max. dimensions
		//
		// Overloaded operators
		//
		T & Matrix::operator()( int i, int j );				// Element
		Matrix Matrix::operator=( const Matrix &m );		// Assignment operators
		bool Matrix::operator==( const Matrix &m );			// Comparison operators
		bool Matrix::operator!=( const Matrix &m );
		Matrix Matrix::operator+( const Matrix &m );		// Arithmetic operators
		Matrix Matrix::operator+=( const Matrix &m );
		Matrix Matrix::operator-( const Matrix &m );
		Matrix Matrix::operator-=( const Matrix &m );
		Matrix Matrix::operator-();
		Matrix Matrix::operator*( double alpha );
		//
		// Especial matrix multiplications
		//
		int Matrix::mmCaABbC( double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::bmmCaABbC( int r, int c, double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::mmCasABbC( int r, double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::mmCaAsBbC( int r, double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::mmCaAUTBPbC( double alpha, double beta, Matrix &A, Matrix &B, int *piv );
		int Matrix::mmCaAATbC( double alpha, double beta, Matrix &A );
		int Matrix::mmCaATAbC( double alpha, double beta, Matrix &A );
		int Matrix::mmCaATBbC( double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::mmCaATBPbC( double alpha, double beta, Matrix &A, Matrix &B, int *piv );
		int Matrix::mmCaABTbC( double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::mmCaABTbC( int r, bool up, double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::bmmCaABTbC( int r, int c, double alpha, double beta, Matrix &A, Matrix &B );
		int Matrix::mmCaIBbC( double alpha, double beta, Matrix &B );
		int Matrix::mmCaIBbC( double alpha, double beta, int *piv, bool rows, Matrix &B );
		int Matrix::mmCaAIbC( double alpha, double beta, Matrix &A );
		int Matrix::mmCaAIbC( double alpha, double beta, Matrix &A, int *piv, bool rows );
		//
		// Solving systems of equations
		//
		int Matrix::utsolve( Matrix &B );					// U X = B, back-substitution
		int Matrix::utsolve( Matrix &B, Matrix &X, int *piv );// U X = B, back-substitution, pivoting
		int Matrix::utsolve( T *b );						// U x = b, back-substitution
		int Matrix::utxsolve( Matrix &B );					// X U = B, back-substitution
		int Matrix::gsolve( Matrix &B );					// A X = B, Gaussian elimination
		int Matrix::gsolve( T *b );							// A x = b, Gaussian elimination
		//
		// Other functions
		//
		double Matrix::fnorm();								// Frobenius norm
		int Matrix::tcfnorm( int nrTC, double *fn );		// vector of Frobenius norms
		//int Matrix::colnorm( T *c );
		int Matrix::colnorm( double *c );
		//int Matrix::colnormdown( int pos, int *piv, T *c );
		int Matrix::colnormdown( int pos, int *piv, double *c );
		int Matrix::cpermutem( int *piv );
		int Matrix::cpermutem( int *piv, bool trans );
		int Matrix::rpermutem( int *piv );
		Matrix Matrix::transpm();
		int Matrix::shift( Matrix &M );
		bool Matrix::isId();
		bool Matrix::isId( double eps );
		bool Matrix::isId( int m1, int m2, int n1, int n2, double eps );
		bool Matrix::isZero();
		bool Matrix::isZero( double eps );
		int Matrix::set2Id();
		int Matrix::set2Id( int m, int n );
		int Matrix::set2Id( int m1, int m2, int n1, int n2 );
		int Matrix::set2zero();
		int Matrix::set2zero( int m, int n );
		int Matrix::set2zero( int m1, int m2, int n1, int n2 );
		int Matrix::set2val( double v );
		int Matrix::set2val( int i, int j, double val );
		int Matrix::trinvm( Matrix &Am );
		int Matrix::trinvm( int r, Matrix &Am );
		bool Matrix::mcompare( Matrix &B );
		bool Matrix::mcompare( Matrix &B, double eps );
		bool Matrix::mcompare( Matrix &B, int r, double eps );
		//
		// Print out functions
		//
		void Matrix::printm( char *str );
		void Matrix::fprintm( FILE * fn, char *str );
		void Matrix::printm( char *str, const char *const color );
		void Matrix::printdm( char *str );
		void Matrix::fprintdm( FILE * fn, char *str );
		void Matrix::printdm( char *str, const char *const color );
		void Matrix::printdm( char *str, const char *const color, double eps );
		void Matrix::fprintdm( FILE * fn, char *str, double eps );
		void Matrix::printm( int *piv, char *str, const char *const color );
		void Matrix::printdm( int *piv, char *str, const char *const color, double eps );
		void Matrix::fprintm( FILE * fn, int *piv, char *str );
		void Matrix::fprintdm( FILE * fn, int *piv, char *str, double eps );
		void Matrix::printtpm( char *str, const char *const color, double eps );
		void Matrix::fprinttpm( FILE * fn, char *str, double eps );
		void Matrix::printtpm( int m1, int m2, int n1, int n2,
				char *str, const char *const color, double eps );
		void Matrix::fprinttpm( FILE * fn, int m1, int m2, int n1, int n2,
				char *str, double eps );
		void Matrix::printtpm( int *piv, char *str, const char *const color, double eps );
		void Matrix::fprinttpm( FILE * fn, int *piv, char *str, double eps );

};

#endif /*MATRIX_H_*/
