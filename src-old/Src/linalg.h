/*!
 * \file linalg.h
 * \brief Header file for \file linalg.cpp .
 * \author D. Monett Diaz
 * \date January, 2007
 * 
 * This is the header of the class that implements some linear algebra algorithms 
 * and related functions.
 * 
 * Included files:
 * 
 */
 
#ifndef LINALG_H_
#define LINALG_H_

#include "Matrix.h"											// Matrix class
#include "adolc/adolc.h"									// for ADOL-C functionalities
#include "TPolyn.h"											// for Taylor arithmetic

//
// MEMBER FUNCTIONS
//

//
// Miscellanea
//
int max( double *v, bool *selected, int *P, int a, int b );
int findIndex( int *v, int t );
template <typename T> T dpc( int r, int m, T *x, T *y );
int transpv( int n, int *x, int *y );
bool vleqeps( int n, double *v, double eps, int & elem );

//
// Householder QR
//
template <typename T> int constructR( Matrix<T> &A, int *piv, Matrix<T> &R, bool permuted );
template <typename T> int hhlh( int pos, T beta, Matrix<T> &A, int *piv, T *v );
template <typename T> int hhlr( int pos, double T, Matrix<T> &A, int *piv, T *v );
template <typename T> int hhr( int rowpos, int colpos, double T, Matrix<T> &A, T *v );
template <typename T> int hm( int r, double eps, Matrix<T> &A, int *P, T *v, T beta );
template <typename T> int hhv( int r, Matrix<T> &A, int *P, T *v, double eps, T &beta );
template <typename T> int hqrfcp( double eps, Matrix<T> &A, int *P, Matrix<T> &Q, int& rank, 
		bool permuted, int printWhat, int printWhere, FILE * fn );
template <typename T> int hqrfcp( double eps, Matrix<T> &A, int *P, Matrix<T> &Q, 
		Matrix<T> &R, int& rank, bool permuted, int printWhat, int printWhere, FILE * fn );
//
// Other operations with/over matrices
//
template <typename T> int pwq( Matrix<T> &G, Matrix<T> &Gm, Matrix<T> &P, 
		Matrix<T> &W, Matrix<T> &Q );
template <typename T> int decompS( int r, Matrix<T> &R, int *P, Matrix<T> &R1, Matrix<T> &R2 );
template <typename T> int decompV( int r, Matrix<T> &V, Matrix<T> &V1, Matrix<T> &V2 );
template <typename T> int decompW( int r, Matrix<T> &W, Matrix<T> &W1, Matrix<T> &W2 );
template <typename T> int formBM( int d, Matrix<T> &BM11, Matrix<T> &BM12, 
		Matrix<T> &BM21, Matrix<T> &BM22, Matrix<T> &BM );
template <typename T> int formBM( int d, bool vert, Matrix<T> &BM1, Matrix<T> &BM2, 
		Matrix<T> &BM );
template <typename T> int filloutS( int r, Matrix<T> &S1, Matrix<T> &S );
template <typename T> int filloutM( int d1, int d2, Matrix<T> &V, Matrix<T> &M );
template <typename T> int computeTlim( Matrix<T> &mi1, Matrix<T> &Si, 
		Matrix<T> &Vitilde, Matrix<T> &Tlim );
template <typename T> int computemTui( Matrix<T> &mi2, Matrix<T> &Si, 
		Matrix<T> &Uitilde, Matrix<T> &mTui );
template <typename T> int computeV( Matrix<T> &M, Matrix<T> &Vi );
template <typename T> int computeV( int *P, bool col, Matrix<T> &M, Matrix<T> &Vi );
template <typename T> int computeV( int *P, Matrix<T> &M, Matrix<T> &Vi );
template <typename T> int computeGm( Matrix<T> &m1, Matrix<T> &m2, 
		Matrix<T> &V, Matrix<T> &R1, Matrix<T> &R1m, Matrix<T> &Um, Matrix<T> &Gm );
template <typename T> int computemi1( Matrix<T> &y, Matrix<T> &U, Matrix<T> &m, 
		Matrix<T> &mi1 );
template <typename T> int computemi1( Matrix<T> &y, Matrix<T> &U, Matrix<T> &mi1 );
template <typename T> int computeDm( int r, Matrix<T> &SA, Matrix<T> &SD, 
		int *PA, int *PD, Matrix<T> &UD, Matrix<T> &Dm );

#endif /*LINALG_H_*/
