/**
FILE NAME:  defs.h
WRITTEN BY: D. Monett D�z
FOR:        DAE-AD HU-project
DUE DATE:   February 14, 2006
PURPOSE: 
	Header containing global constant definitions.
**/

#ifndef _DEFS_H_
#define _DEFS_H_

#include <limits.h>
#include "adolc/adolc.h"					// for ADOL-C functionalities
#include "TPolyn.h"							// Taylor Polynomials class

//
// Global type definitions
//
typedef Matrix<double> DoubleMatrix;
typedef Matrix<int> IntMatrix;
typedef Matrix<short> ShortMatrix;
typedef Matrix<adouble> aDoubleMatrix;
typedef Matrix<TPolyn> TPolynMatrix;

//
// Global constant definitions
//
const double	_NR_COEFF = 50;
const double	_NR_EQ    = 50;
const int		_MAX_NR   = INT_MAX;						// 32767

#endif /*_DEFS_H_*/
