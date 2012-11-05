/**
FILE NAME:  dyn.h
WRITTEN BY: D. Monett D�z
FOR:        DAE-AD HU-project
DUE DATE:   April, 2006
PURPOSE: 
	Header file for dyn.cpp.
MEMBER FUNCTIONS:
	void dyn( adouble *px, adouble t, adouble *pd )
		This function will define the d variable, i.e. the dynamic, as a vector of adoubles.
MEMBER DATA:
	none
**/

#ifndef DYN_H_
#define DYN_H_

#include "adolc/adolc.h"	// for ADOL-C functionalities
#include <math.h>

//
// MEMBER FUNCTIONS
//
//void dyn( adouble *px, adouble t, adouble *pd );
int dyn( adouble *px, adouble t, adouble *pd );

#endif /*DYN_H_*/
