/**
FILE NAME:  dae.h
WRITTEN BY: D. Monett D�z
FOR:        DAE-AD HU-project
DUE DATE:   April, 2006
PURPOSE: 
	Header file for dae.cpp.
MEMBER FUNCTIONS:
	void dae( adouble *py, adouble *px, adouble t, adouble *pf )
		This function will define the f variable, i.e. the DAE, as a vector of doubles. 
MEMBER DATA:
	none
**/

#ifndef DAE_H_
#define DAE_H_

#include "adolc/adolc.h"	// for ADOL-C functionalities
#include <math.h>

//
// MEMBER FUNCTIONS
//
//void dae( adouble *py, adouble *px, adouble t, adouble *pf );
int dae( adouble *py, adouble *px, adouble t, adouble *pf );

#endif /*DAE_H_*/
