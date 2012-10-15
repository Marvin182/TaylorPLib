/**
FILE NAME:  tra.h
WRITTEN BY: D. Monett D�z
FOR:        DAE-AD HU-project
DUE DATE:   April, 2006
PURPOSE: 
	Header file for tra.cpp.
MEMBER FUNCTIONS:
	void tra( double t, adouble *ptra )
		This function will define the x* variable, i.e. the trajectory, as a vector of adoubles.
MEMBER DATA:
	none
**/

#ifndef _TRA_H_
#define _TRA_H_

#include "adolc/adolc.h"	// for ADOL-C functionalities
#include <math.h>

//
// MEMBER FUNCTIONS
//
//void tra( adouble t, adouble *ptra );
int tra( adouble t, adouble *ptra );

#endif /*_TRA_H_*/
