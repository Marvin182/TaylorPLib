/*!
 * \file matrixseq.h
 * \brief Header file for \file matrixseq.cpp .
 * \author D. Monett Diaz
 * \date September, 2006
 * 
 * This is the header of the class that implements the functions needed to 
 * the matrix sequence computation.
 * 
 * Included files:
 *  \file stdio.h for printing out.
 *  \file esccolors.h for colored output.
 * 
 */
 
#ifndef MATRIXSEQ_H_
#define MATRIXSEQ_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include "adolc/adolc.h"			// for ADOL-C functionalities
#include "dae.h"					// to define the DAE, i.e. f
#include "dyn.h"					// to define the dynamic, i.e. d
#include "tra.h"					// to define the trajectory, i.e. x_*
#include "Matrix.h"					// for algebra with matrices
#include "linalg.h"					// for lineal algebra, e.g., operations on matrices.
#include "TPolyn.h"					// for Taylor arithmetic
//#include "time.h"					// for time-related functions and variables
#include "defs.h"					// for global definitions
#include "defaults.h"				// for global definitions
#include "IDOptions.h"				// for global options
#include "IDExample.h"				// for using user example classes

//
// MEMBER FUNCTIONS
//

//
// Main function calls
//
int daeindex( int argc, char *argv[], IDExample * pobj );
int daeindex( int argc, char *argv[], IDExample * pobj, IDOptions op );
//
// Global initializations, allocations/deallocatons
//
int initialize( IDExample * pobj );
int cleanup();
//
// Definitions and calculations involvind ADOL-C active sections
//
int asectra( IDExample *pobj );
int asecdyn( IDExample *pobj );
int asecdae( IDExample *pobj );
//
// Other functions related to the matrix sequence, local ones
//
void shift( int m, int n, double **M1, double **M2 );
int ABD();
int matrixSqc( int& r );
//
// Test functions
//
void testsMatrixClass();

#endif /*MATRIXSEQ_H_*/
