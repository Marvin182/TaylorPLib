/**
FILE NAME:  Defaults.h
WRITTEN BY: D. Monett D�z
FOR:        DAE-AD HU-project
DUE DATE:   September 14, 2006
PURPOSE: 
	Header containing global options' default values.
**/

#ifndef DEFAULTS_H_
#define DEFAULTS_H_

//
// Name identifiers of global options
//
const int		_TP_DEGREE_ = 1;					// Taylor polynomial's degree: name identifier
const int		_QR_EPS_ = 2;						// QR threshold: name identifier
const int		_PRINT_PARAM_WHAT_ = 3;				// what to print out: name identifier
const int		_PRINT_PARAM_WHERE_ = 4;			// where to print out: name identifier
const int		_IO_EPS_ = 5;						// I/O threshold name identifier
//
// Values by default of global options
//
const int 		_TP_DEGREE_VAL = 10;				// degree of Taylor coefficients; highest derivative order;
													// Taylor polynomial's degree
const double 	_QR_EPS_VAL = 1E-15;				// threshold for the Householder QR fact. with column pivoting
const int		_PRINT_PARAM_WHAT_VAL = 0;			// print out parameter for controlling output
													// (e.g., =0 print out nothing)
const int		_PRINT_PARAM_WHERE_VAL = 0;			// output to the screen (=0), to a datafile (=1),
													// or nowhere (=2)
const double 	_IO_EPS_VAL = 1E-15;				// I/O threshold

#endif /*DEFAULTS_H_*/
