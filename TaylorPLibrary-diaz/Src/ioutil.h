/*!
 * \file ioutil.h
 * \brief Header file for \file ioutil.cpp .
 * \author D. Monett D�z
 * \date June, 2006
 * 
 * This is the header of the class that implements some input/output functions for 
 * vectors and matrices.
 * 
 * Included files:
 *  \file stdio.h for printing out.
 * 	\file esccolors.h for colored output.
 * 
 */
 
#ifndef IOUTIL_H_
#define IOUTIL_H_

#include <stdio.h>
#include "esccolors.h"


//
// MEMBER FUNCTIONS
//

// to the screen, with and without colors...
void printdv( int n, double *v, char *str );
void printdv( int n, double *v, char *str, const char *const color );
void printdv( int n, double *v, int *piv, char *str, const char *const color );
void printdv( int a, int b, double *v, char *str, const char *const color );
void printdv( int a, int b, double *v, int *piv, char *str, const char *const color );
void printiv( int n, int *v, char *str, const char *const color );
void printiv( int n, int *v, int *piv, char *str, const char *const color );
void printm( int m, int n, double **M, char *str );
void printm( int m, int n, double **M, char *str, const char *const color );
void printdm( int m, int n, double **M, char *str );
void printdm( int m, int n, double **M, char *str, const char *const color );
void print3ddm( int m, int n, int p, double ***M, char *str, const char *const color );
void printm( int m, int n, double **M, int *piv, char *str, const char *const color );
// to a file...
void fprintdv( FILE * fn, int n, double *v, char *str );
void fprintdv( FILE * fn, int n, double *v, int *piv, char *str );
void fprintdv( FILE * fn, int a, int b, double *v, char *str );
void fprintdv( FILE * fn, int a, int b, double *v, int *piv, char *str );
void fprintiv( FILE * fn, int n, int *v, char *str );
void fprintiv( FILE * fn, int n, int *v, int *piv, char *str );
void fprintm( FILE * fn, int m, int n, double **M, char *str );
void fprintdm( FILE * fn, int m, int n, double **M, char *str );
void fprint3ddm( FILE * fn, int m, int n, int p, double ***M, char *str );
void fprintm( FILE * fn, int m, int n, double **M, int *piv, char *str );

#endif /*IOUTIL_H_*/
