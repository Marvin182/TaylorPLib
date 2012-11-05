/**
FILE NAME:  dyn.cpp
WRITTEN BY: D. Monett Díaz
FOR:        DAE-AD HU-project
DUE DATE:   April, 2006
PURPOSE: 
	To define d for the Tischendorf-Lamour's example.
MEMBER FUNCTIONS:
	void dyn( adouble *px, adouble t, adouble *pd )
		This function will define the d variable, i.e. the dynamic, as a vector of adoubles.
*/

#include "dyn.h"

/**
 * Definition of the dynamic.
 * 		d(x, t): Rm+1 in Rm
 * 
 * \param[in] px The pointer to \a x, an \type adouble vector of algebraic variables.
 * \param[in] t The independent time variable.
 * \param[out] pd An \type adouble vector defining the dynamic. Its dymension is \a m.
 * 
 */
int dyn( adouble *px, adouble t, adouble *pd )
{
	pd[0] = -2 * sqrt( 1 - px[ 0 ] );
	pd[1] = sin( t ) * px[ 1 ];//3*px[ 2 ] + 2*px[ 1 ] + 5 * t;
	
	return 2;
}
