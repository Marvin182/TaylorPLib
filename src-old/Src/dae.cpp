/**
FILE NAME:  dae.cpp
WRITTEN BY: D. Monett Díaz
FOR:        DAE-AD HU-project
DUE DATE:   April, 2006
PURPOSE: 
	To define f for the Tischendorf-Lamour's example.
MEMBER FUNCTIONS:
	void dae( adouble *py, adouble *px, adouble t, adouble *pf )
		This function will define the f variable, i.e. the DAE, as a vector of adoubles.  
**/

#include "dae.h"

/**
 * Auxiliary function.
 *
 * \param[in] x The independent \type adouble variable of the function to evaluate.
 * \return The value of the function on that argument, as an \type adouble.
 */
adouble g( adouble x )
{
	return exp( x ) - 1;
}

/**
 * Definition of the DAE, the system of differential-algebraic equations.
 * Also system of differential equations with additional algebraic constraints
 * 		f(y, x, t): Rm+n+1 in Rn
 * 
 * \param[in] px The pointer to \a x, an \type adouble vector of algebraic variables.
 * \param[in] py The pointer to \a y, an \type adouble vector of differential variables.
 * \param[in] t The independent time variable.
 * \param[out] pf An \type adouble vector defining the DAE. Its dymension is \a n.
 *
 */
int dae( adouble *py, adouble *px, adouble t, adouble *pf )
{
	adouble temp = g(px[ 0 ] - px[ 1 ]);
	pf[0] = py[ 0 ] + px[ 2 ];
	pf[1] = py[ 1 ] - px[ 2 ];
	pf[2] = px[ 2 ] - temp;
	
	return 3;
}
