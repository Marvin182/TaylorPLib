/**
FILE NAME:  tra.cpp
WRITTEN BY: D. Monett Díaz
FOR:        DAE-AD HU-project
DUE DATE:   April, 2006
PURPOSE: 
	To define x_* for the Tischendorf-Lamour's example.
MEMBER FUNCTIONS:
	void tra( double t, adouble *ptra )
		This function will define the x_* variable, i.e. the trajectory, as a vector of adoubles.
*/

#include "tra.h"

/**
 * Definition of the trajectory.
 * 		x_*(t): R in Rn
 * 
 * \param[in] t The independent time variable.
 * \param[out] ptra An \type adouble vector defining the trajectory. Its dymension is \a n.
 * 
 */
int tra( adouble t, adouble *ptra )
{
	ptra[ 0 ] = sin( t );
	ptra[ 1 ] = cos( t );
	ptra[ 2 ] = log( 1 + t ) + 2;
	
	return 3;
}
