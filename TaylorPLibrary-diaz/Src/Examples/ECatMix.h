/*!
 * \file ECatMix.h
 * \brief Example from Lamour. 
 * \author D. Monett Díaz
 * \date September, 2006
 * 
 * This class inherits from the abstract base class IDExample.
 * 
 */
 
#ifndef ECATMIX_H_
#define ECATMIX_H_

#include "../IDExample.h"					/**< Abstract base class. */

class ECatMix : public IDExample		// the class ECatMix inherits from the class IDExample
{
	public:
		/**
		 * Constructor particular to this example. 
		 * 
		 * No parameters are needed, only those of the constructor of the abstract base class:
		 * 
		 * 1st par. : m_dyn = nr. dep. var. from dynamic 
		 * 					= nr. equations in dyn(...)
		 * 2nd par. : m_dae = nr. dep. var. from DAE 
		 * 					= nr. equations in dae(...) 
		 * 					= nr. equations in tra(...)
		 * 3rd par. : indep. var time  = initial time value
		 * 4th par. : class name
		 * 
		 * \param[in] m_dyn The number of dependent variables of the dynamic, as an \type int.
		 * 		It must be equal to the number of equations that define the dynamic.
		 * \param[in] m_dae The number of dependent variables of the DAE, as an \type int.
		 * 		It must be equal to the number of equations that define the DAE, as well as 
		 * 		be equal to the number of equations that define the trajectory.
		 * \param[in] t The initial value for the independent variable time, as a \type double.
		 * \param[in] cn The class name.
		 * 
		 */
		ECatMix( void ) : IDExample( 4, 7, 0.0, "ECatMix" )
		{
		}
		
		/**
		 * Defines the trajectory for this example.
		 * 
		 * \param[in] t The independent variable time, as an \type adouble.
		 * \param[out] ptra A pointer to an \type adouble vector defining the trajectory.
		 * 
		 */
		void tra( adouble t, adouble *ptra )
		{
			ptra[ 0 ] = sin( t );
			ptra[ 1 ] = cos( t );
			ptra[ 2 ] = log( 1 + t );
			ptra[ 3 ] = 2*sin( t );
			ptra[ 4 ] = 3*cos( t );
			ptra[ 5 ] = 2*t;
			ptra[ 6 ] = sin( 2*t ) + 2;
		}
		
		/**
		 * Dynamic for this example.
		 * 
		 * \param[in] px The pointer to an \type adouble vector of algebraic variables.
		 * \param[in] t The independent variable time, as an \type adouble.
		 * \param[out] pd A pointer to an \type adouble vector defining the dynamic.
		 * 
		 */
		void dyn( adouble *px, adouble t, adouble *pd )
		{
			pd[ 0 ] = px[ 0 ];
			pd[ 1 ] = px[ 1 ];
			pd[ 2 ] = px[ 2 ];
			pd[ 3 ] = px[ 3 ];
		}

		/**
		 * Evaluates a function on given parameters. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The first \type adouble variable of the function to evaluate.
		 * \param[in] y The second \type adouble variable to compare with.
		 * \return The maximum value, as an \type adouble.
		 * 
 		 */
		adouble pminus( adouble x, adouble y )
		{			
			return fmax( -x, y );		// function 'fmax' from ADOL-C
		}

		/**
		 * Evaluates a function on given parameters. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The first \type adouble variable of the function to evaluate.
		 * \param[in] y The second \type adouble variable to compare with.
		 * \return The maximum value, as an \type adouble.
		 * 
 		 */
		adouble pplus( adouble x, adouble y )
		{			
			return fmax( x, y );		// function 'fmax' from ADOL-C
		}

		/**
		 * DAE for this example.
		 *
		 * \param[in] py The pointer to the \type adouble vector of differential variables.
		 * \param[in] px The pointer to the \type adouble vector of algebraic variables.
		 * \param[in] t The independent variable time, as an \type adouble.
		 * \param[out] pf A pointer to an \type adouble vector defining the DAE.
		 * 
		 */
		void dae( adouble *py, adouble *px, adouble t, adouble *pf )
		{
			adouble tf = 1; //0.625; //0.4;
			adouble _0d0 = 0.0;
			adouble _1d0 = 1.0;
			
			pf[ 0 ] = py[ 0 ] / tf - ( 10*px[ 1 ] - px[ 0 ] )*px[ 4 ];
			pf[ 1 ] = py[ 1 ] / tf + px[ 1 ] - ( px[ 0 ] - 9*px[ 1 ] )*px[ 4 ];
			pf[ 2 ] = py[ 2 ] / tf - ( px[ 2 ] - px[ 3 ] )*px[ 4 ];
			pf[ 3 ] = py[ 3 ] / tf - ( px[ 3 ] + _1d0 ) 
					- ( 9*px[ 3 ] - 10*px[ 2 ] - _1d0 )*px[ 4 ];
			pf[ 4 ] = ( px[ 2 ] - px[ 3 ] )*( px[ 0 ] - 10*px[ 1 ] )
					- ( 1 + px[ 3 ] )*px[ 1 ] + pminus( px[ 5 ], _0d0 ) 
					- pminus( px[ 6 ], _0d0 );
			pf[ 5 ] = pplus( _0d0, px[ 5 ] ) - px[ 4 ];
			pf[ 6 ] = pplus( _0d0, px[ 6 ] ) - _1d0 + px[ 4 ];
		}
};

#endif /*ECATMIX_H_*/
