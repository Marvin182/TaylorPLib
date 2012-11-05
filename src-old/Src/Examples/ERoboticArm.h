/*!
 * \file ERoboticArm.h
 * \brief Example from Lamour. 
 * \author D. Monett Díaz
 * \date September, 2006
 * 
 * This class inherits from the abstract base class IDExample.
 * 
 */
 
#ifndef EROBOTICARM_H_
#define EROBOTICARM_H_

#include "../IDExample.h"					/**< Abstract base class. */

class ERoboticArm : public IDExample	// the class ERoboticArm inherits from the class IDExample
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
		ERoboticArm( void ) : IDExample( 6, 8, 0.0, "ERoboticArm" )
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
			ptra[ 4 ] = sin( t ) + cos( t );
			ptra[ 5 ] = log( 3 + t );
			ptra[ 6 ] = sin( t ) - cos( t );
			ptra[ 7 ] = 4*cos( t ) + 1;
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
			pd[ 4 ] = px[ 4 ];
			pd[ 5 ] = px[ 5 ];
		}

		/**
		 * Evaluates a function on a given parameter. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on that argument, as an \type adouble.
		 * 
 		 */
		adouble p1( adouble x )
		{
			return cos( exp( x ) - 1 ) + cos( x - 1 );
		}

		/**
		 * Evaluates a function on a given parameter. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on that argument, as an \type adouble.
		 * 
 		 */
		adouble p2( adouble x )
		{
			return sin( 1 - exp( x )) + sin( 1 - x );
		}

		/**
		 * Evaluates a function on a given parameter. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on that argument, as an \type adouble.
		 * 
 		 */
		adouble a( adouble x )
		{
			return 2 / ( 2 - cos( x )*cos( x ));
		}

		/**
		 * Evaluates a function on a given parameter. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on that argument, as an \type adouble.
		 * 
 		 */
		adouble b( adouble x )
		{
			return cos( x ) / ( 2 - cos( x )*cos( x ));
		}

		/**
		 * Evaluates a function on a given parameter. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on that argument, as an \type adouble.
		 * 
 		 */
		adouble c( adouble x )
		{
			return sin( x ) / ( 2 - cos( x )*cos( x ));
		}

		/**
		 * Evaluates a function on given parameters. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on that argument, as an \type adouble.
		 * 
 		 */
		adouble d( adouble x )
		{
			return cos( x )*sin( x ) / ( 2 - cos( x )*cos( x ));
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
			pf[ 0 ] = py[ 0 ] - px[ 3 ];
			pf[ 1 ] = py[ 1 ] - px[ 4 ];
			pf[ 2 ] = py[ 2 ] - px[ 5 ];
			pf[ 3 ] = py[ 3 ] - 2*c( px[ 2 ] )*(px[ 3 ] 
					+ px[ 5 ])*(px[ 3 ] + px[ 5 ])
					- px[ 3 ]*px[ 3 ]*d( px[ 2 ] ) 
					- ( 2*px[ 2 ] - px[ 1 ])*( a( px[ 2 ] )
					+ 2*b( px[ 2 ] )) - a( px[ 2 ] )*px[ 6 ] 
					+ a( px[ 2 ] )*px[ 7 ];
			pf[ 4 ] = py[ 4 ] + 2*c( px[ 2 ] )*(px[ 3 ] 
					+ px[ 5 ])*(px[ 3 ] + px[ 5 ])
					+ px[ 3 ]*px[ 3 ]*d( px[ 2 ] ) 
					- ( 2*px[ 2 ] - px[ 1 ])*( 1 - 3*a( px[ 2 ] )
					- 2*b( px[ 2 ] )) + a( px[ 2 ] )*px[ 6 ] 
					- ( a( px[ 2 ] ) + 1 )*px[ 7 ];
			pf[ 5 ] = py[ 5 ] + 2*c( px[ 2 ] )*(px[ 3 ] 
					+ px[ 5 ])*(px[ 3 ] + px[ 5 ])
					+ px[ 3 ]*px[ 3 ]*d( px[ 2 ] ) 
					- ( 2*px[ 2 ] - px[ 1 ])*( a( px[ 2 ] )
					- 9*b( px[ 2 ] )) + 2*px[ 3 ]*px[ 3 ]*c( px[ 2 ] ) 
					+ (px[ 3 ] + px[ 5 ])*(px[ 3 ] + px[ 5 ])*d( px[ 2 ] )
					+ ( a( px[ 2 ] ) + b( px[ 2 ] ))*px[ 6 ]
					- ( a( px[ 2 ] ) + b( px[ 2 ] ))*px[ 7 ];
			pf[ 6 ] = cos( px[ 0 ] ) + cos( px[ 0 ] + px[ 2 ] ) - p1( t );
			pf[ 7 ] = sin( px[ 0 ] ) + sin( px[ 0 ] + px[ 2 ] ) - p2( t );
		}
};

#endif /*EROBOTICARM_H_*/
