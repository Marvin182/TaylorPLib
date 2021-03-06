/*!
 * \file EHessenberg.h
 * \brief Example from Lamour. 
 * \author D. Monett D�az
 * \date February, 2007
 * 
 * This class inherits from the abstract base class IDExample.
 * 
 */
 
#ifndef EHESSENBERG_H_
#define EHESSENBERG_H_

#include "../IDExample.h"					/**< Abstract base class. */

class EHessenberg : public IDExample	// the class EHessenberg inherits from the class IDExample
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
		EHessenberg( void ) : IDExample( 3, 4, 1.0, "EHessenberg" )
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
			ptra[ 2 ] = sin( 2 * t );
			ptra[ 3 ] = cos( 2 * t );
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
		}

		/**
		 * Evaluates a function on a given parameter. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x1 The first independent \type adouble variable of the function to evaluate.
		 * \param[in] x2 The second independent \type adouble variable of the function to evaluate.
		 * \param[in] x3 The third independent \type adouble variable of the function to evaluate.
		 * \param[in] t The fourth independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on these arguments, as an \type adouble.
		 * 
 		 */
		adouble alpha( adouble x1, adouble x2, adouble x3, adouble t )
		{
			return t;
		}

		/**
		 * Evaluates a function on a given parameter. It is used by dae(...) and is particular 
		 * to this example.
		 * 
		 * \param[in] x The independent \type adouble variable of the function to evaluate.
		 * \return The value of the function on that argument, as an \type adouble.
		 * 
 		 */
		adouble p( adouble x )
		{
			return 1;
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
			pf[ 0 ] = py[ 0 ] + px[ 0 ] + px[ 3 ];
			pf[ 1 ] = py[ 1 ] + alpha( px[ 0 ], px[ 1 ], px[ 2 ], t ) * px[ 3 ];
			pf[ 2 ] = py[ 2 ] + px[ 0 ] + px[ 1 ] + px[ 2 ];
			pf[ 3 ] = px[ 2 ] - p( t );
		}
};

#endif /*EHESSENBERG_H_*/
