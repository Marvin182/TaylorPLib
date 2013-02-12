/*! 
* \file TPolyn.h
 * \brief Header file for the class \file TPolyn.cpp.
 * \author D. Monett Diaz
 * \date July, 2006
 * 
 * This is the header of the class for Taylor arithmetic.
 * 
 * Included files:
 * 	\file IDException.h and <stdexcept> and for handling exceptions.
 * 
 */
#ifndef TPOLYN_H_
#define TPOLYN_H_

#include <stdexcept>
#include <math.h>
#include "IDException.h"

class TPolyn
{
	//
	// Taylor polynomial with derivate degree n (n+1 coefficients):
	//
	// P_n(x) = f(a) + (x-a)*f'(a)/1! + (x-a)^2*f''(a)/2! + (x-a)^3*f'''(a)/3! 
	//			+ ... + (x-a)^n*f^{(n)}(a)/n!
	//
	//		  = sum{k=0}{n} (x-a)^k*f^{(k)}(a)/k!
	//
	// being 'f' the function to be approximated by P_n(x) at point 'a', with its
	// first n derivatives existing on a closed interval I, so that
	//
	// f(x) = P_n(x) + R_n(x)
	//
	// the remainder term being R_n(x) = (x-a)^{n+1}*f^{(n+1)}(c)/(n+1)! for some
	// 'c' between 'x' and 'a'.
	//
	// Another often used form:
	//
	// f(x0+h) = f(x0) + h*f'(x0)/1! + h^2*f''(x0)/2! + h^3*f'''(x0)/3! 
	//			+ ... + h^n*f^{(n)}(x0)/n!
	//
	private:
		int itsOrder;				// The derivative order of the Taylor polynomial
									// Nr. of coeff = itsGrade + 1
    	double *itsCoeff;			// The coefficients
    	bool _constant;				// whether it is a constant Taylor polynomial or not

	public:
		//
		// Constructors, destructor
		//
		TPolyn::TPolyn();
		TPolyn::TPolyn( int order );
		TPolyn::TPolyn( const TPolyn &tp );
		TPolyn::~TPolyn();
		//
		// Accessing properties
		//
		double* coeffs() { return itsCoeff; }				// Returns the ptr to the alloc. memory
		int order() { return itsOrder; }
		int nrcoeff() { return itsOrder + 1; }
		char* typeName() { return "TPolyn"; }
		//
		// Overloaded operators for Taylor arithmetic
		//
		double & TPolyn::operator[]( int index );
		TPolyn TPolyn::operator=( const TPolyn &tp );
		int TPolyn::operator==( const TPolyn &tp );
		int TPolyn::operator!=( const TPolyn &tp );
		int TPolyn::operator<( const TPolyn &tp );
		int TPolyn::operator<=( const TPolyn &tp );
		int TPolyn::operator>( const TPolyn &tp );
		int TPolyn::operator>=( const TPolyn &tp );
		TPolyn TPolyn::operator+( const TPolyn &w );
		TPolyn TPolyn::operator+=( const TPolyn &w );
		TPolyn TPolyn::operator-( const TPolyn &w );
		TPolyn TPolyn::operator-=( const TPolyn &w );
		TPolyn TPolyn::operator-();
		TPolyn TPolyn::operator*( const TPolyn &w );
		TPolyn TPolyn::operator*=( const TPolyn &w );
		TPolyn TPolyn::operator*( double alpha );
		TPolyn TPolyn::operator/( const TPolyn &w );
		TPolyn TPolyn::operator/=( const TPolyn &w );
		//
		// Special operators
		//
		TPolyn TPolyn::sqr();
		TPolyn TPolyn::sqrt();
		//
		// Other functions
		//
		void TPolyn::print();
		void TPolyn::print( FILE * fn );
		double TPolyn::eval( double x, double alpha );
		double TPolyn::feval();
		void TPolyn::shift();
		bool TPolyn::isConst();
		bool TPolyn::isConst( double eps );
		bool TPolyn::isId();
		bool TPolyn::isId( double eps );
		bool TPolyn::isZero();
		bool TPolyn::isZero( double eps );
		int TPolyn::set2zero();
		int TPolyn::set2zero( int ord );
		int TPolyn::set2const( double c );
		int TPolyn::setCoeffs( double *c );
};

#endif /*TPOLYN_H_*/
