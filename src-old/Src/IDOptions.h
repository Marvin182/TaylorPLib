/*!
 * \file IDOptions.h
 * \brief Header file for the class \file IDOptions.cpp .
 * \author D. Monett D�z
 * \date September, 2006
 * 
 * This is the header of the class for options definitions.
 * 
 * Included files:
 * 	\file IDException.h for handling exceptions.
 *  \file defs.h for global definitions.
 * 
 */
 
#ifndef IDOPTIONS_H_
#define IDOPTIONS_H_

#include "IDException.h"	/**< For handling exceptions. */
#include "defaults.h"		/**< For working with default values. */

class IDOptions
{
	public:
		IDOptions();		/**< Default constructor. */
		IDOptions( const IDOptions &op );	/**< Copy constructor. */
		~IDOptions();		/**< Destructor. */
		IDOptions daeindexset( int nameid, double value );	/**< Sets or alter options values. */
		IDOptions daeindexset( int size, int *nameid, double *value );	/**< Sets or alter options values. */
		int getDegree();	/**< Returns the degree of Taylor coefficients; highest derivative degree. */
		double getQREps();	/**< Returns the threshold for the Householder QR factorization with column pivoting. */
		int getPrintWhat();/**< Returns what to print out. */
		int getPrintWhere();/**< Returns where to print out. */
		double getIOEps();	/**< Returns the I/O threshold. */
	private:
		int degree;			/**< The degree of Taylor coefficients; highest derivative degree. */
		double qreps;		/**< The threshold for the Householder QR factorization with column pivoting. */
		int printWhat;		/**< The print out parameter what to print out. */
		int printWhere;		/**< The print out parameter where to print out. */
		double ioeps;		/**< The I/O threshold. */
};

#endif /*IDOPTIONS_H_*/
