/*!
 * \file IDException.h
 * \brief Header file for the class \file IDException.cpp .
 * \author D. Monett D�z
 * \date July, 2006
 * 
 * This is the header of the class for exception handling.
 * 
 * Included files:
 * 	\file stdio.h and \file esccolors.h for printing out.
 * 
 */
 
#ifndef IDEXCEPTION_H_
#define IDEXCEPTION_H_

#include <stdio.h>
#include "esccolors.h"

class IDException
{
	public:
		IDException();
		IDException( char* msg );
		IDException( char* msg, int err );
		IDException( const IDException& e );
		~IDException();
		void addReason( char* msg );
		void addReason( char* msg, int err );
		void report();
		int getErrCode();

	private:
		char* reason;
		int errCode;
};

#endif /*IDEXCEPTION_H_*/
