/*!
 * \file IDException.cpp
 * \brief Class for handling exceptions.
 * \author D. Monett Díaz
 * \date July, 2006
 * 
 * This class handles exceptions that can occur when determining the index in systems of 
 * ordinary differential equations (ODEs) using Automatic/Algorithmic Differentiation (AD).
 * 
 * Included files:
 * 	\file IDException.h header file.
 * 
 */

#include "IDException.h"

/**
 * Default constructor for the class. Creates the object.
 * 
 */
IDException::IDException()
{
}

/**
 * Constructor for the class with a message as parameter. Creates the object.
 * 
 * \param[in] msg The message.
 * 
 */
IDException::IDException( char* msg )
{
	reason = msg;
}

/**
 * Constructor for the class with a message and its corresponding error code as parameters. Creates the object.
 * 
 * \param[in] msg The message.
 * \param[in] err The error code.
 * 
 */
IDException::IDException( char* msg, int err )
{
	reason = msg;
	errCode = err;
}

/**
 * Copy constructor.
 * 
 * \param[in] e The exception object to copy from.
 * 
 */
IDException::IDException( const IDException& e )
{
	reason = e.reason;
	errCode = e.errCode;
}

/**
 * Destructor. Cleans up the object.
 * 
 */
IDException::~IDException()
{
}

/**
 * Adds a reason for the exception.
 * 
 * \param[in] msg The reason to be added to.
 * 
 */
void IDException::addReason( char* msg )
{
	reason = msg;
}

/**
 * Adds a reason for the exception with its corresponding error code.
 * 
 * \param[in] msg The reason to be added to.
 * \param[in] err The error code.
 * 
 */
void IDException::addReason( char* msg, int err )
{
	reason = msg;
	errCode = err;
}

/**
 * Prints out a message with a corresponding error code.
 * 
 */
void IDException::report()
{
	printf( red );
	printf( "\n***The exception number %d occurred.", errCode );
	printf( "\n***%s", reason );
	printf( black );
};

/**
 * Returns the error code of the exception.
 * 
 * \return The corresponding error code.
 * 
 */
int IDException:: getErrCode()
{
	return errCode;
}
