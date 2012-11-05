/*!
 * \file IDOptions.cpp
 * \brief Class for managing global options to the method.
 * \author D. Monett Díaz
 * \date September, 2006
 * 
 * This class manages the options that can be set when determining the index in systems of 
 * ordinary differential equations (ODEs) using Automatic/Algorithmic Differentiation (AD).
 * 
 * Included files:
 * 	\file IDOptions.h header file.
 * 
 */

#include "IDOptions.h"

//not effective for arrays passed in as function parameters!!!
//#define _ARRAY_SIZE( x ) ( sizeof( x ) / sizeof( *( x )))


/**
 * Default constructor for the class. Creates the object with default values for the properties.
 * 
 */
IDOptions::IDOptions()
{
	degree = _TP_DEGREE_VAL;
	qreps = _QR_EPS_VAL;
	printWhat = _PRINT_PARAM_WHAT_VAL;
	printWhere = _PRINT_PARAM_WHERE_VAL;
	ioeps = _IO_EPS_VAL;
}

/**
 * Copy constructor.
 * 
 * \param[in] op The options object to copy from.
 * 
 */
IDOptions::IDOptions( const IDOptions &op )
{
	degree = op.degree;
	qreps = op.qreps;
	printWhat = op.printWhat;
	printWhere = op.printWhere;
	ioeps = op.ioeps;
}

/**
 * Returns the value of the degree of Taylor coefficients or highest derivative degree.
 * 
 * \return The degree.
 * 
 */
int IDOptions::getDegree()
{
	return degree;
}

/**
 * Returns the value of the threshold for the Householder QR factoritation with column pivoting.
 * 
 * \return The threshold.
 * 
 */
double IDOptions::getQREps()
{
	return qreps;
}

/**
 * Returns the value of the print out parameter that indicates what to print out.
 * 
 * \return The print out parameter value.
 * 
 */
int IDOptions::getPrintWhat()
{
	return printWhat;
}

/**
 * Returns the value of the print out parameter that indicates where to print out
 * (i.e. to a data file, to the screen, or nowhere).
 * 
 * \return The print out parameter value.
 * 
 */
int IDOptions::getPrintWhere()
{
	return printWhere;
}

/**
 * Returns the value of the I/O threshold.
 * 
 * \return The threshold.
 * 
 */
double IDOptions::getIOEps()
{
	return ioeps;
}

/**
 * Sets or alters options for the main function indexdet.
 *
 */
IDOptions IDOptions::daeindexset( int nameid, double value )
{
	try
	{
		switch( nameid )
		{
			case _TP_DEGREE_ :
				degree = (int)value;
				break;
			case _QR_EPS_ : 
				qreps = value;
				break;
			case _PRINT_PARAM_WHAT_ : 
				printWhat = (int)value;
				break;
			case _PRINT_PARAM_WHERE_ : 
				printWhere = (int)value;
				break;
			case _IO_EPS_ : 
				ioeps = value;
				break;
			default: 
				throw IDException( "Wrong property name when setting options structure.", 33 );
				break;
		}
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	
	return *this;
}

/**
 * Sets or alters options for the main function indexdet.
 *
 */
IDOptions IDOptions::daeindexset( int size, int *nameid, double *value )
{
	//printf( "\nsize = %d", size );
	try
	{
		for( int i = 0; i < size; i++ )
		{
			//printf( "\nnameid[i] = %d", nameid[ i ] );
			//printf( "\nvalue[i] = %.8lg\n", value[ i ] );
			switch( nameid[ i ] )
			{
				case _TP_DEGREE_ :
					degree = (int)value[ i ];
					break;
				case _QR_EPS_ : 
					qreps = value[ i ];
					break;
				case _PRINT_PARAM_WHAT_ : 
					printWhat = (int)value[ i ];
					break;
				case _PRINT_PARAM_WHERE_ : 
					printWhere = (int)value[ i ];
					break;
				case _IO_EPS_ : 
					ioeps = value[ i ];
					break;
				default: 
					throw IDException( "Wrong property name when setting options structure.", 33 );
					break;
			}
		}
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	
	return *this;
}

/**
 * Destructor. Cleans up the object.
 * 
 */
IDOptions::~IDOptions()
{
}
