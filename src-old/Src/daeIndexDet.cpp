/**
FILE NAME:  indexdet.cpp
WRITTEN BY: D. Monett Diaz
FOR:        DAE-AD HU-project
DUE DATE:   September 06, 2006
PURPOSE: 
	This is the main program that tests the index determination functions.
INCLUDED FILES:
	"matrixseq.h" for index determination.
	Other .h headers containing user examples (the trajectory, the dynamic, and the DAE).
	<iostream>, standard header for reading from and writing to the standard streams.
	<typeinfo>, standard header for getting information about both static and dynamic types.
**/

//
// I N D E X D E T   L Y B R A R Y   H E A D E R
//
#include "matrixseq.h"							// for the index determination as well
//
// U S E R   H E A D E R S
//
//#include "Examples/Example1.h"							// (1) user example header
//#include "Examples/EPend3.h"							// (2) user example header
//#include "Examples/ECircuit.h"							// (3) user example header
//#include "Examples/ECatMix.h"							// (4) user example header
#include "Examples/ERoboticArm.h"						// (5) user example header
//#include "Examples/EHessenberg.h"						// (6) user example header
//#include "Examples/EHessenberg1.h"						// (7) user example header
//#include "Examples/EHessenberg2.h"						// (8) user example header
//#include "Examples/ESchumilina.h"						// (9) user example header
//#include "Examples/ESchumilina1.h"						// (10) user example header
//
// S T A N D A R D   C++   L Y B R A R Y   H E A D E R S
//
#include <iostream>
#include <typeinfo>

//using namespace std;

/**
 * Main function call.
 * 
 */
int main( int argc, char *argv[] )
{
	int errcode;
	
	//
	// M O D I F Y I N G   D E F A U L T   O P T I O N S
	//
	
	IDOptions options;
	options.daeindexset( _TP_DEGREE_, 5 );		// new value for the Taylor polynomial's degree
	options.daeindexset( _QR_EPS_, 1E-10 );		// new value for the threshold
	options.daeindexset( _IO_EPS_, 1E-12 );		// new value for the I/O threshold
	options.daeindexset( _PRINT_PARAM_WHAT_, 2 );// what to print (all: 22)
	options.daeindexset( _PRINT_PARAM_WHERE_, 0 );// where to print
												// =0 to the screen; =1 to a file; =2 nowhere
	//
	// C L A S S E S   A N D   O B J E C T   D E F I N I T I O N S
	// ( U S E R   E X A M P L E S   D E P E N D E N T )
	//
	// IMPORTANT NOTE:
	// The corresponding header file should be previously included 
	// ("#include" section above)
	//
	
/*(1) Simple example */
	//Example1 ex1;								// instantiate object of class Example1.
	//IDExample * pobj = &ex1;					// pobj is a pointer to the object ex1
/*(2) Pendulum */
	//EPend3 ep3;								// instantiate object of class EPend3.
	//IDExample * pobj = &ep3;					// pobj is a pointer to the object ep3
/*(3) Electric Circuit */
	//ECircuit ecirc;							// instantiate object of class ECircuit.
	//IDExample * pobj = &ecirc;				// pobj is a pointer to the object ecirc
/*(4) ?? Example */
	//ECatMix ecm;								// instantiate object of class ECatMix.
	//IDExample * pobj = &ecm;					// pobj is a pointer to the object ecm
/*(5) Robotic Arm */
	ERoboticArm era;							// instantiate object of class ERoboticArm.
	IDExample * pobj = &era;					// pobj is a pointer to the object era
/*(6) Hessenberg DAE */
	//EHessenberg ehess;						// instantiate object of class EHessenberg.
	//IDExample * pobj = &ehess;				// pobj is a pointer to the object ehess
/*(7) Hessenberg DAE, modified (I)*/
	//EHessenberg1 ehess1;						// instantiate object of class EHessenberg1.
	//IDExample * pobj = &ehess1;				// pobj is a pointer to the object ehess1
/*(8) Hessenberg DAE, modified (II)*/
	//EHessenberg2 ehess2;						// instantiate object of class EHessenberg2.
	//IDExample * pobj = &ehess2;				// pobj is a pointer to the object ehess2
/*(9) DAE from Schumilina's thesis */
	//ESchumilina eschum;						// instantiate object of class ESchumilina.
	//IDExample * pobj = &eschum;				// pobj is a pointer to the object eschum
/*(10) DAE from Schumilina's thesis, modified (I)*/
	//ESchumilina1 eschum1;						// instantiate object of class ESchumilina1.
	//IDExample * pobj = &eschum1;				// pobj is a pointer to the object eschum1

	//
	// M A I N   F U N C T I O N   C A L L
	//
	// NOTE: uncomment just one of the following two versions
	//
	
	//
	// ...when no new parameter value is set or altered, 
	// i.e., it works with default values for the options
	//
	//errcode = daeindex( argc, argv, pobj );
	
	//
	// ...when the default value of at least one option is altered
	//
	errcode = daeindex( argc, argv, pobj, options );
	
	//
	// S O M E   S U P P O R T I N G   N O T E S
	//

	/*
	How to declare and initialize variables related to the global options:
	
	IDOptions options;							// variable of 'IDOptions' type
	int nameid[] = { _TP_DEGREE_, _QR_EPS_ };	// array of name ids of the properties to be set or altered
												// (recommended when more than one property is going to be set;
												// these identifiers are unique and are defined in "defs.h")
	double value[] = { 8, 0.00001 };			// array of new values for these properties
												// (recommended when more than one property is going to be set;
												// maintaining the same order as in the list of name identifiers
												// is crucial!!)
	
	Examples of use:
	
	options.daeindexset( _TP_DEGREE_, 12 );		// sets up one property
	options.daeindexset( _QR_EPS_, 0.001 );		// sets up one property
	options.daeindexset( 2, nameid, value );	// sets up at least one property
	
	Corresponding function calls:
	
	daeindex( argc, argv, options );			// in case at least one option was set or altered
	daeindex( argc, argv );						// otherwise, i.e., work with default values
	*/
	
	//
	// O T H E R   P O S S I B L E   F U N C T I O N   C A L L S
	//
	
	//
	// Versions without considering the functions tra, dyn and dae as parameters but hard coded
	// in the called program, library, or class using them as independent .cpp/.h files:
	//
	//errcode = daeindex( argc, argv );			// no new value is set or altered;
												// works with default values for the options
	//errcode = daeindex( argc, argv, options );	// includes new values for the options

	//
	// Versions in case no clases of examples exist, but a unique definition of tra, dyn and dae
	// as independent .cpp/.h files:
	//
	//errcode = daeindex( argc, argv, tra, dyn, dae );
												// no new value is set or altered;
												// works with default values for the options
												// functions tra, dyn, and dae as pointer functions
	//errcode = daeindex( argc, argv, tra, dyn, dae, options );
												// includes new values for the options
												// functions tra, dyn, and dae as pointer functions

	//
	// R E T U R N E D   P R O G R A M   C O D E
	//
	// (  = 0; successfull termination
	//   <> 0; exception dependent )
	//
	
	return errcode;
}
