/*!
 * \file matrixseq.cpp
 * \brief Functions related to the matrix sequence computation with matrices of Taylor polynomials.
 * \author D. Monett Diaz
 * \date January, 2007
 * 
 * This class implements the functions needed to the matrix sequence computation, included the 
 * ones related to ADOL-C's active sections that define the trajectory, the dynamic, and the DAE.
 * Matrices of Taylor polynomials are used.
 * 
 * Included files:
 * 	\file matrixseq.h header file.
 * 
 */

#include "matrixseq.h"

using namespace std;

//
// G L O B A L   V A R I A B L E   D E C L A R A T I O N S
//
IDOptions options;	/**< The global options the function 'daeindex' should work with. */
int	printWhat,		/**< The print out parameter what to print out. */
	printWhere,		/**< The variable controling where to print out. */
	deg,			/**< The degree of Taylor coefficients; highest derivative order. */
	r,				/**< The rank that is calculated. */
	success = 0;	/**< The success variable. */
double eps,			/**< The threshold for the Householder QR factoritation with column pivoting. */
	ioeps;			/**< The threshold for I/O. */
FILE * pof;			/**< The output datafile. */
//
// ...for the TRAJECTORY
//
int	m_tra,			/**< The number of dependent variables of the trajectory. */ // = m_dae!!
	n_tra;			/**< The number of independent variables of the trajectory. */ // = 1
double t;			/**< The independent variable of the trajectory; time, as double. */
double **X_tra;		/**< The vector of Taylor coefficients X[n_tra][deg+1] for the trajectory.*/
double **Y_tra;		/**< The vector of Taylor coefficients Y[m_tra][deg+1] for the trajectory. */
//DoubleMatrix X_tra,	/**< The vector of Taylor coefficients X[n_tra][deg+1] for the trajectory.*/
//	Y_tra;			/**< The vector of Taylor coefficients Y[m_tra][deg+1] for the trajectory. */
//
// ...for the DYNAMIC
//
int m_dyn,			/**< The number of dependent variables of the dynamic. */
	n_dyn;			/**< The number of independent variables of the dynamic. */ // = m_dae!!
double **Y_ddyn;	/**< The vector of Taylor coefficients Y'[m_dyn][deg+1] for the dynamic. */
//DoubleMatrix Y_ddyn;/**< The vector of Taylor coefficients Y'[m_dyn][deg+1] for the dynamic. */
double ***Z_dyn; 	/**< The vector of adjoints Z[m_dyn][n_dyn][deg+1] for the dynamic. */
//
// ...for the DAE
//
int m_dae,			/**< The number of dependent variables of the DAE. */
	n_dae;			/**< The number of independent variables of the DAE. */ // n_dae = m_dyn + m_dae
double **Y_dae; 	/**< The vector of Taylor coefficients Y[m_dae][deg+1] for the DAE. */
//DoubleMatrix Y_dae; /**< The vector of Taylor coefficients Y[m_dae][deg+1] for the DAE. */
double ***Z_dae; 	/**< The vector of adjoints Z[m_dae][n_dae][deg+1] for the DAE. */
//
// Matrices to be used by the matrix sequence
//
double ***A; 	/**< The matrix A. */
double ***B; 	/**< The matrix B. */
double ***D; 	/**< The matrix D. */

//
// M E M B E R   F U N C T I O N S
//

/**
 * Converts integer into string.
 *
 * param[in] x The integer value to be converted to.
 * 
 */
static inline string int2str( int x )
{
	std::ostringstream o;
	
	if ( !( o << x ) )
		printf( "Error in conversion from int to string \n" );
	return o.str();
}

/**
 * Process data, define variables, and initialize them.
 * 
 * Input from and output to data files are not necessary any more, since the variable passed in
 * as parameter encapsulates similar functionalities.
 * 
 */
int initialize( IDExample * pobj )
{
	int i, j, year, month, day;
	time_t now;
	tm * ptm;
	char y[ 5 ], m[ 3 ], d[ 3 ], filename[ 50 ] = "";
	
	try
	{
		//
		// DEFINE AND INITIALIZE GLOBAL VARIABLES
		//
		printWhere = options.getPrintWhere();				// print out to a file/screen/nowhere?
		deg = options.getDegree();							// degree of Taylor coefficients; 
															// highest derivative order
		eps = options.getQREps();							// threshold for the Househ. QR fact. 
															// with column pivoting
		printWhat = options.getPrintWhat();					// print out parameter to control output
		ioeps = options.getIOEps();							// threshold for I/O

		if( deg <= 0 ) 
			throw IDException( "Wrong highest derivative order (it should be greater than zero).", 42 );
		if( printWhere < 0   ||  printWhere > 2 ) 
			throw IDException( "Wrong print out parameter (where to print). Its value should be either =0 (screen), =1 (file), or =2 (nowhere).", 43 );
		if( printWhat < 0 ) 
			throw IDException( "Wrong print out parameter (what to print). Its value should not be negative.", 44 );
		if( eps < 0 ) 
			throw IDException( "Wrong threshold value. It should not be negative.", 45 );
		if( ioeps < 0 ) 
			throw IDException( "Wrong threshold for I/O. It should not be negative.", 46 );

		t = pobj->getTimeVal();								// time initialization, at starting point
		n_tra = pobj->getIndepVarTra();						// Nr. indep. var. for the trajectory
		m_tra = pobj->getDepVarTra();						// Nr. dep. var. for the trajectory
		n_dyn = pobj->getIndepVarDyn();						// Nr. indep. var. for the dynamic;
															// n_dyn = m_tra
		m_dyn = pobj->getDepVarDyn();						// Nr. dep. var. for the dynamic
		n_dae = pobj->getIndepVarDAE();						// Nr. indep. var. for the DAE;
															// n_dae = m_dyn + m_tra
		m_dae = pobj->getDepVarDAE();						// Nr. dep. var. for the DAE
		//printf( "\n\nClass name = %s", pobj->name() );	
		
		if( printWhere == 1 )
		{
			//
			// Form output file name from current date and problem class name
			//
			time( &now );
			ptm = gmtime( &now );
			year = ptm->tm_year + 1900;
			month = ptm->tm_mon + 1;
			day = ptm->tm_mday;
			
			strcpy( y, int2str( year ).c_str() );			// itoa( year, y, 10 );
			strcpy( m, int2str( month ).c_str() );			// itoa( month, m, 10 );
			strcpy( d, int2str( day ).c_str() );			// itoa( day, d, 10 );
			/*printf( "\nYear: %d\n", year );
			printf( "Month:  %d\n", month );
			printf( "Day:    %d\n", day );
			printf( "\nYear: %s\n", y );
			printf( "Month:  %s\n", m );
			printf( "Day:    %s\n", d );*/
			strcat( filename, pobj->name() );
			strcat( filename, "-" );
			strcat( filename, y );
			strcat( filename, "-" );
			strcat( filename, m );
			strcat( filename, "-" );
			strcat( filename, d );
			strcat( filename, ".log" );
			printf( "\n   (Results filename: %s)\n\n", filename );

			pof = fopen( filename, "a" );					// open output file and 
															// append data at the end of the file.
															// The file is created if it doesn't exist
			if( pof == NULL ) 
				throw IDException( "Error opening output file.", 2 );
			
			ptm = localtime( &now );
			fprintf( pof, "\n-----------------------------------------------" );
			fprintf( pof, "\nIndex determination in DAEs using AD techniques" );
			fprintf( pof, "\n-----------------------------------------------\n" );
			fprintf( pof, "\n*************************************\n\n" );
			fprintf( pof, "File name: %s\n", filename );
			fprintf( pof, "Date and time: %s\n", asctime( ptm ) );
			fprintf( pof, "Global parameters:\n\n" );
			fprintf( pof, "     deg = %d (Degree of Taylor coefficients)\n", deg );
			fprintf( pof, "     eps = %.16lg (Threshold for the Householder QR factorization with column pivoting)\n", eps );
			fprintf( pof, "     ioeps = %.16lg (Threshold for I/O)\n", ioeps );
			fprintf( pof, "     printWhat = %d (Print out parameter for controlling output)\n\n", printWhat );
			fprintf( pof, "*************************************\n\n" );
			fprintf( pof, "Initial time and problem dimensions:\n" );
			fprintf( pof, "     t = %.16lg (Independent var. time)\n", t );
			fprintf( pof, "     n_tra = %d (Nr. indep. var. trajectory)\n", n_tra );
			fprintf( pof, "     m_tra = %d (Nr. dep. var. trajectory)\n", m_tra );
			fprintf( pof, "     n_dyn = %d (Nr. indep. var. dynamic)\n", n_dyn );
			fprintf( pof, "     m_dyn = %d (Nr. dep. var. dynamic)\n", m_dyn );
			fprintf( pof, "     n_dae = %d (Nr. indep. var. DAE)\n", n_dae );
			fprintf( pof, "     m_dae = %d (Nr. dep. var. DAE)\n\n", m_dae );
			fprintf( pof, "*************************************\n" );
		}
		else if( printWhere == 0 )
		{
			printf( cyan );
			printf( "\n-----------------------------------------------" );
			printf( "\nIndex determination in DAEs using AD techniques" );
			printf( "\n-----------------------------------------------\n\n" );
			printf( normal );
			if( printWhat > 0 )
			{
				printf( black );
				printf( "deg = %d", deg );
				printf( "\neps = %.16lg", eps );
				printf( "\nioeps = %.16lg", ioeps );
				printf( "\nprintWhat = %d", printWhat );
				printf( "\nt = %.16lg", t );
				printf( "\nn_tra = %d", n_tra );
				printf( "\nm_tra = %d", m_tra );
				printf( "\nn_dyn = %d", n_dyn );
				printf( "\nm_dyn = %d", m_dyn );
				printf( "\nn_dae = %d", n_dae );
				printf( "\nm_dae = %d\n", m_dae );
				printf( normal );
			}
		}
		//
		// ...for the TRAJECTORY
		//
		//X_tra.redim( n_tra, deg + 1 );						// Taylor coeff. X[n_tra][deg+1]
		X_tra = new double*[ n_tra ];						// Taylor coeff. X[n_tra][deg+1]
		if( (X_tra == 0)  ||  (X_tra == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < n_tra; i++ )
		{
			X_tra[ i ] = new double[ deg + 1 ];
			if( (X_tra[ i ] == 0)  ||  (X_tra[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}
		//Y_tra.redim( m_tra, deg + 1 );						// Taylor coeff. Y[m_tra][deg+1]
		Y_tra = new double*[ m_tra ];						// Taylor coeff. Y[m_tra][deg+1]
		if( (Y_tra == 0)  ||  (Y_tra == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_tra; i++ )
		{
			Y_tra[ i ] = new double[ deg + 1 ];
			if( (Y_tra[ i ] == 0)  ||  (Y_tra[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}
		//
		// ...for the DYNAMIC and also matrix D
		//
		//Y_ddyn.redim( m_dyn, deg + 1 );						// Taylor coeff. Y'[m_dyn][deg+1]
		Y_ddyn = new double*[ m_dyn ];						// Taylor coeff. Y'[m_dyn][deg+1]
		if( (Y_ddyn == 0)  ||  (Y_ddyn == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dyn; i++ )
		{
			Y_ddyn[ i ] = new double[ deg + 1 ];
			if( (Y_ddyn[ i ] == 0)  ||  (Y_ddyn[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}
		Z_dyn = new double**[ m_dyn ];						// adjoints Z[m_dyn][n_dyn+1][deg+1]
		if( (Z_dyn == 0)  ||  (Z_dyn == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dyn; i++ )
		{
			Z_dyn[ i ] = new double*[ n_dyn + 1 ];
			if( (Z_dyn[ i ] == 0)  ||  (Z_dyn[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
			for( j = 0; j < n_dyn + 1; j++ )
			{
				Z_dyn[ i ][ j ] = new double[ deg + 1 ];
				if( (Z_dyn[ i ][ j ] == 0)  ||  (Z_dyn[ i ][ j ] == NULL) )
					throw IDException( "Memory allocation failure.", 3 );
			}
		}
		D = new double**[ m_dyn ];							// matrix D
		if( (D == 0)  ||  (D == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dyn; i++ )
		{
			D[ i ] = new double*[ n_dyn ];
			if( (D[ i ] == 0)  ||  (D[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
			for( j = 0; j < n_dyn; j++ )
			{
				D[ i ][ j ] = new double[ deg ];
				if( (D[ i ][ j ] == 0)  ||  (D[ i ][ j ] == NULL) )
					throw IDException( "Memory allocation failure.", 3 );
			}
		}
		//
		// ...for the DAE and also matrices A and B
		//
		//Y_dae.redim( m_dae, deg + 1 );						// Taylor coeff. Y[m_dae][deg+1]
		Y_dae = new double*[ m_dae ];						// Taylor coeff. Y[m_dae][deg+1]
		if( (Y_dae == 0)  ||  (Y_dae == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dae; i++ )
		{
			Y_dae[ i ] = new double[ deg + 1 ];
			if( (Y_dae[ i ] == 0)  ||  (Y_dae[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}
		Z_dae = new double**[ m_dae ];						// adjoints Z[m_dae][n_dae+1][deg+1]
		if( (Z_dae == 0)  ||  (Z_dae == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dae; i++ )
		{
			Z_dae[ i ] = new double*[ n_dae + 1 ];
			if( (Z_dae[ i ] == 0)  ||  (Z_dae[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
			for( j = 0; j < n_dae + 1; j++ )
			{
				Z_dae[ i ][ j ] = new double[ deg + 1 ];
				if( (Z_dae[ i ][ j ] == 0)  ||  (Z_dae[ i ][ j ] == NULL) )
					throw IDException( "Memory allocation failure.", 3 );
			}
		}
		A = new double**[ m_dae ];						// matrix A
		if( (A == 0)  ||  (A == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dae; i++ )
		{
			A[ i ] = new double*[ m_dyn ];
			if( (A[ i ] == 0)  ||  (A[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
			for( j = 0; j < m_dyn; j++ )
			{
				A[ i ][ j ] = new double[ deg ];
				if( (A[ i ][ j ] == 0)  ||  (A[ i ][ j ] == NULL) )
					throw IDException( "Memory allocation failure.", 3 );
			}
		}
		B = new double**[ m_dae ];						// matrix B
		if( (B == 0)  ||  (B == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dae; i++ )
		{
			B[ i ] = new double*[ m_tra ];
			if( (B[ i ] == 0)  ||  (B[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
			for( j = 0; j < m_tra; j++ )
			{
				B[ i ][ j ] = new double[ deg ];
				if( (B[ i ][ j ] == 0)  ||  (B[ i ][ j ] == NULL) )
					throw IDException( "Memory allocation failure.", 3 );
			}
		}
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
		if( printWhere == 1 ) 
			fclose( pof );									// Close data files
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
		if( printWhere == 1 ) 
			fclose( pof );									// Close data files
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		if( printWhere == 1 ) 
			fclose( pof );									// Close data files
		throw e.getErrCode();
	}
	catch(...)												// others, including the rethrowed ones
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		if( printWhere == 1 ) 
			fclose( pof );									// Close data files
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Performs global cleaning operations, in a last-in-first-out fashion.
 * 
 */
int cleanup()
{	
	int i, j;
	
	try
	{
		printf( "\nPerforming cleaning operations...\n" );
		for( i = 0; i < m_dae; i++ )
		{
			for( j = 0; j < m_tra; j++ )
				delete [] B[ i ][ j ];
			delete [] B[ i ];
		}
		delete [] B;
		//printf( "done for B...\n" );
		for( i = 0; i < m_dae; i++ )
		{
			for( j = 0; j < m_dyn; j++ )
				delete [] A[ i ][ j ];
			delete [] A[ i ];
		}
		delete [] A;
		//printf( "done for A...\n" );
		/*for( i = 0; i < m_dae; i++ )
		{
			for( j = 0; j < n_dae + 1; j++ )
				delete [] Z_dae[ i ][ j ];
			delete [] Z_dae[ i ];
		}*/
		delete [] Z_dae;
		//printf( "done for Z_dae...\n" );
		//rintf( "done for Y_dae...\n" );
		/*for( i = 0; i < m_dyn; i++ )
		{
			for( j = 0; j < n_dyn; j++ )
				delete [] D[ i ][ j ];
			delete [] D[ i ];
		}*/
		delete [] D;
		//printf( "done for D...\n" );
		/*for( i = 0; i < m_dyn; i++ )
		{
			for( j = 0; j < n_dyn + 1; j++ )
				delete [] Z_dyn[ i ][ j ];
			delete [] Z_dyn[ i ];
		}*/
		delete [] Z_dyn;
		//printf( "done for Z_dyn...\n" );
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// all other exceptions
	{
		IDException e( "Error when deallocating global variables.", 41 );
        e.report();
		throw e.getErrCode();
	}
	printf( "Done.\n" );
	return 0;
}

/**
 * Shift operator to calculate the derivative of Taylor polynomials.
 * 
 * y(t) = sum_{j=0}^{d} y_j * t^j + O(t^d+1)
 * 		= y_0 + y_1*t + y_2*t^2 + ... + y_d*t^d
 * 
 * y'(t) = y_1 + 2*y_2*t + 3*y_3*t^2 + ... + d*y_d*t^d-1
 * 
 * The rows of the input matrix represent the different Taylor polinomials.
 * The columns of the input matrix represent the coefficients of the corresponding 
 * Taylor polinomial.
 * 
 * The output matrix has the same number of columns as the input matrix. Its last column, 
 * however, contains only zeros. The resting coefficients are shifted to the left.
 * 
 * \param[in] m The number of rows of the matrix.
 * \param[in] n The number of columns of the matrix.
 * \param[in] M1 The pointer to \a M1, a matrix of \type doubles.
 * \param[out] M2 The pointer to \a M2, the resulting matrix.
 *
 */
void shift( int m, int n, double **M1, double **M2 )
{
	for( int i = 0; i < m; i++ )
	{
		//for( int j = 1; j < n - 1; j++ )
		//	M2[ i ][ j - 1 ] = ( j - 1 )* M1[ i ][ j ];
		for( int j = 1; j < n; j++ )
			M2[ i ][ j - 1 ] = j * M1[ i ][ j ];
		M2[ i ][ n - 1 ] = 0.0;								// Last column is zeroed
	}	
}

/**
 * Defines the active section for the trajectory.
 * 
 * Calculates x_* for a given t, where:
 * 
 * 	x_*(t): R in Rn (see tra.cpp/.h)
 * 
 * \param[in] pobj The pointer to the object that encapsulates the user example.
 * 
 */
int asectra( IDExample *pobj )
{
	int i;
	adouble at;			/**< The indep. variable of the trajectory; time, as adouble. */
	double *px;			/**< The pointer to px, the dep. vector of doubles for the trajectory. */
	adouble *pax;		/**< The pointer to pax, the dep. vector of adoubles for the trajectory. */
	double *indep_tra;	/**< The indep. variable of the trajectory; time, as a vector of doubles. */
	double *dep_tra;	/**< The dep. variable of the trajectory as a vector of doubles. */
	
	try
	{
		px = new double[ m_tra ];							// dep. vector with doubles 
		pax = new adouble[ m_tra ];							// dep. vector with adoubles
		indep_tra = new double[ n_tra ];					// indep. vector with doubles
		dep_tra = new double[ m_tra ];						// dep. vector with doubles
		if( (px == 0)  ||  (px == NULL)  ||  (pax == 0)  ||  (pax == NULL)  ||
			(indep_tra == 0)  ||  (indep_tra == NULL)  ||  (dep_tra == 0)  ||  (dep_tra == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
	
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\n-------------------" );
				fprintf( pof, "\nT R A J E C T O R Y" );
				fprintf( pof, "\n-------------------" );
			}
			else if( printWhere == 0 ) 
			{
				printf( cyan );
				printf( "\n-------------------" );
				printf( "\nT R A J E C T O R Y" );
				printf( "\n-------------------" );
				printf( normal );
			}
		}
		//
		// ACTIVE SECTION 0
		//
		trace_on( 0 );										// mark active section for the 
															// trajectory (tag = 0)
		at <<= t;											// 1 indep. var.: 
															// at is active, t is passive 
		(pobj->tra)( at, pax );								// compute the trajectory via 
															// function pointer to 'tra'
		for( i = 0; i < m_tra; i++ )
			pax[ i ] >>= px[ i ];							// m_tra dep. var.: pax is active, 
															//px is passive
		if( printWhat > 1 )
		{
			if( printWhere == 1 ) 
				fprintdv( pof, m_tra, px, "\nx(t) =\n" );
			else if( printWhere == 0 ) 
				printdv( m_tra, px, "\nx(t) =\n", black );
		}
		trace_off();										// end active section
															// trace_off( 1 ); if save to file
		//
		// INITIALIZE TAYLOR COEFFICIENTS
		//
		/*X_tra.set2val( 0, 0, at.value() );				// the trajectory depends on...
		X_tra.set2val( 0, 1, 1.0 );							// ...one parameter only
		for( i = 2; i < deg + 1; i++ )
			X_tra.set2val( 0, i, 0.0 );*/
		X_tra[ 0 ][ 0 ] = at.value();						// the trajectory depends on...
		X_tra[ 0 ][ 1 ] = 1.0 ;								// ...one parameter only
		for( i = 2; i < deg + 1; i++ )
			X_tra[ 0 ][ i ] = 0.0 ;
		//
		// TEST:  int function( tag, m, n, x, y )...
		//
		/*
		indep_tra[ 0 ] = at.value();						// ...with the indep. var. as a vector
		function( 0, m_tra, n_tra, indep_tra, dep_tra );	// tag = 0
		printdv( m_tra, dep_tra, "\n...testing 'function', y[1] = dependent vector \n" );
		*/
		// Forward mode, tag = 0
		forward( 0, m_tra, n_tra, deg, deg+1, X_tra, Y_tra );
		//forward( 0, m_tra, n_tra, deg, deg+1, X_tra.data(), Y_tra.data() );
		if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				fprintdm( pof, n_tra, deg + 1, X_tra, "\nX_tra[n_tra][d+1] (Taylor coeff. X of indep. var. t) = \n" );
				fprintdm( pof, m_tra, deg + 1, Y_tra, "\nY_tra[m_tra][d+1] (Taylor coeff. Y of dep. var. x(t)) = \n" );
				//X_tra.fprintdm( pof, "\nX_tra[n_tra][d+1] (Taylor coeff. X of indep. var. t) = \n" );
				//Y_tra.fprintdm( pof, "\nY_tra[m_tra][d+1] (Taylor coeff. Y of dep. var. x(t)) = \n" );
			}
			else if( printWhere == 0 )
			{
				printdm( n_tra, deg + 1, X_tra, "\nX_tra[n_tra][d+1] (Taylor coeff. X of indep. var. t) = \n", blue );
				printdm( m_tra, deg + 1, Y_tra, "\nY_tra[m_tra][d+1] (Taylor coeff. Y of dep. var. x(t)) = \n", blue );
				//X_tra.printdm( "\nX_tra[n_tra][d+1] (Taylor coeff. X of indep. var. t) = \n", blue );
				//Y_tra.printdm( "\nY_tra[m_tra][d+1] (Taylor coeff. Y of dep. var. x(t)) = \n", blue );
			}
		}
		delete [] dep_tra;
		delete [] indep_tra;
		delete [] pax;
		delete [] px;
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including rethrowed ones 
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Defines the active section for the dynamic.
 * 
 * Calculates d for given x(t), where:
 * 
 * 		d( x(t), t ): Rn+1 in Rm (see dyn.cpp/.h)
 * 
 * 		D = d_x = sum_{j=0}^{d} D_j * t^j  in Rmxn
 * 
 * \param[in] m_dyn The number of dependent variables of the dynamic.
 * \param[in] n_dyn The number of independent variables of the dynamic.
 * \param[in] pax The pointer to pax, the dependent vector of \type adoubles for the trajectory.
 * \param[in] t The independent variable of the trajectory; time, as a \type double.
 * \param[out] pad The pointer to pad, the dependent vector of \type adoubles for the dynamic.
 * \param[in] Y_tra The vector of Taylor coefficients \a Y[m_tra][deg+1] for the trajectory.
 * \param[out] Y_ddyn The vector of Taylor coefficients \a Y[m_dyn][deg+1] for the dynamic.
 * 
 */
int asecdyn( IDExample *pobj )
{
	int i, j;
	adouble atd;		/**< The indep. variable of the dynamic; time, as adouble. */
	double *pd;			/**< The pointer to pd, the dep. vector of doubles for the dynamic. */
	adouble *pad;		/**< The pointer to pad, the dep. vector of adoubles for the dynamic. */
	adouble *padyn;		/**< The pointer to padyn, the dep. vector of adoubles for the dynamic. */
	double *dep_dyn;	/**< The dep. variable of the dynamic as a vector of doubles. */
	double **X_dyn;		/**< The vector of Taylor coefficients X[n_dyn+1][deg+1] for the dynamic.*/
	double **Y_dyn; 	/**< The vector of Taylor coefficients Y[m_dyn][deg+1] for the dynamic. */
	//DoubleMatrix X_dyn,	/**< The vector of Taylor coefficients X[n_dyn+1][deg+1] for the dynamic.*/
	//	Y_dyn; 			/**< The vector of Taylor coefficients Y[m_dyn][deg+1] for the dynamic. */
	short **nz_dyn;		/**< The nonzero pattern nz[m_dyn][n_dyn] of Z for the dynamic. */

	try
	{
		pd = new double[ m_dyn ];							// dep. vector with doubles
		pad = new adouble[ n_dyn ];							// dep. vector with adoubles
		padyn = new adouble[ m_dyn ];
		dep_dyn = new double[ m_dyn ];						// dep. vector with doubles
		if( (pd == 0)  ||  (pd == NULL)  ||  (pad == 0)  ||  (pad == NULL)  ||  
			(padyn == 0)  ||  (padyn == NULL)  ||  (dep_dyn == 0)  ||  (dep_dyn == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		//X_dyn.redim( n_dyn + 1, deg + 1 );					// Taylor coeff. X[n_dyn+1][deg+1]
		X_dyn = new double*[ n_dyn + 1 ];					// Taylor coeff. X[n_dyn+1][deg+1]
		if( (X_dyn == 0)  ||  (X_dyn == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < n_dyn + 1; i++ )
		{
			X_dyn[ i ] = new double[ deg + 1 ];
			if( (X_dyn[ i ] == 0)  ||  (X_dyn[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}
		//Y_dyn.redim( m_dyn, deg + 1 );						// Taylor coeff. Y[m_dyn][deg+1]
		Y_dyn = new double*[ m_dyn ];						// Taylor coeff. Y[m_dyn][deg+1]
		if( (Y_dyn == 0)  ||  (Y_dyn == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dyn; i++ )
		{
			Y_dyn[ i ] = new double[ deg + 1 ];
			if( (Y_dyn[ i ] == 0)  ||  (Y_dyn[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}
		nz_dyn = new short*[ m_dyn ];						// nonzero pattern nz[m_dyn][n_dyn+1] of Z
		if( (nz_dyn == 0)  ||  (nz_dyn == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dyn; i++ )
		{
			nz_dyn[ i ] = new short[ n_dyn + 1 ];
			if( (nz_dyn[ i ] == 0)  ||  (nz_dyn[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}

		if( printWhat > 1 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\n-------------" );
				fprintf( pof, "\nD Y N A M I C" );
				fprintf( pof, "\n-------------" );
			}
			else if( printWhere == 0 )
			{
				printf( cyan );
				printf( "\n-------------" );
				printf( "\nD Y N A M I C" );
				printf( "\n-------------" );
				printf( normal );
			}
		}
		//
		// ACTIVE SECTION 1
		//
		trace_on( 1 );										// mark active section for the 
															// dynamic (tag = 1)
		for( i = 0; i < n_dyn; i++ )
			pad[ i ] <<= Y_tra[ i ][ 0 ];						// independent active vector
			//pad[ i ] <<= Y_tra( i, 0 );
		atd <<= t;											// t is indep. active var.
		/*if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\nDynamic \n" );
				for( i = 0; i < n_dyn; i++ )
					fprintf( pof, "%c(%c%.16lg%c)\n", '\t','\t', pad[ i ].value(), '\t' );
				fprintf( pof, "%c(%c%.16lg%c)\n", '\t','\t', atd.value(), '\t' );
			}
			else if( printWhere == 0 )
			{
				printf( black );
				//printf( "\natd=%.16lg\n", atd.value() );
				//printf( "\n-------n_dyn = %d \n", n_dyn );
				printf( "\nDynamic \n" );
				for( i = 0; i < n_dyn; i++ )
					printf( "%c(%c%.16lg%c)\n", '\t','\t', pad[ i ].value(), '\t' );
				printf( "%c(%c%.16lg%c)\n", '\t','\t', atd.value(), '\t' );
				printf( normal );
			}
		}*/
		(pobj->dyn)( pad, atd, padyn );						// compute dynamic via 
															// function pointer to 'dyn'
		//
		// Active variables: Dependants
		//
		for( i = 0; i < m_dyn; i++ )
			padyn[ i ] >>= pd[ i ];							// m_dyn dep. var.: 
															// padyn is active, pd is passive
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintdv( pof, m_dyn, pd, "\nd(x(t), t) =\n" );
			else if( printWhere == 0 )
				printdv( m_dyn, pd, "\nd(x(t), t) =\n", black );
		}
		trace_off();										// end active section
															// trace_off( 1 ); if save to file

		//
		// INITIALIZE TAYLOR COEFFICIENTS
		//
		// X_dyn = [ Y_tra, X_tra ]  !!!!!
		//
		//for( j = 0; j < deg + 1; j++ )
		//{
		//	for( i = 0; i < n_dyn; i++ )					// first: Taylor coef. of x(t)
		//		X_dyn.set2val( i, j, Y_tra( i, j ) );
		//	X_dyn.set2val( n_dyn, j, X_tra( 0, j ) );		// second: Taylor coef. of t
		//}
		for( j = 0; j < deg + 1; j++ )
		{
			for( i = 0; i < n_dyn; i++ )					// first: Taylor coef. of x(t)
				X_dyn[ i ][ j ] = Y_tra[ i ][ j ];
			X_dyn[ n_dyn ][ j ] = X_tra[ 0 ][ j ];			// second: Taylor coef. of t
		}
		//printdm( n_dyn + 1, deg + 1, X_dyn, "\nX_dyn[n_dyn+1][d+1] -----after Taylor init.----- = \n", cyan );
		// Forward mode, tag = 1
		forward( 1, m_dyn, n_dyn + 1, deg, deg+1, X_dyn, Y_dyn );
		//forward( 1, m_dyn, n_dyn + 1, deg, deg+1, X_dyn.data(), Y_dyn.data() );
		reverse( 1, m_dyn, n_dyn + 1, deg, Z_dyn, nz_dyn );
		shift( m_dyn, deg + 1, Y_dyn, Y_ddyn );				// Calculate d' => SHIFT OPERATOR on Y_dyn
		//shift( m_dyn, deg + 1, Y_dyn.data(), Y_ddyn.data() );// Calculate d' => SHIFT OPERATOR on Y_dyn

		if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				fprintdm( pof, n_dyn + 1, deg + 1, X_dyn, "\nX_dyn[n_dyn+1][d+1] (Taylor coeff. X of indep. var. x(t) and t) = \n" );
				fprintdm( pof, m_dyn, deg + 1, Y_dyn, "\nY_dyn[m_dyn][d+1] (Taylor coeff. Y of dep. var. d(x(t),t) = \n" );
				fprint3ddm( pof, m_dyn, n_dyn + 1, deg + 1, Z_dyn, "\nZ_dyn[m_dyn][n_dyn+1][d] (Adjoints) = \n" );
				fprintdm( pof, m_dyn, deg + 1, Y_ddyn, "\nY_ddyn[m_dyn][d+1] after shift op. (Taylor coeff. of dep. var. d'(x(t),t) = \n" );
			}
			else if( printWhere == 0 )
			{
				printdm( n_dyn + 1, deg + 1, X_dyn, "\nX_dyn[n_dyn+1][d+1] (Taylor coeff. X of indep. var. x(t) and t) = \n", blue );
				printdm( m_dyn, deg + 1, Y_dyn, "\nY_dyn[m_dyn][d+1] (Taylor coeff. Y of dep. var. d(x(t),t) = \n", blue );
				//X_dyn.printdm( "\nX_dyn[n_dyn+1][d+1] (Taylor coeff. X of indep. var. x(t) and t) = \n", blue );
				//Y_dyn.printdm( "\nY_dyn[m_dyn][d+1] (Taylor coeff. Y of dep. var. d(x(t),t) = \n", blue );
				print3ddm( m_dyn, n_dyn + 1, deg + 1, Z_dyn, "\nZ_dyn[m_dyn][n_dyn+1][d] (Adjoints) = \n", green );
				printdm( m_dyn, deg + 1, Y_ddyn, "\nY_ddyn[m_dyn][d+1] after shift op. (Taylor coeff. of dep. var. d'(x(t),t) = \n", blue );
				//Y_ddyn.printdm( "\nY_ddyn[m_dyn][d+1] after shift op. (Taylor coeff. of dep. var. d'(x(t),t) = \n", blue );
			}
		}

		delete [] dep_dyn;									// Deallocating local variables...
		delete [] padyn;
		delete [] pad;	
		delete [] pd;	
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including the rethrowed ones
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Defines the active section for the DAE.
 * 
 * Calculates the DAE for given y(t), x(t), and t where:
 * 
 * 		f( y(t), x(t), t ): Rm+n+1 in Rn (see dae.cpp/.h)
 * 
 * 		A = f_y = sum_{j=0}^{d} A_j * t^j  in Rnxm
 * 		B = f_x = sum_{j=0}^{d} B_j * t^j  in Rnxn
 * 
 * \param[in] m_dae The number of dependent variables of the DAE.
 * \param[in] n_dae The number of independent variables of the DAE.
 * \param[in] pax The pointer to pax, the dependent vector of \type adoubles for the trajectory.
 * \param[in] pad The pointer to pad, the dependent vector of \type adoubles for the dynamic.
 * \param[in] t The independent variable of the trajectory; time, as a \type double.
 * \param[in] Y_ddyn The vector of Taylor coefficients \a Y[m_dyn][deg+1] for the dynamic.
 * \param[out] Y_dae The vector of Taylor coefficients \a Y[m_dae][deg+1] for the DAE.
 * \param[out] Z_dae The vector of adjoints \a Z[m_dae][n_dae][deg+1] for the DAE.
 * 
 */
int asecdae( IDExample *pobj )
{
	int i, j;
	adouble atdae;		/**< The indep. variable of the DAE; time, as adouble. */
	double *pf;			/**< The pointer to pf, the dep. vector of doubles for the DAE. */
	adouble *paf;		/**< The pointer to paf, the dep. vector of adoubles for the DAE. */
	adouble *padae1;	/**< The pointer to padae1, the first indep. vector of adoubles for the DAE. */
	adouble *padae2;	/**< The pointer to padae2, the second indep. vector of adoubles for the DAE. */
	double *padae;		/**< The pointer to padae, the indep. vector of doubles for the DAE. */
	double *dep_dae;	/**< The dependent variable of the DAE as a vector of doubles. */
	double **X_dae; 	/**< The vector of Taylor coefficients X[n_dae][deg+1] for the DAE.*/
	//DoubleMatrix X_dae; /**< The vector of Taylor coefficients X[n_dae][deg+1] for the DAE.*/
	short **nz_dae;		/**< The nonzero pattern nz[m_dae][n_dae] of Z for the DAE. */

	try
	{
		pf = new double[ m_dae ];							// dep. vector with doubles
		paf = new adouble[ m_dae ];							// dep. vector with adoubles
		padae1 = new adouble[ m_dyn ];
		padae2 = new adouble[ m_tra ];
		padae = new double[ n_dae ];						// indep. vector with doubles
		dep_dae = new double[ m_dae ];						// dep. vector with doubles
		if( (pf == 0)  ||  (pf == NULL)  ||  (paf == 0)  ||  (paf == NULL)  ||
			(padae1 == 0)  ||  (padae1 == NULL)  ||  (padae2 == 0)  ||  (padae2 == NULL)  ||
			(padae == 0)  ||  (padae == NULL)  ||  (dep_dae == 0)  ||  (dep_dae == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		//X_dae.redim( n_dae + 1, deg + 1 );					// Taylor coeff. X[n_dae+1][deg+1]
		X_dae = new double*[ n_dae + 1 ];					// Taylor coeff. X[n_dae+1][deg+1]
		if( (X_dae == 0)  ||  (X_dae == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < n_dae + 1; i++ )
		{
			X_dae[ i ] = new double[ deg + 1 ];
			if( (X_dae[ i ] == 0)  ||  (X_dae[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}
		nz_dae = new short*[ m_dae ];						// nonzero pattern nz[m_dae][n_dae+1] of Z
		if( (nz_dae == 0)  ||  (nz_dae == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < m_dae; i++ )
		{
			nz_dae[ i ] = new short[ n_dae + 1 ];
			if( (nz_dae[ i ] == 0)  ||  (nz_dae[ i ] == NULL) )
				throw IDException( "Memory allocation failure.", 3 );
		}

		if( printWhat > 1 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\n-----" );
				fprintf( pof, "\nD A E" );
				fprintf( pof, "\n-----" );
			}
			else if( printWhere == 0 )
			{
				printf( cyan );
				printf( "\n-----" );
				printf( "\nD A E" );
				printf( "\n-----" );
				printf( normal );
			}
		}
		//
		// ACTIVE SECTION 2
		//
		trace_on( 2 );										// mark active section for the DAE 
															// (tag = 2)
		for( i = 0; i < m_dyn; i++ )
			padae1[ i ] <<= Y_ddyn[ i ][ 0 ];				// first part of indep. vector
			//padae1[ i ] <<= Y_ddyn( i, 0 );
		for( i = 0; i < m_tra; i++ )
			padae2[ i ] <<= Y_tra[ i ][ 0 ];				// second part of indep. vector
			//padae2[ i ] <<= Y_tra( i, 0 );
		atdae <<= t;										// third part of indep. vector
		/*if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\nIndependents \n" );
				for( i = 0; i < m_dyn; i++ )
					fprintf( pof, "%c(%c%.16lg%c)\n", '\t','\t', padae1[ i ].value(), '\t' );
				for( i = 0; i < m_tra; i++ )
					fprintf( pof, "%c(%c%.16lg%c)\n", '\t','\t', padae2[ i ].value(), '\t' );
				fprintf( pof, "%c(%c%.16lg%c)\n", '\t','\t', atdae.value(), '\t' );
			}
			else if( printWhere == 0 )
			{
				printf( black );
				printf( "\nIndependents \n" );
				for( i = 0; i < m_dyn; i++ )
					printf( "%c(%c%.16lg%c)\n", '\t','\t', padae1[ i ].value(), '\t' );
				for( i = 0; i < m_tra; i++ )
					printf( "%c(%c%.16lg%c)\n", '\t','\t', padae2[ i ].value(), '\t' );
				printf( "%c(%c%.16lg%c)\n", '\t','\t', atdae.value(), '\t' );
				printf( normal );
			}
		}*/
		(pobj->dae)( padae1, padae2, atdae, paf );			// compute the DAE via function 
															// pointer to 'dae'
		//
		// Active variables: Dependants
		//
		for( i = 0; i < m_dae; i++ )
			paf[ i ] >>= pf[ i ];							// m_dae dep. var.: 
															// paf is active, pf is passive
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintdv( pof, m_dae, pf, "\nf(d'(x(t),t), x(t), t) =\n" );
			else if( printWhere == 0 )
				printdv( m_dae, pf, "\nf(d'(x(t),t), x(t), t) =\n", black );
		}
		trace_off();										// end active section
															// trace_off( 1 ); if save to file
		//
		// INITIALIZE TAYLOR COEFFICIENTS
		//
		// X_dae = [ Y_ddyn, Y_tra, X_tra ]  !!!!!
		//
		//for( j = 0; j < deg + 1; j++ )
		//{
		//	for( i = 0; i < m_dyn; i++ )					// first: m_dyn components of Y_ddyn
		//		X_dae.set2val( i, j, Y_ddyn( i, j ) );
		//	for( i = m_dyn; i < n_dae; i++ )				// second: m_tra components of Y_tra
		//		X_dae.set2val( i, j, Y_tra( i - m_dyn, j ) );
		//	X_dae.set2val( n_dae, j, X_tra( 0, j ) );		// third: Taylor coef. of t
		//}
		for( j = 0; j < deg + 1; j++ )
		{
			for( i = 0; i < m_dyn; i++ )					// first: m_dyn components of Y_ddyn
				X_dae[ i ][ j ] = Y_ddyn[ i ][ j ];
				//X_dae.set2val( i, j, Y_ddyn( i, j ) );
			for( i = m_dyn; i < n_dae; i++ )				// second: m_tra components of Y_tra
				X_dae[ i ][ j ] = Y_tra[ i - m_dyn ][ j ];
				//X_dae.set2val( i, j, Y_tra( i - m_dyn, j ) );
			X_dae[ n_dae ][ j ] = X_tra[ 0 ][ j ];			// third: Taylor coef. of t
			//X_dae.set2val( n_dae, j, X_tra( 0, j ) );		// third: Taylor coef. of t
		}
		//
		//
		// TEST:  int function( tag, m, n, x, y )...
		//
		/*
		for( i = 0; i < m_dyn; i++ )						// y = [ z[], x[] ]
			padae[ i ] = padae1[ i ].value();
		for( i = m_dyn; i < n_dae; i++ )
			padae[ i ] = padae2[ i - m_dyn ].value();
		function( 2, m_dae, n_dae, padae, dep_dae );		// tag = 2
		printdv( m_dae, dep_dae, "\n...testing 'function', y[1] = dependent vector \n" );
		*/
		// Forward mode, tag = 2
		forward( 2, m_dae, n_dae + 1, deg, deg+1, X_dae, Y_dae );
		//forward( 2, m_dae, n_dae + 1, deg, deg+1, X_dae.data(), Y_dae.data() );
		reverse( 2, m_dae, n_dae + 1, deg, Z_dae, nz_dae );

		if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				fprintdm( pof, n_dae + 1, deg + 1, X_dae, "\nX_dae[n_dae+1][d+1] (Taylor coefficients X of indep. var. [z(t) x(t)]) = \n" );
				fprintdm( pof, m_dae, deg + 1, Y_dae, "\nY_dae[m_dae][d+1] (Taylor coefficients Y of dep. var. f(z(t),x(t),t) = \n" );
				fprint3ddm( pof, m_dae, n_dae + 1, deg + 1, Z_dae, "\nZ_dae[m_dae][n_dae+1][d] (Adjoints) = \n" );
			}
			else if( printWhere == 0 )
			{
				printdm( n_dae + 1, deg + 1, X_dae, "\nX_dae[n_dae+1][d+1] (Taylor coefficients X of indep. var. [z(t) x(t)]) = \n", blue );
				printdm( m_dae, deg + 1, Y_dae, "\nY_dae[m_dae][d+1] (Taylor coefficients Y of dep. var. f(z(t),x(t),t) = \n", blue );
				//X_dae.printdm( "\nX_dae[n_dae+1][d+1] (Taylor coefficients X of indep. var. [z(t) x(t)]) = \n", blue );
				//Y_dae.printdm( "\nY_dae[m_dae][d+1] (Taylor coefficients Y of dep. var. f(z(t),x(t),t) = \n", blue );
				print3ddm( m_dae, n_dae + 1, deg + 1, Z_dae, "\nZ_dae[m_dae][n_dae+1][d] (Adjoints) = \n", green );
			}
		}

		delete [] dep_dae;									// Deallocating local variables...
		delete [] padae;
		delete [] padae2;	
		delete [] padae1;	
		delete [] paf;	
		delete [] pf;	
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including the rethrowed ones
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Constructs the matrices A, B, and D that will be used by the matrix sequence.
 * 
 */
int ABD()
{	
	int i, j, k;
	
	try
	{
		for( int i = 0; i < m_dyn; i++  )
			for( int k = 0; k < deg; k++ )
				for( int j = 0; j < n_dyn; j++ )
					D[ i ][ j ][ k ] = Z_dyn[ i ][ j ][ k ];// D(t)
		for( int i = 0; i < m_dae; i++  )
			for( int k = 0; k < deg; k++ )
			{
				for( int j = 0; j < m_dyn; j++ )
					A[ i ][ j ][ k ] = Z_dae[ i ][ j ][ k ];// A(t)
				for( int j = m_dyn; j < n_dae; j++ )
					B[ i ][ j - m_dyn ][ k ] = Z_dae[ i ][ j ][ k ];// B(t)
			}
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "------------------------" );
				fprintf( pof, "\nMatrices  D,  A,  and  B" );
				fprintf( pof, "\n------------------------" );
				fprint3ddm( pof, m_dyn, n_dyn, deg, D, "\nD = \n" );
				fprint3ddm( pof, m_dae, m_dyn, deg, A, "\nA = \n" );
				fprint3ddm( pof, m_dae, n_dae - m_dyn, deg, B, "\nB = \n" );
			}
			else if( printWhere == 0 )
			{
				printf( cyan );
				printf( "------------------------" );
				printf( "\nMatrices  D,  A,  and  B" );
				printf( "\n------------------------" );
				printf( normal );
				print3ddm( m_dyn, n_dyn, deg, D, "\nD = \n", magenta );
				print3ddm( m_dae, m_dyn, deg, A, "\nA = \n", magenta );
				print3ddm( m_dae, n_dae - m_dyn, deg, B, "\nB = \n", magenta );
			}
		}
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// all other exceptions
	{
		IDException e( "Error when constructing matrices A, B, or D.", 31 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Process data, define variables, and initialize them for the the matrix sequence calculation.
 * Also firsts calculations are computed.
 * 
 */
int matrixSqc( int& r )
{
	int i, j, r_A, r_D, r_i, r_Bi, r_Zt, dim, niter, elem,
		ncolb = 0, ncolZ = 0, TPgrade;
	double fnorm;
	bool cond;

	try
	{
		TPgrade = deg - 1;									// till now was 'deg' the actual nr.
															// of Taylor coefficents!! Thus, 
															// decrement here for being compatible
															// and consistent with the meaning
															
		TPolynMatrix AA( m_dae, m_dyn, TPgrade ),			// AA = matrix of Taylor polyn, i.e. A
			DD( m_dyn, m_dae, TPgrade ),					// DD = first submatrix from D
			DT( m_dae, m_dyn, TPgrade ),					// D^T is the transposed matrix of DD
			Dm( m_dae, m_dyn, TPgrade ),					// inverse of D^T 
			S_A( m_dae, m_dyn, TPgrade ),					// = R after QR decomp. of AA
			U_A( m_dae, m_dae, TPgrade ),					// = Q after QR decomp. of AA
			S_D( m_dae, m_dyn, TPgrade ),					// = R after QR decomp. of D^T
			U_D( m_dae, m_dae, TPgrade ),					// = Q after QR decomp. of D^T
			Ui( m_dae, m_dae, TPgrade ),
			Uim( m_dae, m_dae, TPgrade ),					// U_0^T = U_0^-1 = U_0m
			Gi( m_dae, m_dae, TPgrade ),					// G0 = A*D
			G0( m_dae, m_dae, TPgrade ),					// G0, copy of the original matrix
			Gim( m_dae, m_dae, TPgrade ),					// Gi^-1
			S_Gi( m_dae, m_dae, TPgrade ),					// = R after QR decomposition of Gi
			Vi( m_dae, m_dae, TPgrade ),					// Vi = [Vi1 | Vi2]
			Vi1( m_dae, m_dae, TPgrade ), 
			Vi2( m_dae, m_dae, TPgrade ),
			Ri1( m_dae, m_dae, TPgrade ), 
			Ri2( m_dae, m_dae, TPgrade ),					// decomposing R
			Ri1m( m_dae, m_dae, TPgrade ),					// inverse of Ri1
			b( m_dae, m_dae, TPgrade ),						// b[m_dae][m_dae]
			Si( m_dae, m_dae, TPgrade ), 
			Sim( m_dae, m_dae, TPgrade ),
			B0( m_dae, m_dae, TPgrade ),					// B0 = first submatrix from B
			Bi( m_dae, m_dae, TPgrade ), 
			Bi_old( m_dae, m_dae, TPgrade ),
			Bibar( m_dae, m_dae, TPgrade ),					// Bibar
			Bi22( m_dae, m_dae, TPgrade ),					// Bi22
			Utilde( m_dae, m_dae, TPgrade ), 
			Stilde( m_dae, m_dae, TPgrade ), 
			Vtilde( m_dae, m_dae, TPgrade ),
			Z( m_dae, m_dae, TPgrade ), 
			Zi( m_dae, m_dae, TPgrade ), 
			Ztilde( m_dae, m_dae, TPgrade ),
			SZtilde( m_dae, m_dae, TPgrade ),				// = R after QR decomp. of Ztilde
			UZtilde( m_dae, m_dae, TPgrade ),				// = Q after QR decomp. of Ztilde
			UZtildeT( m_dae, m_dae, TPgrade ),				// UZtilde^T
			Rtilde( m_dae, m_dae, TPgrade ), 
			Rtildem( m_dae, m_dae, TPgrade ),
			mi1( m_dae, m_dae, TPgrade ),					// m01, m11, m21,...
			mi2( m_dae, m_dae, TPgrade ),					// m02, m12, m22,...
			mtilde( m_dae, m_dae, TPgrade ),
			M11( m_dae, m_dae, TPgrade ),
			M12( m_dae, m_dae, TPgrade ),
			Wi( m_dae, m_dae, TPgrade ),
			R( m_dyn, m_dyn, TPgrade ), 					// projector
			Pi( m_dae, m_dae, TPgrade ),
			Pi_old( m_dae, m_dae, TPgrade ), 
			DerivPPi( m_dyn, m_dyn, TPgrade ),
			PPi( m_dyn, m_dyn, TPgrade ), 
			PPi_old( m_dyn, m_dyn, TPgrade ),				// PP_i = PP_{i-1}*D*P_i*Dm
			DPiDm( m_dyn, m_dyn, TPgrade ), 				// D*P_i*Dm
			dDPiDm( m_dyn, m_dyn, TPgrade ), 
			dDPiDm_old( m_dyn, m_dyn, TPgrade ),			// (D*P_i*Dm)'
			aux1( m_dyn, m_dyn, TPgrade ), 					// auxiliary matrices
			aux2( m_dyn, m_dae, TPgrade ), 
			aux3( m_dae, m_dyn, TPgrade ), 
			aux4( m_dae, m_dae, TPgrade ), 
			aux5( m_dae, m_dae, TPgrade ), 
			aux6( m_dae, m_dyn, TPgrade ), 
			aux7( m_dyn, m_dyn, TPgrade ),
			aux8( m_dae, m_dae, TPgrade ),
			I( m_dae, m_dae, TPgrade ),
			Tlim( m_dae, m_dae, TPgrade ),					// Tli * matrix
			mTui( m_dae, m_dae, TPgrade ),					// matrix * Tui
			Uim_new( m_dae, m_dae, TPgrade ),
			Si_new( m_dae, m_dae, TPgrade ), 
			Sim_new( m_dae, m_dae, TPgrade ),
			Wihat( m_dae, m_dae, TPgrade ), 
			Wihat1( m_dae, m_dae, TPgrade ), 
			Wihat2( m_dae, m_dae, TPgrade );
		TPolynMatrix *Qi = new TPolynMatrix[ m_dae ];		// array of Qi's
		for( i = 0; i < m_dae; i++ )
			Qi[ i ] = TPolynMatrix ( m_dae, m_dae, TPgrade );
		int *PI_A = new int[ m_dyn ];
		int *PI_D = new int[ m_dyn ];						// defined for D^T!
		int *PI_G0 = new int[ m_dae ];
		int *PItilde = new int[ m_dae ];					// for Bi22
		int *PIZtilde = new int[ m_dae ];					// for Ztilde
		double *fn = new double[ TPgrade + 1 ];				// array of Frobenius norms, for all
															// Taylor coefficients!
		if( (PI_A == 0)  ||  (PI_A == NULL)  ||  (PI_D == 0)  ||  (PI_D == NULL)  ||
			(PI_G0 == 0)  ||  (PI_G0 == NULL)  ||  (PItilde == 0)  ||  (PItilde == NULL)  ||
			(PIZtilde == 0)  ||  (PIZtilde == NULL)  ||  (fn == 0)  ||  (fn == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );

		//
		// I M P O R T A N T   I N I T I A L I Z A T I O N S
		//
		//printf( "\nInitializations for the matrix sequence...\n" );
		for( i = 0; i < m_dae; i++ )
		{
			for( j = 0; j < m_dyn; j++ )
				AA( i, j ).setCoeffs( A[ i ][ j ] );		// initialize AA
			for( j = 0; j < m_dae; j++ )
			{
				Bi( i, j ).setCoeffs( B[ i ][ j ] );		// initialize Bi (B0 = B)
				B0( i, j ).setCoeffs( B[ i ][ j ] );		// initialize B0
			}
		}
		for( i = 0; i < m_dyn; i++ )
			for( j = 0; j < m_dae; j++ )
				DD( i, j ).setCoeffs( D[ i ][ j ] );		// initialize DD
		Gi.mmCaABbC( 1.0, 1.0, AA, DD );					// G0 = A * D
		//Gi.printtpm( "\nG0 = A * D\n", black, ioeps );
		//
		// Q R   F A C T O R I Z A T I O N   F O R   M A T R I C E S   AA ,  DD   A N D   Gi
		// A N D   I N I T I A L   T E S T S
		//
		if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\n---------------------------------------------------" );
				fprintf( pof, "\nInitializations for the matrix sequence computation" );
				fprintf( pof, "\n---------------------------------------------------" );
				AA.fprinttpm( pof, "\nA = \n", ioeps );
				DD.fprinttpm( pof, "\nD = \n", ioeps );
				Gi.fprinttpm( pof, "\nG0 = A * D\n", ioeps );
				fprintf( pof, "\n---------------------------------------------------------------" );
				fprintf( pof, "\nHouseholder QR Factorization with Column Pivoting for matrix AA" );
				fprintf( pof, "\n       A = Q * R * P    <=>    A = U_A * S_A * PI_A^T" );
				fprintf( pof, "\n---------------------------------------------------------------\n" );
			}
			else if( printWhere == 0 )
			{
				printf( cyan );
				printf( "\n---------------------------------------------------" );
				printf( "\nInitializations for the matrix sequence computation" );
				printf( "\n---------------------------------------------------" );
				AA.printtpm( "\nA = \n", black, ioeps );
				DD.printtpm( "\nD = \n", black, ioeps );
				Gi.printtpm( "\nG0 = A * D\n", black, ioeps );
				printf( cyan );
				printf( "\n---------------------------------------------------------------" );
				printf( "\nHouseholder QR Factorization with Column Pivoting for matrix AA" );
				printf( "\n       A = Q * R * P    <=>    A = U_A * S_A * PI_A^T" );
				printf( "\n---------------------------------------------------------------\n" );
				printf( normal );
			}
		}
		//
		// 1st   H O U S E H O L D E R   Q R   F A C T O R I Z A T I O N
		//
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintf( pof, "Applying Householder to A...\n" );
			else if( printWhere == 0 )
				printf( "Applying Householder to A...\n" );
		}
/*(HH1)*/hqrfcp( eps, AA, PI_A, U_A, S_A, r_A, 0, printWhat,// QR fact. for AA, non-permuted
			printWhere, pof );								// S_A is an extra out-parameter
		//constructR( AA, PI_A, S_A, 0 );
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintf( pof, "\n\nRank(A) = %d\n", r_A );
			else if( printWhere == 0 )
				printf( "\n\nRank(A) = %d\n", r_A );
		}
		if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				//S_A.fprinttpm( pof, "\nS_A = \n", ioeps );
				//fprintiv( pof, m_dyn, PI_A, "\nPI_A = \n" );
				fprintf( pof, "\n---------------------------------------------------------------" );
				fprintf( pof, "\nHouseholder QR Factorization with Column Pivoting for matrix DD" );
				fprintf( pof, "\n       D = Q * R * P    <=>    D = U_D * S_D * PI_D^T" );
				fprintf( pof, "\n       (D is first transposed!)" );
				fprintf( pof, "\n---------------------------------------------------------------\n" );
			}
			else if( printWhere == 0 )
			{
				//S_A.printtpm( "\nS_A = \n", black, ioeps );
				//printiv( m_dyn, PI_A, "\nPI_A = \n", black );
				printf( cyan );
				printf( "\n---------------------------------------------------------------" );
				printf( "\nHouseholder QR Factorization with Column Pivoting for matrix DD" );
				printf( "\n       D = Q * R * P    <=>    D = U_D * S_D * PI_D^T " );
				printf( "\n       (D is first transposed!)" );
				printf( "\n---------------------------------------------------------------\n" );
				printf( normal );
			}
		}
		DT = DD.transpm();									// Transpose matrix DD to 
															// obtain D^T[m_dae][m_dyn]
		if( printWhat > 21 )
		{
			if( printWhere == 1 )
				DT.fprinttpm( pof, "\nD^T = \n", ioeps );
			else if( printWhere == 0 )
				DT.printtpm( "\nD^T = \n", black, ioeps );
		}
		//
		// 2nd   H O U S E H O L D E R   Q R   F A C T O R I Z A T I O N
		//
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintf( pof, "\nApplying Householder to D^T...\n" );
			else if( printWhere == 0 )
				printf( "\nApplying Householder to D^T...\n" );
		}
/*(HH2)*/hqrfcp( eps, DT, PI_D, U_D, S_D, r_D, 0, printWhat,// QR fact. for DT, non-permuted
			printWhere, pof );								// S_D is an extra out-parameter, 
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\n\nRank(D^T) = %d\n", r_D );
				//S_D.fprinttpm( pof, "\nS_D = \n", ioeps );
				//fprintiv( pof, m_dyn, PI_D, "\nPI_D = \n" );
			}
			else if( printWhere == 0 )
			{
				printf( "\n\nRank(D^T) = %d\n", r_D );
				//S_D.printtpm( "\nS_D = \n", black, ioeps );
				//printiv( m_dyn, PI_D, "\nPI_D = \n", black );
			}
		}	
		if( r_A != r_D )									// Testing whether properly 
															// stated leading term or not
			throw IDException( "Not properly stated leading term:  rank(AA) <> rank(DD^T).", 6 );
		if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				fprintf( pof, "\n---------------------------------------------------------------" );
				fprintf( pof, "\nHouseholder QR Factorization with Column Pivoting for matrix G0" );
				fprintf( pof, "\n     G0 = Q * R * P     <=>     G0 = U_G0 * S_G0 * PI_G0^T" );
				fprintf( pof, "\n---------------------------------------------------------------\n" );
			}
			else if( printWhere == 0 )
			{
				printf( cyan );
				printf( "\n---------------------------------------------------------------" );
				printf( "\nHouseholder QR Factorization with Column Pivoting for matrix G0" );
				printf( "\n     G0 = Q * R * P     <=>     G0 = U_G0 * S_G0 * PI_G0^T" );
				printf( "\n---------------------------------------------------------------\n" );
				printf( normal );
			}
		}	
		//
		// 3rd   H O U S E H O L D E R   Q R   F A C T O R I Z A T I O N
		//
		G0 = Gi;											// save original matrix
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintf( pof, "\nApplying Householder to G0=A*D...\n" );
			else if( printWhere == 0 )
				printf( "\nApplying Householder to G0=A*D...\n" );
		}
/*(HH3)*/hqrfcp( eps, Gi, PI_G0, Ui, S_Gi, r_i, 0, printWhat,// QR fact. for G0, non-permuted
			printWhere, pof );								// S_Gi is an extra out-parameter, 

		Gi = G0;											// reload original matrix
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintf( pof, "\n\nRank(G0=A*D) = %d\n", r_i );
			else if( printWhere == 0 )
				printf( "\n\nRank(G0=A*D) = %d\n", r_i );
		}
		if( r_A != r_i )									// Testing whether properly 
															// stated leading term or not
			throw IDException( "Not properly stated leading term:  rank(A) <> rank(G0).", 7 );
		r = r_A;											// Computed rank until now!!
		if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				S_Gi.fprinttpm( pof, "\nS_G0 = \n", ioeps );
				Ui.fprinttpm( pof, "\nUi = \n", ioeps );
				fprintiv( pof, m_dae, PI_G0, "\nPI_G0 = \n" );
				fprintf( pof, "\n------------------" );
				fprintf( pof, "\nOther calculations" );
				fprintf( pof, "\n------------------" );
			}
			else if( printWhere == 0 )
			{
				S_Gi.printtpm( "\nS_G0 = \n", black, ioeps );
				printiv( m_dae, PI_G0, "\nPI_G0 = \n", black );
				printf( cyan );
				printf( "\n------------------" );
				printf( "\nOther calculations" );
				printf( "\n------------------" );
				printf( normal );
			}
		}
		if( r == m_dae )									// we are done! "this" is the rank!
		{
			if( printWhat > 0 )
			{
				if( printWhere == 1 )
				{
					fprintf( pof, "\n\nThe rank is = %d\n", r );
					fprintf( pof, "\n\nStopping calculations...\n" );
				}
				else if( printWhere == 0 )
				{
					printf( "\n\nThe rank is = %d\n", r );
					printf( "\n\nStopping calculations...\n" );
				}
			}
		}
		else
		{
			//
			// I N I T I A L I Z A T I O N S   F O R   T H E   M A T R I X   S E Q U E N C E
			//
/*(S.4)*/
			//
			// C O M P U T E   Dm
			//
			if( printWhat > 1 )
			{
				if( printWhere == 1 )
					fprintf( pof, "\nComputing Dm..." );
				else if( printWhere == 0 )
					printf( "\nComputing Dm..." );
			}
			computeDm( r, S_A, S_D, PI_A, PI_D, U_D, Dm );	// compute Dm
			if( printWhat > 21 )
				if( printWhere == 1 )
					Dm.fprinttpm( pof, "\nDm = \n", ioeps );
				else if( printWhere == 0 )
					Dm.printtpm( "\nDm = \n", black, ioeps );
			//
			// T E S T:  D = D*Dm*D ???
			//
			aux1.mmCaABbC( 1.0, 1.0, DD, Dm );				// aux1 = D*Dm
			R = aux1;										// projector R
			aux2.mmCaABbC( 1.0, 1.0, aux1, DD );			// aux2 = D*Dm*D
			if( printWhat > 1 )
			{
				if( printWhere == 1 )
					R.fprinttpm( pof, "\nR = D*Dm =\n", ioeps );
				else if( printWhere == 0 )
					R.printtpm( "\nR = D*Dm =\n", black, ioeps );
			}
			if( printWhat > 21 )
			{
				if( printWhere == 1 )
				{
					aux2.fprinttpm( pof, "\nD*Dm*D = \n", ioeps );
					fprintf( pof, "\nT E S T:   D == D*Dm*D ???  (eps=%.16lg)\n", eps );
				}
				else if( printWhere == 0 )
				{
					aux2.printtpm( "\nD*Dm*D = \n", black, ioeps );
					printf( "\nT E S T:   D == D*Dm*D ???  (eps=%.16lg)\n", eps );
				}
			}
			aux2 -= DD;								// D*Dm*D - D  (should be zero...)
			if( printWhat == 1 )
			{
				fnorm = aux2.fnorm();				// Frobenius norm, first Taylor coeff.
				if( fnorm > eps )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\nWARNING  -->  |D*Dm*D - D| > eps (for eps=%.16lg)\n", eps );
						fprintf( pof, "\nFrobenius norm of D*Dm*D - D = %.16lg > eps", fnorm );
						fprintf( pof, "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
					}
					else if( printWhere == 0 )
					{
						printf( red );
						printf( "\nWARNING  -->  |D*Dm*D - D| > eps (for eps=%.16lg)\n", eps );
						printf( "\nFrobenius norm of D*Dm*D - D = %.16lg > eps", fnorm );
						printf( "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
						printf( normal );
					}
				}
			}
			else if( printWhat > 1 )
			{
				aux2.tcfnorm( TPgrade + 1, fn );	// Frobenius norms, all Taylor coeff.
				cond = vleqeps( TPgrade + 1, fn, eps, elem );
				if( !cond )							// at least one element is greater
				{									// than the epsilon
					if( printWhere == 1 )
					{
						fprintf( pof, "\nWARNING  -->  |D*Dm*D - D| > eps (for eps=%.16lg)\n", eps );
						fprintf( pof, "\nAt least for one Taylor coefficient holds: Frobenius norm of D*Dm*D - D > eps" );
						fprintdv( pof, TPgrade + 1, fn, "\nFrobenius norms =\n" );
						if( printWhat > 21 )
							aux2.fprinttpm( pof, "\nD*Dm*D - D = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						printf( red );
						printf( "\nWARNING  -->  |D*Dm*D - D| > eps (for eps=%.16lg)\n", eps );
						printf( "\nAt least for one Taylor coefficient holds: the Frobenius norm of D*Dm*D - D > eps" );
						printf( normal );
						printdv( TPgrade + 1, fn, "\nFrobenius norms =\n", red );
						if( printWhat > 21 )
							aux2.printtpm( "\nD*Dm*D - D = \n", red, ioeps );
					}
				}
			}
			//
			// T E S T:  Dm = Dm*D*Dm ???
			//
			aux3.mmCaABbC( 1.0, 1.0, Dm, aux1 );			// aux3 = Dm*D*Dm
			if( printWhat > 21 )
			{
				if( printWhere == 1 )
				{
					aux3.fprinttpm( pof, "\nDm*D*Dm = \n", ioeps );
					fprintf( pof, "\nT E S T:   Dm == Dm*D*Dm ???  (eps=%.16lg)\n", eps );
				}
				else if( printWhere == 0 )
				{
					aux3.printtpm( "\nDm*D*Dm = \n", black, ioeps );
					printf( "\nT E S T:   Dm == Dm*D*Dm ???  (eps=%.16lg)\n", eps );
				}
			}
			aux3 -= Dm;								// Dm*D*Dm - Dm  (should be zero...)
			if( printWhat == 1 )
			{
				fnorm = aux3.fnorm();				// Frobenius norm, first Taylor coeff.
				if( fnorm > eps )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\nWARNING  -->  |Dm*D*Dm - Dm| > eps (for eps=%.16lg)\n", eps );
						fprintf( pof, "\nFrobenius norm of Dm*D*Dm - Dm = %.16lg > eps", fnorm );
						fprintf( pof, "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
					}
					else if( printWhere == 0 )
					{
						printf( red );
						printf( "\nWARNING  -->  |Dm*D*Dm - Dm| > eps (for eps=%.16lg)\n", eps );
						printf( "\nFrobenius norm of Dm*D*Dm - Dm = %.16lg > eps", fnorm );
						printf( "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
						printf( normal );
					}
				}
			}
			else if( printWhat > 1 )
			{
				aux3.tcfnorm( TPgrade + 1, fn );	// Frobenius norms, all Taylor coeff.
				cond = vleqeps( TPgrade + 1, fn, eps, elem );
				if( !cond )							// at least one element is greater
				{									// than the epsilon
					if( printWhere == 1 )
					{
						fprintf( pof, "\nWARNING  -->  |Dm*D*Dm - Dm| > eps (for eps=%.16lg)\n", eps );
						fprintf( pof, "\nAt least for one Taylor coefficient holds: Frobenius norm of Dm*D*Dm - Dm > eps" );
						fprintdv( pof, TPgrade + 1, fn, "\nFrobenius norms =\n" );
						if( printWhat > 21 )
							aux3.fprinttpm( pof, "\nDm*D*Dm - Dm = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						printf( red );
						printf( "\nWARNING  -->  |Dm*D*Dm - Dm| > eps (for eps=%.16lg)\n", eps );
						printf( "\nAt least for one Taylor coefficient holds: the Frobenius norm of Dm*D*Dm - Dm > eps" );
						printf( normal );
						printdv( TPgrade + 1, fn, "\nFrobenius norms =\n", red );
						if( printWhat > 21 )
							aux3.printtpm( "\nDm*D*Dm - Dm = \n", red, ioeps );
					}
				}
			}
			//
			// D E C O M P O S E   S_Gi
			//
			Ri1.redim( r, r, TPgrade );
			Ri1m.redim( r, r, TPgrade );
			Ri1m.set2zero();								// inverse of Ri1
			Ri2.redim( r, m_dae - r, TPgrade );
			if( printWhat > 1 )
			{
				if( printWhere == 1 )
					fprintf( pof, "\nComputing R01m, -R01m*R02, V0, U0, ..." );
				else if( printWhere == 0 )
					printf( "\nComputing R01m, -R01m*R02, V0, U0, ..." );
			}
			decompS( r, S_Gi, PI_G0, Ri1, Ri2 );
			if( r > 0 )										// Ri1 exists
			{
				//
				// T R I A N G U L A R   I N V E R S E   of   Ri1m
				// (Ri1m is upper triangular; a especial inverse calculation procedure is used)
				//
				Ri1.trinvm( Ri1m );
				aux5.redim( r, r, TPgrade );
				aux5.mmCaABbC( 1.0, 1.0, Ri1, Ri1m );		// Ri1*Ri1m 
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						//Ri1m.fprinttpm( pof, "\nR01m after trinvm = \n", ioeps );
						fprintf( pof, "\nT E S T:   Ri1*Ri1m == I ???\n" );
					}
					else if( printWhere == 0 )
					{
						//Ri1m.printtpm( "\nR01m after trinvm = \n", black, ioeps );
						printf( "\nT E S T:   Ri1*Ri1m == I ???\n" );
					}
				}
				if( !aux5.isId( eps ) )						// R01*R01m == I ?
					if( printWhat > 0 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\nWARNING  -->  R01*R01m <> I\n" );
						else if( printWhere == 0 )
						{
							printf( red );					
							printf( "\nWARNING  -->  R01*R01m <> I\n" );
							printf( normal );
						}
					}
				Ri1m = -Ri1m;
				if( r < m_dae )								// there exists also Ri2
				{
					aux5.setncols( m_dae - r );
					aux5.mmCaABbC( 1.0, 1.0, Ri1m, Ri2 );	// aux5=-Ri1m*Ri2
					//aux5.fprinttpm( pof, "\n-Ri1m*Ri2 = \n", ioeps );
					//
					// C O M P U T E   V0
					//
					computeV( PI_G0, 0, aux5, Vi );			// rows permut.
				}
				else
				{
					Vi.set2Id();							// Vi is set to the identity
					Vi.rpermutem( PI_G0 );					// row permutation
				}
			}
			else											// Ri1 does not exist but Ri2
			{
				Vi.set2Id();								// Vi is set to the identity
				Vi.rpermutem( PI_G0 );						// row permutation
			}
			Uim = Ui.transpm();								// transpose U_G0[m_dae][m_dae]
															// U0T = U0m
			if( printWhat > 21 )
			{
				if( printWhere == 1 )
				{
					Ui.fprinttpm( pof, "\nU0 = \n", ioeps );
					Uim.fprinttpm( pof, "\nU0m = U0T = \n", ioeps );
				}
				else if( printWhere == 0 )
				{
					Ui.printtpm( "\nU0 = \n", black, ioeps );
					Uim.printtpm( "\nU0m = U0T = \n", black, ioeps );
				}
			}
			if( printWhat > 1 )
			{
				if( printWhere == 1 )
					fprintf( pof, "\nComputing G0m..." );
				else if( printWhere == 0 )
					printf( "\nComputing G0m..." );
			}
			if( r > 0 )										// Ri1 and Ri1m exist
			{
				//
				// C O M P U T E   G0m
				//
				mi1.redim( m_dae - r, r, TPgrade );
				mi2.redim( r, m_dae - r, TPgrade );
				Ri1m = -Ri1m;
				mi1.set2zero();								// mi1 = 0, free parameter
				mi2.set2zero();								// mi2 = 0, free parameter
				computeGm( mi1, mi2, Vi, Ri1, Ri1m, Uim, Gim );
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
						Gim.fprinttpm( pof, "\nG0m = \n", ioeps );
					else if( printWhere == 0 )
						Gim.printtpm( "\nG0m = \n", black, ioeps );
				}
				Si.redim( r, r, TPgrade );
				Si = Ri1;									// S0 = R_01
			}
			else
			{
				Gim.mmCaABbC( 1.0, 1.0, Vi, Uim );			// Gim=PI_G0*Vi*Uim
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
						Gim.fprinttpm( pof, "\nG0m = \n", ioeps );
					else if( printWhere == 0 )
						Gim.printtpm( "\nG0m = \n", black, ioeps );
				}
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nComputing V01 and V02..." );
					else if( printWhere == 0 )
						printf( "\nComputing V01 and V02..." );
				}
/*??*/			//Si = ???
			}
			Vi1.setncols( r );
			Vi2.setncols( m_dae - r );
			decompV( r, Vi, Vi1, Vi2 );						// Vi = [ Vi1 | Vi2 ]
			ncolb += m_dae - r;								// update nr. of columns from b
			b.setncols( ncolb );
			filloutM( ncolb - m_dae + r, m_dae - r, Vi2, b );// b = V02
			
			dim = r;										// update global counters
			niter = -1;
			if( printWhat > 1 )
			{
				if( printWhere == 1 )
					fprintf( pof, "\n\nniter = %d,  dim = %d\n", niter, dim );
				else if( printWhere == 0 )
					printf( "\n\nniter = %d,  dim = %d\n", niter, dim );
			}
			PPi.set2Id();									// PPi = I
			Bi_old = B0;
			//
			//
			//  M   A   T   R   I   X       S   E   Q   U   E   N   C   E       L   O   O   P
			//
			//
/*(S.4)*/
			if( printWhat > 1 )
			{
				if( printWhere == 1 )
				{
					fprintf( pof, "\n---------------------------------------" );
					fprintf( pof, "\nM A T R I X   S E Q U E N C E   L O O P" );
					fprintf( pof, "\n---------------------------------------" );
				}
				else if( printWhere == 0 )
				{
					printf( cyan );
					printf( "\n---------------------------------------" );
					printf( "\nM A T R I X   S E Q U E N C E   L O O P" );
					printf( "\n---------------------------------------" );
					printf( normal );
				}
			}
			while( (dim < m_dae)  &&  (niter <= m_dae) )
			{
				//
				// I N I T I A L   C O M P U T A T I O N S
				//
				niter++;									// increment global counter
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\n\ni = %d,  dim = %d", niter, dim );
						fprintf( pof, "\n----------------\n" );
					}
					else if( printWhere == 0 )
					{
						printf( "\n\ni = %d,  dim = %d", niter, dim );
						printf( "\n----------------\n" );
					}
				}
				if( niter > 0 )
					Qi[ niter - 1 ] = Qi[ niter ];
				Pi_old = Pi;
				pwq( Gi, Gim, Pi, Wi, Qi[ niter ] );		// calculate Pi, Wi, Qi
				if( printWhat > 2 )
				{
					if( printWhere == 1 )
					{
						Gim.fprinttpm( pof, "\nGim = \n", ioeps );
						Gi.fprinttpm( pof, "\nGi = \n", ioeps );
						Pi.fprinttpm( pof, "\nPi = \n", ioeps );
						Qi[ niter ].fprinttpm( pof, "\nQi = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						Gim.printtpm( "\nGim = \n", black, ioeps );
						Gi.printtpm( "\nGi = \n", black, ioeps );
						Pi.printtpm( "\nPi = \n", black, ioeps );
						Qi[ niter ].printtpm( "\nQi = \n", black, ioeps );
					}
				}
				//
				// T E S T:  Q_i+1*Qi == 0 ??
				//
				aux4.redim( m_dae, m_dae, TPgrade );
				if( niter > 0 )
				{
					if( printWhat > 21 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\nT E S T:   Qi*Qj == 0 ???  (eps=%.16lg)\n", eps );
						else if( printWhere == 0 )
							printf( "\nT E S T:   Qi*Qj == 0 ???  (eps=%.16lg)\n", eps );
					}
					for( j = 0; j < niter; j++ )
					{
						aux4.mmCaABbC( 1.0, 1.0, Qi[ niter ], Qi[ j ] );// Qi*Qj
						if( printWhat > 21 )
						{
							if( printWhere == 1 )
							{
								fprintf( pof, "\nFor j = %d:", j );
								aux4.fprinttpm( pof, "\nQi*Qj = \n", ioeps );
							}
							else if( printWhere == 0 )
							{
								printf( "\nFor j = %d:", j );
								aux4.printtpm( "\nQi*Qj = \n", black, ioeps );
							}
						}
						if( printWhat == 1 )
						{
							fnorm = aux4.fnorm();				// Frobenius norm, first Taylor coeff.
							if( fnorm > eps )
							{
								if( printWhere == 1 )
								{
									fprintf( pof, "\nWARNING  -->  |Qi*Qj| > eps (for eps=%.16lg)\n", eps );
									fprintf( pof, "\nFrobenius norm of Qi*Qj = %.16lg > eps", fnorm );
									fprintf( pof, "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
								}
								else if( printWhere == 0 )
								{
									printf( red );
									printf( "\nWARNING  -->  |Qi*Qj| > eps (for eps=%.16lg)\n", eps );
									printf( "\nFrobenius norm of Qi*Qj = %.16lg > eps", fnorm );
									printf( "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
									printf( normal );
								}
							}
						}
						else if( printWhat > 1 )
						{
							aux4.tcfnorm( TPgrade + 1, fn );	// Frobenius norms, all Taylor coeff.
							cond = vleqeps( TPgrade + 1, fn, eps, elem );
							if( !cond )							// at least one element is greater
							{									// than the epsilon
								if( printWhere == 1 )
								{
									fprintf( pof, "\nWARNING  -->  |Qi*Qj| > eps (for eps=%.16lg)\n", eps );
									fprintf( pof, "\nAt least for one Taylor coefficient holds: Frobenius norm of Qi*Qj > eps" );
									fprintdv( pof, TPgrade + 1, fn, "\nFrobenius norms =\n" );
									if( printWhat > 21 )
									{
										fprintf( pof, "\nFor j = %d:", j );
										aux4.fprinttpm( pof, "\nQi*Qj = \n", ioeps );
									}
								}
								else if( printWhere == 0 )
								{
									printf( red );
									printf( "\nWARNING  -->  |Qi*Qj| > eps (for eps=%.16lg)\n", eps );
									printf( "\nAt least for one Taylor coefficient holds: the Frobenius norm of Qi*Qj > eps" );
									printf( normal );
									printdv( TPgrade + 1, fn, "\nFrobenius norms =\n", red );
									if( printWhat > 21 )
									{
										printf( "\nFor j = %d:", j );
										aux4.printtpm( "\nQi*Qj = \n", red, ioeps );
									}
								}
							}
						}
					}
				}
				//
				// T E S T:  Qi^2 == Qi ???
				//
				aux4.mmCaABbC( 1.0, 1.0, Qi[ niter ], Qi[ niter ] );// Qi^2 * Qi
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\nT E S T:   Qi^2 == Qi ???  (eps=%.16lg)\n", eps );
						//aux4.fprinttpm( pof, "\nQi^2 = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						printf( "\nT E S T:   Qi^2 == Qi ???  (eps=%.16lg)\n", eps );
						//aux4.printtpm( "\nQi^2 = \n", black, ioeps );
					}
				}
				/*if( Qi.mcompare( aux4, eps ) )			// old test
					printf( "\nQi^2 = Qi\n" );
				else printf( "\nQi^2 <> Qi\n" );*/
				aux4 -= Qi[ niter ];						// Qi^2 - Qi
				if( printWhat == 1 )
				{
					fnorm = aux4.fnorm();				// Frobenius norm, first Taylor coeff.
					if( fnorm > eps )
					{
						if( printWhere == 1 )
						{
							fprintf( pof, "\nWARNING  -->  |Qi^2 - Qi| > eps (for eps=%.16lg)\n", eps );
							fprintf( pof, "\nThe calculation of the projectors might have problems..." );
							fprintf( pof, "\nFrobenius norm of Qi^2 - Qi = %.16lg > eps", fnorm );
							fprintf( pof, "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
						}
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\nWARNING  -->  |Qi^2 - Qi| > eps (for eps=%.16lg)\n", eps );
							printf( "\nThe calculation of the projectors might have problems..." );
							printf( "\nFrobenius norm of Qi^2 - Qi = %.16lg > eps", fnorm );
							printf( "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
							printf( normal );
						}
					}
				}
				else if( printWhat > 1 )
				{
					aux4.tcfnorm( TPgrade + 1, fn );	// Frobenius norms, all Taylor coeff.
					cond = vleqeps( TPgrade + 1, fn, eps, elem );
					if( !cond )							// at least one element is greater
					{									// than the epsilon
						if( printWhere == 1 )
						{
							fprintf( pof, "\nWARNING  -->  |Qi^2 - Qi| > eps (for eps=%.16lg)\n", eps );
							fprintf( pof, "\nThe calculation of the projectors might have problems..." );
							fprintf( pof, "\nAt least for one Taylor coefficient holds: Frobenius norm of Qi^2 - Qi > eps" );
							fprintdv( pof, TPgrade + 1, fn, "\nFrobenius norms =\n" );
							if( printWhat > 21 )
								aux4.fprinttpm( pof, "\nQi^2 - Qi = \n", ioeps );
						}
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\nWARNING  -->  |Qi^2 - Qi| > eps (for eps=%.16lg)\n", eps );
							printf( "\nThe calculation of the projectors might have problems..." );
							printf( "\nAt least for one Taylor coefficient holds: the Frobenius norm of Qi^2 - Qi > eps" );
							printf( normal );
							printdv( TPgrade + 1, fn, "\nFrobenius norms =\n", red );
							if( printWhat > 21 )
								aux4.printtpm( "\nQi^2 - Qi = \n", red, ioeps );
						}
					}
				}
				//
				// T E S T:  Gi*Qi == 0 ??
				//
				aux4.mmCaABbC( 1.0, 1.0, Gi, Qi[ niter ] );	// Gi * Qi
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nT E S T:   Gi*Qi == 0 ???  (i=%d, eps=%.16lg)\n", niter, eps );
					else if( printWhere == 0 )
						printf( "\nT E S T:   Gi*Qi == 0 ???  (i=%d, eps=%.16lg)\n", niter, eps );
				}
				if( printWhat == 1 )
				{
					fnorm = aux4.fnorm();				// Frobenius norm, first Taylor coeff.
					if( fnorm > eps )
					{
						if( printWhere == 1 )
						{
							fprintf( pof, "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
							fprintf( pof, "\nFrobenius norm of Gi*Qi = %.16lg > eps", fnorm );
							fprintf( pof, "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
						}
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
							printf( "\nFrobenius norm of Gi*Qi = %.16lg > eps", fnorm );
							printf( "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
							printf( normal );
						}
					}
				}
				else if( printWhat > 1 )
				{
					aux4.tcfnorm( TPgrade + 1, fn );	// Frobenius norms, all Taylor coeff.
					cond = vleqeps( TPgrade + 1, fn, eps, elem );
					if( !cond )							// at least one element is greater
					{									// than the epsilon
						if( printWhere == 1 )
						{
							fprintf( pof, "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
							fprintf( pof, "\nAt least for one Taylor coefficient holds: the Frobenius norm of Gi*Qi > eps" );
							fprintdv( pof, TPgrade + 1, fn, "\nFrobenius norms =\n" );
							if( printWhat > 21 )
								aux4.fprinttpm( pof, "\nGi*Qi = \n", ioeps );
						}
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
							printf( "\nAt least for one Taylor coefficient holds: the Frobenius norm of Gi*Qi > eps" );
							printf( normal );
							printdv( TPgrade + 1, fn, "\nFrobenius norms =\n", red );
							if( printWhat > 21 )
								aux4.printtpm( "\nGi*Qi = \n", red, ioeps );
						}
					}
				}
				if( niter == 0 )							// only the first time
				{
					//
					// T E S T:  P0 == Dm*D ??
					//
					aux4.mmCaABbC( 1.0, 1.0, Dm, DD );		// Dm*D
					aux4 -= Pi;								// Dm*D - P0  (should be zero...)
					if( printWhat > 21 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\nT E S T:   P0 == Dm*D ???  (eps=%.16lg)\n", eps );
						else if( printWhere == 0 )
							printf( "\nT E S T:   P0 == Dm*D ???  (eps=%.16lg)\n", eps );
					}
					if( printWhat == 1 )
					{
						fnorm = aux4.fnorm();				// Frobenius norm, first Taylor coeff.
						if( fnorm > eps )
						{
							if( printWhere == 1 )
							{
								fprintf( pof, "\nWARNING  -->  |Dm*D - P0| > eps (for eps=%.16lg)\n", eps );
								fprintf( pof, "\nFrobenius norm of Dm*D - P0 = %.16lg > eps", fnorm );
								fprintf( pof, "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
							}
							else if( printWhere == 0 )
							{
								printf( red );
								printf( "\nWARNING  -->  |Dm*D - P0| > eps (for eps=%.16lg)\n", eps );
								printf( "\nFrobenius norm of Dm*D - P0 = %.16lg > eps", fnorm );
								printf( "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
								printf( normal );
							}
						}
					}
					else if( printWhat > 1 )
					{
						aux4.tcfnorm( TPgrade + 1, fn );	// Frobenius norms, all Taylor coeff.
						cond = vleqeps( TPgrade + 1, fn, eps, elem );
						if( !cond )							// at least one element is greater
						{									// than the epsilon
							if( printWhere == 1 )
							{
								fprintf( pof, "\nWARNING  -->  |Dm*D - P0| > eps (for eps=%.16lg)\n", eps );
								fprintf( pof, "\nAt least for one Taylor coefficient holds: Frobenius norm of Dm*D - P0 > eps" );
								fprintdv( pof, TPgrade + 1, fn, "\nFrobenius norms =\n" );
								if( printWhat > 21 )
									aux4.fprinttpm( pof, "\nDm*D - P0 = \n", ioeps );
							}
							else if( printWhere == 0 )
							{
								printf( red );
								printf( "\nWARNING  -->  |Dm*D - P0| > eps (for eps=%.16lg)\n", eps );
								printf( "\nAt least for one Taylor coefficient holds: the Frobenius norm of Dm*D - P0 > eps" );
								printf( normal );
								printdv( TPgrade + 1, fn, "\nFrobenius norms =\n", red );
								if( printWhat > 21 )
									aux4.printtpm( "\nDm*D - P0 = \n", red, ioeps );
							}
						}
					}
				}
				//
				// C O M P U T E   D*Pi*Dm  and  PPi
				//
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nComputing D*Pi*Dm and updating PPi..." );
					else if( printWhere == 0 )
						printf( "\nComputing D*Pi*Dm and updating PPi..." );
				}
				aux2.mmCaABbC( 1.0, 1.0, DD, Pi );			// D[m_dyn][m_dae]*Pi[m_dae][m_dae]
				//aux2.printtpm( "\nD*Pi = \n", black, ioeps );
				DPiDm.mmCaABbC( 1.0, 1.0, aux2, Dm );		// D*Pi * Dm[m_dae][m_dyn]
				aux7.mmCaABbC( 1.0, 1.0, PPi, DPiDm );		// PPi_old[m_dyn][m_dyn] * DPiDm
				PPi_old = PPi;
				PPi = aux7;									// update PPi
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
					{
						DPiDm.fprinttpm( pof, "\nD*Pi*Dm = \n", ioeps );
						PPi.fprinttpm( pof, "\nupdated PPi = PP_i-1*D*Pi*Dm = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						DPiDm.printtpm( "\nD*Pi*Dm = \n", black, ioeps );
						PPi.printtpm( "\nupdated PPi = PP_i-1*D*Pi*Dm = \n", black, ioeps );
					}
				}
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nComputing Bibar, Bibar22, ..." );
					else if( printWhere == 0 )
						printf( "\nComputing Bibar, Bibar22, ..." );
				}
				//
				// C O M P U T E   Bibar
				//
				aux4.redim( m_dae - dim, B0.ncols(), TPgrade );
				aux4.mmCasABbC( m_dae - dim, 1.0, 1.0, Uim, B0 );// special A
				//printf( "\nOnly last -%d- rows from Uim are of interest:", m_dae - dim );
				//aux4.printtpm( "\nUim*B0 = \n", black, ioeps );
				Bi22.redim( m_dae - dim, m_dae - dim, TPgrade );
				Bi22.mmCaABbC( 1.0, 1.0, aux4, Vi2 );
				//Bi22.printtpm( "\nBi22 = \n", black, ioeps );
				//
				// 4th   H O U S E H O L D E R   Q R   F A C T O R I Z A T I O N
				//
				Stilde.redim( m_dae - dim, m_dae - dim, TPgrade );
				Utilde.redim( m_dae - dim, m_dae - dim, TPgrade );
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
					{
						//Bi22.fprinttpm( pof, "\nBibar22 = right, lower part of Uim*B0*Vi2 = \n", ioeps );
						fprintf( pof, "\nApplying Householder to Bibar22...\n" );
					}
					else if( printWhere == 0 )
					{
						//Bi22.printtpm( "\nBibar22 = right, lower part of Uim*B0*Vi2 = \n", black, ioeps );
						printf( "\nApplying Householder to Bibar22...\n" );
					}
				}
/*(HH4)*/		hqrfcp( eps, Bi22, PItilde, Utilde, Stilde,	// QR fact. for Bi22, non-permuted
					r_Bi, 0, printWhat, printWhere, pof );	// Stilde is an extra out-parameter

				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\n\nIntermediate rank value (rank of Bbar_i,22) = %d\n", r_Bi );
					else if( printWhere == 0 )
					{
						printf( magenta );
						printf( "\n\nIntermediate rank value (rank of Bbar_i,22) = %d\n", r_Bi );
						printf( normal );
					}
				}
				dim += r_Bi;								// update global counter
/*done?*/		if( dim == m_dae )							// we are done!
				{
					if( printWhat > 1 )
						if( printWhere == 1 )
						{
							fprintf( pof, "\ni = %d, dim = %d, m_dae = %d", niter, dim, m_dae );
							fprintf( pof, "\nThe rank is = %d", dim );
							fprintf( pof, "\nWe are done! Terminating while-cycle...\n\n" );
						}
						else if( printWhere == 0 )
						{
							printf( "\ni = %d, dim = %d, m_dae = %d", niter, dim, m_dae );
							printf( "\nThe rank is = %d", dim );
							printf( "\nWe are done! Terminating while-cycle...\n\n" );
						}
					break;									// terminate while-cycle
				}
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\ndim = %d < m_dae = %d", dim, m_dae );
						fprintf( pof, "\nContinue with the while cycle...\n", dim, m_dae );
					}
					else if( printWhere == 0 )
					{
						printf( "\ndim = %d < m_dae = %d", dim, m_dae );
						printf( "\nContinue with the while cycle...\n", dim, m_dae );
					}
				}
				//
				// D E C O M P O S E   Stilde[m_dae-r][m_dae-r]
				//
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nDecomposing Stilde and calculating Ri1m, Ri1m*Ri2..." );
					else if( printWhere == 0 )
						printf( "\nDecomposing Stilde and calculating Ri1m, Ri1m*Ri2..." );
				}
				if( r_Bi > 0 )
				{
					Ri1.redim( r_Bi, r_Bi, TPgrade );
					Ri2.redim( r_Bi, Stilde.ncols() - r_Bi, TPgrade );
				}
				decompS( r_Bi, Stilde, PItilde, Ri1, Ri2 );
				Vtilde.redim( Stilde.nrows(), Stilde.ncols(), TPgrade );
				if( r_Bi > 0 )								// Ri1 exists
				{
					Ri1m.redim( r_Bi, r_Bi, TPgrade );
					Ri1.trinvm( Ri1m );						// inverse of Ri1, upper triang.
					aux4.redim( r_Bi, r_Bi, TPgrade );
					aux4.mmCaABbC( 1.0, 1.0, Ri1, Ri1m );	//to test: Ri1*Ri1m = I ?? 
					if( 0 < Stilde.ncols() - r_Bi )			// there exists also Ri2
					{
						aux5.redim( r_Bi, Stilde.ncols() - r_Bi, TPgrade );
						aux5.mmCaABbC( 1.0, 1.0, Ri1m, Ri2 );// aux5 = Ri1m*Ri2
						aux5 = -aux5;
						//
						// C O M P U T E   Vtilde
						//
						computeV( PItilde, 0, aux5, Vtilde );// rows permutation
					}
					else
					{
						Vtilde.set2Id();					// Vtilde is set to the identity
						Vtilde.rpermutem( PItilde );		// rows permutation
					}
				}
				else										// Ri1 does not exist but Ri2
				{
					Vtilde.set2Id();						// Vtilde is set to the identity
					Vtilde.rpermutem( PItilde );			// rows permutation
				}
/*(S.6b)*/ 
				aux4.redim( m_dae, m_dae, TPgrade );
				aux5.redim( m_dae, m_dae, TPgrade );
				//
				// C O M P U T E   B_i+1
				//
				if( niter > 0 )
				{
					if( printWhat > 1 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\nApplying shift operator to DPiDm to calculate (DPiDm)'..." );
						else if( printWhere == 0 )
							printf( "\nApplying shift operator to DPiDm to calculate (DPiDm)'..." );
					}
					//
					// SHIFT OPERATOR on DPiDm  (DPiDm[m_dyn][m_dyn])
					//
					if( niter == 1 )						// initializing or saving old variable
						R.shift( dDPiDm_old );
					else
						dDPiDm_old = dDPiDm;
					DPiDm.shift( dDPiDm );
					if( printWhat > 1 )
					{
						if( printWhere == 1 )
						{
							DPiDm.fprinttpm( pof, "\nDPiDm = \n", ioeps );
							dDPiDm.fprinttpm( pof, "\ndDPiDm = DPiDm shifted = (DPiDm)' = \n", ioeps );
							PPi_old.fprinttpm( pof, "\nPPi_{i-1} = \n", ioeps );
							dDPiDm_old.fprinttpm( pof, "\n(DP_{i-1}Dm)' = \n", ioeps );
						}
						else if( printWhere == 0 )
						{
							DPiDm.printtpm( "\nDPiDm = \n", black, ioeps );
							dDPiDm.printtpm( "\ndDPiDm = DPiDm shifted = (DPiDm)' = \n", black, ioeps );
							PPi_old.printtpm( "\nPPi_{i-1} = \n", black, ioeps );
							dDPiDm_old.printtpm( "\n(DP_{i-1}Dm)' = \n", black, ioeps );
						}
					}
					//
					// Rule to derivate the product:
					// 		DerivPPi = dDPiDm_old * DPiDm  +  PPi_old * dDPiDm
					// <=>
					//		(PP_{i-1}DPiDm)' = (DP_{i-1}Dm)' * DPiDm  +  PPi_{i-1} * (DPiDm)'
					//
					DerivPPi.mmCaABbC( 1.0, 1.0, dDPiDm_old, DPiDm );
					aux1.mmCaABbC( 1.0, 1.0, PPi_old, dDPiDm );//aux1[m_dyn][m_dyn]
					DerivPPi += aux1;
					/*ONLY PROVISIONALLY, for testing purposes!!!!!!!!!*/ //DerivPPi.set2zero();
					if( printWhat > 1 )
					{
						if( printWhere == 1 )
							DerivPPi.fprinttpm( pof, "\nDerivPPi = (PP_{i-1}DPiDm)' = \n", ioeps );
						else if( printWhere == 0 )
							DerivPPi.printtpm( "\nDerivPPi = (PP_{i-1}DPiDm)' = \n", black, ioeps );
					}
					aux4.mmCaABbC( 1.0, 1.0, Bi_old, Pi_old );// aux4[m_dae][m_dae]
					aux3.mmCaABbC( 1.0, 1.0, Gi, Dm );		// aux3[m_dae][m_dyn]
					aux6.mmCaABbC( 1.0, 1.0, aux3, DerivPPi );// aux6[m_dae][m_dyn]
					aux3.mmCaABbC( 1.0, 1.0, aux6, PPi_old );// aux3[m_dae][m_dyn]
					aux5.mmCaABbC( 1.0, 1.0, aux3, DD );	// aux5[m_dae][m_dae]
					Bi_old = Bi;
					Bi = aux4 - aux5;						// update Bi
					if( printWhat > 21 )
					{
						if( printWhere == 1 )
							Bi.fprinttpm( pof, "\nBi = B_i-1*P_i-1 - Gi*Dm*DerivPPi*PP_i-1*D\n", ioeps );
						else if( printWhere == 0 )
							Bi.printtpm( "\nBi = B_i-1*P_i-1 - Gi*Dm*DerivPPi*PP_i-1*D\n", black, ioeps );
					}
				}
/*(S.5)*/
				aux4.mmCaABbC( 1.0, 1.0, Bi, Qi[ niter ] );	// Bi * Qi
				//Bi.printtpm( "\nBi = \n", green, ioeps );
				//Qi[ niter ].printtpm( "\nQi = \n", blue, ioeps );
				//aux4.printtpm( "\nBi * Qi = \n", red, ioeps );
				aux5.mmCaABbC( 1.0, 1.0, Gim, aux4 );		// Gim * Bi*Qi
				I.set2Id();
				aux5 = I - aux5;							// Fim = I - Gim*Bi*Qi
				//
				// C O M P U T E   n e w   Gi
				//
				Gi += aux4;									// Gi = Gi + Bi*Qi
				mi2.redim( dim - r_Bi, m_dae - dim + r_Bi, TPgrade );
				mi2.set2zero();
				computemTui( mi2, Si, Utilde, mTui );		// matrix * Tui
				Uim_new.mmCaABbC( 1.0, 1.0, mTui, Uim );	// U_i+1m
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
					{
						//aux5.fprinttpm( pof, "\nFim = I - Gim*Bi*Qi = \n", ioeps );
						//Gi.fprinttpm( pof, "\nG_i+1 = Gi + Bi*Qi\n", ioeps );
						//Uim_new.fprinttpm( pof, "\nU_i+1m = matrix * Tui * Uim\n", ioeps );
						fprintf( pof, "\nComputing Zi..." );
					}
					else if( printWhere == 0 )
					{
						//aux5.printtpm( "\nFim = I - Gim*Bi*Qi = \n", black, ioeps );
						//Gi.printtpm( "\nG_i+1 = Gi + Bi*Qi\n", black, ioeps );
						//Uim_new.printtpm( "\nU_i+1m = matrix * Tui * Uim\n", black, ioeps );
						printf( "\nComputing Zi..." );
					}
				}
				//
				// C O M P U T E   Zi
				//
				aux4.redim( m_dae, m_dae, TPgrade );
				aux4.mmCaABbC( 1.0, 1.0, Uim_new, Bi );		// U_i+1m * Bi
				Zi.setncols( Vi2.ncols() );
				Zi.mmCaABbC( 1.0, 1.0, aux4, Vi2 );			// U_i+1m*Bi * Vi2
				ncolZ += Vi2.ncols();						// update nr. of columns from Z
				Z.setncols( ncolZ );
				filloutM( ncolZ - Vi2.ncols(), Vi2.ncols(), Zi, Z );// Z = [Z|Zi]
				//
				// C O M P U T E   V_i+1
				//
				computeTlim( mi1, Si, Vtilde, Tlim );		// Tli * matrix
				aux4.mmCaABbC( 1.0, 1.0, aux5, Vi );		// Fim * Vi
				Vi.mmCaABbC( 1.0, 1.0, aux4, Tlim );		// Fim*Vi * Tlim
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						Zi.fprinttpm( pof, "\nZi = U_i+1m*Bi*Vi2 = \n", ioeps );
						Z.fprinttpm( pof, "\nZ = [Z|Zi] =\n", ioeps );
					}
					else if( printWhere == 0 )
					{
						Zi.printtpm( "\nZi = U_i+1m*Bi*Vi2 = \n", black, ioeps );
						Z.printtpm( "\nZ = [Z|Zi] =\n", black, ioeps );
					}
				}
				//
				// C O M P U T E  Si_new
				//
				if( niter == 0 )
				{
					Sim.setdim( dim - r_Bi, dim - r_Bi, TPgrade );
					Si_new.setdim( dim - r_Bi, dim - r_Bi, TPgrade );
					Sim_new.setdim( dim - r_Bi, dim - r_Bi, TPgrade );
					filloutS( 0, Si, Si_new );				// initialize Si_new
					Si.trinvm( Sim );						// inverse of Si, upper triang.
					filloutS( 0, Sim, Sim_new );			// initialize Sim_new
				}
				else
				{
					if( r_Bi > 0 )							// Ri1, Ri1m exists
					{
						filloutS( r_Bi, Ri1, Si_new );
						filloutS( r_Bi, Ri1m, Sim_new );
					}
				}
				//
				// C O M P U T E   Wihat
				// (using Gaussian elimination with scalling and column pivoting)
				//
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nComputing Wihat using Gaussian elimin. with scalling and col. pivoting..." );
					else if( printWhere == 0 )
						printf( "\nComputing Wihat using Gaussian elimin. with scalling and col. pivoting..." );
				}
				Wihat.setncols( ncolb );
				Wihat = b;									// b would be overwritten in 'gsolve'
				aux5 = Vi;									// V_i+1 would be overwritten in 'gsolve'
				aux5.gsolve( Wihat );						// V_i+1*Wihat = b
				//
				// T E S T:   V_i+1*Wihat - b == 0 ???
				//
				aux4.redim( Vi.nrows(), Wihat.ncols(), TPgrade );
				aux4.mmCaABbC( 1.0, 1.0, Vi, Wihat );
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\nT E S T:   V_i+1*Wihat - b == 0 ???  (eps=%.16lg)\n", eps );
						//aux4.fprinttpm( pof, "\nV_i+1*Wihat = \n", ioeps );
						//b.fprinttpm( pof, "\nb = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						printf( "\nT E S T:   V_i+1*Wihat - b == 0 ???  (eps=%.16lg)\n", eps );
						//aux4.printtpm( "\nV_i+1*Wihat = \n", black, ioeps );
						//b.printtpm( "\nb = \n", black, ioeps );
					}
				}
				aux4 -= b;
				if( !aux4.isZero( eps ) )
					if( printWhat > 0 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\n\nWARNING  -->  V_i+1*Wihat <> b (for eps=%.16lg)\n", eps );
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\n\nWARNING  -->  V_i+1*Wihat <> b (for eps=%.16lg)\n", eps );
							printf( normal );
						}
					}
				//
				// D E C O M P O S I N G   V   a n d   U P D A T  I N G   b
				//
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nDecomposing V_i+1=[V_i+1,1|V_i+1,2] and updating b=[b|V_i+1,2]..." );
					else if( printWhere == 0 )
						printf( "\nDecomposing V_i+1=[V_i+1,1|V_i+1,2] and updating b=[b|V_i+1,2]..." );
				}
				Vi1.setncols( dim );
				Vi2.setncols( m_dae - dim );
				decompV( dim, Vi, Vi1, Vi2 );
				ncolb += m_dae - dim;						// update nr. of columns from b
				b.setncols( ncolb );
				filloutM( ncolb - m_dae + dim, m_dae - dim, Vi2, b );// update b: b=[b|V_i+1,2]
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						Vi1.fprinttpm( pof, "\nV_i+1,1 = \n", ioeps );
						Vi2.fprinttpm( pof, "\nV_i+1,2 = \n", ioeps );
						b.fprinttpm( pof, "\nb updated = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						Vi1.printtpm( "\nV_i+1,1 = \n", black, ioeps );
						Vi2.printtpm( "\nV_i+1,2 = \n", black, ioeps );
						b.printtpm( "\nb updated = \n", black, ioeps );
					}
				}
/*(S.6)*/
/* First part (calculation of Zi and Z): see above!!*/
/*(S.6a)*/
				Ztilde = Z;
				Ztilde.setnrows( dim );						// Z = ( Ztilde )
															//     (    0   )  <- not interesting
				//
				// T E S T:   What_1 == Sim_new*Ztilde ???
				//
				aux5.redim( Sim_new.nrows(), Ztilde.ncols(), TPgrade );
				aux5.mmCaABbC( 1.0, 1.0, Sim_new, Ztilde );	// Sim_new*Ztilde
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\nT E S T:   What_1 == Sim_new*Ztilde ???  (eps=%.16lg)\n", eps );
						//aux5.fprinttpm( pof, "\nSim_new*Ztilde = \n", ioeps );
						//Wihat.fprinttpm( pof, "\nWihat = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						printf( "\nT E S T:   What_1 == Sim_new*Ztilde ???  (eps=%.16lg)\n", eps );
						//aux5.printtpm( "\nSim_new*Ztilde = \n", black, ioeps );
						//Wihat.printtpm( "\nWihat = \n", black, ioeps );
					}
				}
				if( !aux5.mcompare( Wihat, aux5.nrows(), eps ) )// Sim_new*Ztilde == What_1 ?
					if( printWhat > 0 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\n\nWARNING  -->  What_1 <> Sim_new*Ztilde (for eps=%.16lg)\n", eps );
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\n\nWARNING  -->  What_1 <> Sim_new*Ztilde (for eps=%.16lg)\n", eps );
							printf( normal );
						}
					}
				if( printWhat > 1 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nApplying Householder to Ztilde...\n" );
					else if( printWhere == 0 )
						printf( "\nApplying Householder to Ztilde...\n" );
				}
				//
				// 5th   H O U S E H O L D E R   Q R   F A C T O R I Z A T I O N
				//
				SZtilde.redim( Ztilde.nrows(), Ztilde.ncols(), TPgrade );
				UZtilde.redim( Ztilde.nrows(), Ztilde.nrows(), TPgrade );
				//Ztilde.printtpm( "\nZtilde before QR\n", black, ioeps );
/*(HH5)*/		hqrfcp( eps, Ztilde, PIZtilde, UZtilde, SZtilde,// QR fact. for Ztilde, non-permuted
					r_Zt, 0, printWhat, printWhere, pof );	// SZtilde is an extra out-parameter
				//Ztilde.printtpm( "\nZtilde after QR\n", black, ioeps );
				//UZtilde.printtpm( "\nUZtilde after QR\n", black, ioeps );
				//SZtilde.printtpm( "\nSZtilde after QR\n", black, ioeps );
				//
				// T E S T:   r_Zt == Ztilde.ncols() ???
				//
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\n\nT E S T:   r_Zt == Ztilde.ncols() ???  (eps=%.16lg)\n", eps );
						//fprintf( pof, "\nr_Zt = %d\n", r_Zt );
					}
					else if( printWhere == 0 )
					{
						printf( "\n\nT E S T:   r_Zt == Ztilde.ncols() ???  (eps=%.16lg)\n", eps );
						//printf( "\nr_Zt = %d\n", r_Zt );
					}
				}
				if( r_Zt != Ztilde.ncols() )
				{
					if( printWhat > 0 )
					{
						if( printWhere == 1 )
						{
							fprintf( pof, "\nERROR  -->  Rank of Ztilde = %d <> %d = Nr. columns of Ztilde", r_Zt, Ztilde.ncols() );
							fprintf( pof, "\nSTOP.\n" );
						}
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\nERROR  -->  Rank of Ztilde = %d <> %d = Nr. columns of Ztilde", r_Zt, Ztilde.ncols() );
							printf( "\nSTOP.\n" );
							printf( normal );
						}
					}
					success = 1;
					break;
				}
				Rtilde.redim( SZtilde.nrows(), SZtilde.ncols(), TPgrade );
				Rtilde = SZtilde;
				Rtilde.setnrows( r_Zt );					// only its upper part is interesting
				Rtilde.cpermutem( PIZtilde );
				//
				// D E C O M P O S E   Wihat
				//
				Wihat1.redim( SZtilde.nrows(), Wihat.ncols(), TPgrade );
				Wihat2.redim( Wihat.nrows() - SZtilde.nrows(), Wihat.ncols(), TPgrade );
				decompW( Ztilde.nrows(), Wihat, Wihat1, Wihat2 );
				//
				// C O M P U T E   y
				//
				Wihat2.cpermutem( PIZtilde );
				aux4.redim( Wihat2.nrows(), Wihat2.ncols(), TPgrade );
				aux5.redim( Rtilde.nrows(), Rtilde.ncols(), TPgrade );
				aux4 = Wihat2;								// save for later use
				aux5 = Rtilde;
				Rtilde.utxsolve( Wihat2 );					// y * Rtilde = Wihat2 (X*A = B)
				//aux4.printtpm( "\noriginal Wihat2 = \n", black, ioeps );
				//aux5.printtpm( "\noriginal Rtilde = \n", black, ioeps );
				//
				// T E S T:   Wihat2 - y*Rtilde == 0 ???
				//
				aux8.redim( Wihat2.nrows(), aux5.ncols(), TPgrade );
				aux8.mmCaABbC( 1.0, 1.0, Wihat2, aux5 );
				//Wihat2.printtpm( "\ny = \n", black, ioeps );
				//aux8.printtpm( "\ny*Rtilde = \n", black, ioeps );
				aux4 -= aux8;
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\nT E S T:   Wihat2 - y*Rtilde == 0 ???  (eps=%.16lg)\n", eps );
						//aux4.fprinttpm( pof, "\nWihat2 - y*Rtilde = \n", ioeps );
					}
					else if( printWhere == 0 )
					{
						printf( "\nT E S T:   Wihat2 - y*Rtilde == 0 ???  (eps=%.16lg)\n", eps );
						//aux4.printtpm( "\nWihat2 - y*Rtilde = \n", black, ioeps );
					}
				}
				if( !aux4.isZero( eps ) )
					if( printWhat > 0 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\n\nWARNING  -->  Wihat2 <> y*Rtilde (for eps=%.16lg)", eps );
						else if( printWhere == 0 )
						{
							printf( red );
							printf( "\n\nWARNING  -->  Wihat2 <> y*Rtilde (for eps=%.16lg)", eps );
							printf( normal );
						}
					}
				//
				// C O M P U T E   m_i+1,1
				//
				mi1.redim( Wihat2.nrows(), UZtilde.nrows(), TPgrade );
				if( UZtilde.ncols() - Wihat2.ncols() > 0 )	// mtilde exists
				{
					mtilde.redim( Wihat2.nrows(), UZtilde.ncols() - Wihat2.ncols(), TPgrade );
					mtilde.set2zero();						// mtilde = 0, free parameter
					computemi1( Wihat2, UZtilde, mtilde, mi1 );
				}
				else										// use previous dimensions
				{
					computemi1( Wihat2, UZtilde, mtilde, mi1 );
				}
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
						fprintf( pof, "\nComputing Gim_new...\n" );
					else if( printWhere == 0 )
						printf( "\nComputing Gim_new...\n" );
				}
/*(S.6 cont.)*/		
				//
				// C O M P U T E  the new Gm
				//
				mi2.redim( Sim_new.nrows(), m_dae - Sim_new.ncols(), TPgrade );
				mi2.set2zero();								// mi2 = 0, free parameter
				computeGm( mi1, mi2, Vi, Si_new, Sim_new, Uim, Gim );
				if( printWhat > 21 )
				{
					if( printWhere == 1 )
						Gim.fprinttpm( pof, "\nGim = \n", ioeps );
					else if( printWhere == 0 )
						Gim.printtpm( "\nGim = \n", black, ioeps );
				}
				Si = Si_new;								// update for next time
				Sim = Sim_new;
				if( (dim == m_dae)  ||  (niter > m_dae) )
				{
					//
					// T E S T:  Gi*Qi == 0 ??
					//
					aux4.mmCaABbC( 1.0, 1.0, Gi, Qi[ niter ] );	// Gi * Qi
					if( printWhat > 21 )
					{
						if( printWhere == 1 )
							fprintf( pof, "\nT E S T:   Gi*Qi == 0 ???  (i=%d, eps=%.16lg)\n", niter, eps );
						else if( printWhere == 0 )
							printf( "\nT E S T:   Gi*Qi == 0 ???  (i=%d, eps=%.16lg)\n", niter, eps );
					}
					if( printWhat == 1 )
					{
						fnorm = aux4.fnorm();				// Frobenius norm, first Taylor coeff.
						if( fnorm > eps )
						{
							if( printWhere == 1 )
							{
								fprintf( pof, "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
								fprintf( pof, "\nFrobenius norm of Gi*Qi = %.16lg > eps", fnorm );
								fprintf( pof, "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
							}
							else if( printWhere == 0 )
							{
								printf( red );
								printf( "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
								printf( "\nFrobenius norm of Gi*Qi = %.16lg > eps", fnorm );
								printf( "\n(Frob. norm of the function evaluation, i.e., considering only the first Taylor coeff.)\n" );
								printf( normal );
							}
						}
					}
					else if( printWhat > 1 )
					{
						aux4.tcfnorm( TPgrade + 1, fn );	// Frobenius norms, all Taylor coeff.
						cond = vleqeps( TPgrade + 1, fn, eps, elem );
						if( !cond )							// at least one element is greater
						{									// than the epsilon
							if( printWhere == 1 )
							{
								fprintf( pof, "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
								fprintf( pof, "\nAt least for one Taylor coefficient holds: the Frobenius norm of Gi*Qi > eps" );
								fprintdv( pof, TPgrade + 1, fn, "\nFrobenius norms =\n" );
								if( printWhat > 21 )
									aux4.fprinttpm( pof, "\nGi*Qi = \n", ioeps );
							}
							else if( printWhere == 0 )
							{
								printf( red );
								printf( "\nWARNING  -->  |Gi*Qi| > eps (for eps=%.16lg)\n", eps );
								printf( "\nAt least for one Taylor coefficient holds: the Frobenius norm of Gi*Qi > eps" );
								printf( normal );
								printdv( TPgrade + 1, fn, "\nFrobenius norms =\n", red );
								if( printWhat > 21 )
									aux4.printtpm( "\nGi*Qi = \n", red, ioeps );
							}
						}
					}
				}
			} // END OF while-cycle
			if( success == 0 )
			{
				if( printWhat > 0 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "**************\n" );
						fprintf( pof, "The INDEX is = %d\n", niter + 1 );
						fprintf( pof, "**************\n" );
					}
					else if( printWhere == 0 )
					{
						printf( magenta );
						printf( "**************\n" );
						printf( "The INDEX is = %d\n", niter + 1 );
						printf( "**************\n" );
						printf( normal );
					}
				}
			}
			else if( success == 1 )							// Rank of Ztilde <> Nr. col. of Ztilde
			{
				if( printWhat > 0 )
				{
					if( printWhere == 1 )
					{
						fprintf( pof, "\nThe Matrix Sequence Loop finished." );
						fprintf( pof, "\n\n*************\n" );
						fprintf( pof, "The RANK is = %d\n", dim );
						fprintf( pof, "*************\n" );
						fprintf( pof, "****************************\n" );
						fprintf( pof, "A singular point was found!!\n" );
						fprintf( pof, "****************************\n" );
					}
					else if( printWhere == 0 )
					{
						printf( "\nThe Matrix Sequence Loop finished." );
						printf( magenta );
						printf( "\n\n*************\n" );
						printf( "The RANK is = %d\n", dim );
						printf( "*************\n" );
						printf( "****************************\n" );
						printf( "A singular point was found!!\n" );
						printf( "****************************\n" );
						printf( normal );
					}
				}
			}
		
		} // END OF big else, when r < m_dae

		delete [] Qi;										// Deallocating local variables...
		delete [] PIZtilde;
		delete [] PItilde;
		delete [] PI_G0;
		delete [] PI_D;
		delete [] PI_A;
	} // try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including the rethrowed ones
	{
		IDException e( "Error when initializing the matrix sequence.", 8 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Defines an example of active section for testing Taylor Arithmetic later.
 * 
 */
int asectpol( int m, int n, int d, double x, double **X_f, double **Y_f )
{
	int i;
	adouble ax;
	double *px;
	adouble *pax;
	
	try
	{
		px = new double[ m ];								// dep. vector with doubles 
		pax = new adouble[ m ];								// dep. vector with adoubles
		if( (px == 0)  ||  (px == NULL)  ||  (pax == 0)  ||  (pax == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		//
		// ACTIVE SECTION 4
		//
		trace_on( 4 );										// mark active section
		ax <<= x;											// 1 indep. var.: 
															// ax is active, x is passive 
		pax[ 0 ] = sin( ax );
		pax[ 1 ] = cos( ax );
		pax[ 2 ] = sin( ax ) + cos( ax );
		pax[ 3 ] = sin( ax ) - cos( ax );
		pax[ 4 ] = sin( ax ) * sin( ax );
		pax[ 5 ] = cos( ax ) * cos( ax );
		pax[ 6 ] = sin( ax ) * sin( ax ) + cos( ax ) * cos( ax );
		pax[ 7 ] = sin( ax ) / cos( ax );
		pax[ 8 ] = sin( ax ) * 2;
		pax[ 9 ] = sqrt( sin( ax ) );

		for( i = 0; i < m; i++ )
			pax[ i ] >>= px[ i ];							// m_tra dep. var.: pax is active, 
															// px is passive
		trace_off();										// end active section
		printdv( m, px, "\nf(x): R in Rn\n", black );
		//
		// INITIALIZE TAYLOR COEFFICIENTS
		//
		X_f[ 0 ][ 0 ] = ax.value();							// the function depends on...
		X_f[ 0 ][ 1 ] = 1.0;								// ...one parameter only, i.e., x
		for( i = 2; i < d + 1; i++ )
			X_f[ 0 ][ i ] = 0.0;
		//
		// FORWARD MODE
		//
		forward( 4, m, n, d, d+1, X_f, Y_f );				// tape = 4
		printdm( n, d + 1, X_f, "\nX_f[n][d+1] (Taylor coeff. X of indep. var. x) = \n", yellow );
		printdm( m, d + 1, Y_f, "\nY_f[m][d+1] (Taylor coeff. Y of dep. var. f(x)) = \n", yellow );
		
		delete [] pax;
		delete [] px;
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including rethrowed ones 
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Defines an example of active section for testing Taylor Arithmetic later.
 * 
 */
int asectpol2( int m, int n, int d, double x, double **X_f, double **Y_f )
{
	int i;
	adouble ax;
	double *px;
	adouble *pax;
	
	try
	{
		px = new double[ m ];								// dep. vector with doubles 
		pax = new adouble[ m ];								// dep. vector with adoubles
		if( (px == 0)  ||  (px == NULL)  ||  (pax == 0)  ||  (pax == NULL) ) 
			throw IDException( "Memory allocation failure.", 3 );
		//
		// ACTIVE SECTION 4
		//
		trace_on( 5 );										// mark active section
		ax <<= x;											// 1 indep. var.: 
															// ax is active, x is passive 
		pax[ 0 ] = 1;
		pax[ 1 ] = cos( ax );
		pax[ 2 ] = sin( ax );
		pax[ 3 ] = 0;
		pax[ 4 ] = 0;
		pax[ 5 ] = sqrt( cos( ax ) );

		for( i = 0; i < m; i++ )
			pax[ i ] >>= px[ i ];							// m_tra dep. var.: pax is active, 
															// px is passive
		trace_off();										// end active section
		printdv( m, px, "\nf(x): R in Rn\n", black );
		//
		// INITIALIZE TAYLOR COEFFICIENTS
		//
		X_f[ 0 ][ 0 ] = ax.value();							// the function depends on...
		X_f[ 0 ][ 1 ] = 1.0;								// ...one parameter only, i.e., x
		for( i = 2; i < d + 1; i++ )
			X_f[ 0 ][ i ] = 0.0;
		//
		// FORWARD MODE
		//
		forward( 5, m, n, d, d+1, X_f, Y_f );				// tape = 5
		printdm( n, d + 1, X_f, "\nX_f[n][d+1] (Taylor coeff. X of indep. var. x) = \n", yellow );
		printdm( m, d + 1, Y_f, "\nY_f[m][d+1] (Taylor coeff. Y of dep. var. f(x)) = \n", yellow );
		
		delete [] pax;
		delete [] px;
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including rethrowed ones 
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Tests the Taylor arithmetic.
 * Tests operations over/with Taylor polynomials (TPolyn class).
 * 
 */
int testTPolyn()
{
	TPolyn u( 10 ), u2( 10 ), uu( 10 ), v( 10 ), v2( 10 ), vv( 10 ),
		w1( 10 ), w2( 10 ), usv( 10 ), umv( 10 );
	int i, j, m, n, d;
	double x;
	double **X_f;											// Taylor coeff. X[n][d+1], indep. var.
	double **Y_f;											// Taylor coeff. Y[m][d+1], dep. var.
	
	try
	{
		m = 10;												// 10 equations in 'asectpol(...)'
		n = 1;												// 1 indep. variable, i.e., x
		x = 1.0;											// initialization for x
		d = 10;												// derivative degree
		//
		// Variables declarations (Taylor coefficients)
		//
		X_f = new double*[ n ];
		for( i = 0; i < n; i++ )
			X_f[ i ] = new double[ d + 1 ];
		Y_f = new double*[ m ];
		for( i = 0; i < m; i++ )
			Y_f[ i ] = new double[ d + 1 ];
		//
		// Active section: using ADOL-C
		//
		asectpol( m, n, d, x, X_f, Y_f );					// active section for testing
		//
		// Variables declarations (Taylor polynomials)
		//
		printf( "\n\nTaylor coefficients of..." );
		printf( "\n\nu = sin(x) =\n" );
		u.setCoeffs( Y_f[ 0 ] );							// u = sin(x)
		u.print();
		printf( "\n\nv = cos(x) =\n" );
		v.setCoeffs( Y_f[ 1 ] );							// v = cos(x)
		v.print();
		//
		// Testing Taylor arithmetic
		//
		printf( "\n\nuu = u;  uu += cos(x) =\n" );
		uu = u;
		uu += v;											// uu += cos(x)  <=>  sin(x)+cos(x)
		uu.print();
		printf( "\n\nvv = v;  vv += sin(x) =\n" );
		vv = v;
		vv += u;											// vv += sin(x)  <=>  cos(x)+sin(x)
		vv.print();
		printf( "\n\nusv = sin(x)+cos(x) =\n" );
		usv = u + v;										// usv = sin(x)+cos(x)
		usv.print();
		printf( "\n\nuu = u;  uu -= cos(x) =\n" );
		uu = u;
		uu -= v;											// uu -= cos(x)  <=>  sin(x)-cos(x)
		uu.print();
		printf( "\n\nvv = v;  vv -= sin(x) =\n" );
		vv = v;
		vv -= u;											// vv -= sin(x)  <=>  cos(x)-sin(x)
		vv.print();
		printf( "\n\numv = sin(x)-cos(x) =\n" );
		umv = u - v;										// usv = sin(x)+cos(x)
		umv.print();
		printf( "\n\nu2 = sin(x)*sin(x) =\n" );
		u2 = u * u;											// u2 = sin(x)*sin(x)
		u2.print();
		printf( "\n\nu2 = u;  u2 *= sin(x) =\n" );
		u2 = u;
		u2 *= u;											// u2 *= sin(x)  <=>  sin(x)*sin(x)
		u2.print();
		printf( "\n\nv2 = cos(x)*cos(x) =\n" );
		v2 = v * v;											// v2 = cos(x)*cos(x)
		v2.print();
		printf( "\n\nw1 = sin(x)*sin(x) + cos(x)*cos(x) == 1 ??\n" );
		w1 = u2 + v2;										// w1 = sin(x)*sin(x) + cos(x)*cos(x)
		w1.print();
		printf( "\n\nu2 = sin^2 =\n" );
		u2 = u.sqr();										// u2 = sin^2
		u2.print();
		printf( "\n\nv2 = cos^2 =\n" );
		v2 = v.sqr();										// v2 = cos^2
		v2.print();
		printf( "\n\nw1 = sin^2 + cos^2 == 1 ??\n" );
		w1 = u2 + v2;										// w1 = sin^2 + cos^2
		w1.print();
		printf( "\n\nw2 = sin/cos =\n" );
		w2 = u / v;											// w2 = sin(x)/cos(x)
		w2.print();
		printf( "\n\nw2 = u;  w2 /= cos(x) =\n" );
		w2 = u;
		w2 /= v;											// w2 /= cos(x)  <=>  sin(x)/cos(x)
		w2.print();
		printf( "\n\nuu = sin(x)*2 =\n" );
		uu = u * 2;											// uu = sin(x)*2
		uu.print();
		printf( "\n\nuu = sqrt( sin(x) ) =\n" );
		uu = u.sqrt();										// uu = sqrt( sin(x) )
		uu.print();
		printf( "\n\nuu = sqrt( sin(x) )^2 =\n" );
		uu = uu * uu;										// uu = sqrt^2( sin(x) )
		uu.print();
		printf( "\n\nuu = -sin(x) =\n" );
		uu = -u;											// uu = -sin(x)
		uu.print();
		printf( "\n\nvv = -cos(x) =\n" );
		vv = -v;											// vv = -cos(x)
		vv.print();
		printf( "\n\nf[0]= %.16lg", vv.feval() );
		printf( "\n\n" );
		
	} // try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including rethrowed ones 
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Tests the QR factorization of specific TPolyn matrices, i.e., matrices of Taylor polynomials,
 * as well as functions for Taylor arithmetic.
 * 
 */
int testQRTPolynMatrix()
{
	int i, j, m, n, d, r_Zt;
	double coeffs00[ 10 ] = { 0.0 };						// all 10 values initialized to zero
	double coeffs01[ 10 ] = { 0.0 };
	double coeffs10[ 10 ] = { 0.0 };
	double coeffs11[ 10 ] = { 0.0 };
	double coeffs20[ 10 ] = { 0.0 };
	double coeffs21[ 10 ] = { 0.0 };
	double coeffs30[ 10 ] = { -2, 0, 2, -2, -0.83333333, 3.6666667, -2.6111111, -2.3333333, 
			5.8234127, -2.515873 };
	double coeffs31[ 10 ] = { 2, 0, -2, 2, 0.83333333, -3.6666667, 2.6111111, 2.3333333, 
			-5.8234127, 2.515873 };
	double coeffs40[ 10 ] = { 2, 0, -2, 2, 0.83333333, -3.6666667, 2.6111111, 2.3333333, 
			-5.8234127, 2.515873 };
	double coeffs41[ 10 ] = { -3, 0, 2, -2, -0.83333333, 3.6666667, -2.6111111, -2.3333333, 
			5.8234127, -2.515873 };
	double coeffs50[ 10 ] = { 3, 0, -3.5, 3.5, 1.3333333, -6.1666667, 4.3611111, 4.0833333, 
			-10.045139, 4.4027778 };
	double coeffs51[ 10 ] = { -3, 0, 3.5, -3.5, -1.3333333, 6.1666667, -4.3611111, -4.0833333, 
			10.045139, -4.4027778 };
	
	try
	{
		m = 6;												// matrix dimensions...
		n = 2;
		d = 10;												// nr. of Taylor coefficients

		TPolynMatrix Ztilde( m, n, d - 1 ),
			RZtilde( m, n, d - 1 ),							// = R after QR decomp. of Ztilde
			QZtilde( m, m, d - 1 );							// = Q after QR decomp. of Ztilde
		int *PZtilde = new int[ m ];						// permutations vector

		//
		// Initializations
		//
		Ztilde( 0, 0 ).setCoeffs( coeffs00 );				
		Ztilde( 0, 1 ).setCoeffs( coeffs01 );				
		Ztilde( 1, 0 ).setCoeffs( coeffs10 );				
		Ztilde( 1, 1 ).setCoeffs( coeffs11 );				
		Ztilde( 2, 0 ).setCoeffs( coeffs20 );				
		Ztilde( 2, 1 ).setCoeffs( coeffs21 );				
		Ztilde( 3, 0 ).setCoeffs( coeffs30 );				
		Ztilde( 3, 1 ).setCoeffs( coeffs31 );				
		Ztilde( 4, 0 ).setCoeffs( coeffs40 );				
		Ztilde( 4, 1 ).setCoeffs( coeffs41 );				
		Ztilde( 5, 0 ).setCoeffs( coeffs50 );				
		Ztilde( 5, 1 ).setCoeffs( coeffs51 );				
		Ztilde.printtpm( "\nZtilde = \n", black, ioeps );
		//
		// QR of Ztilde
		//
		hqrfcp( eps, Ztilde, PZtilde, QZtilde, RZtilde, r_Zt, 0, 0, printWhere, pof );

		//RZtilde.printtpm( "\nRZtilde = \n", black, ioeps );
		//QZtilde.printtpm( "\nQZtilde = \n", black, ioeps );
		//printiv( 2, PZtilde, "\nPZtilde = \n", black );
		
	} // try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including rethrowed ones 
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Tests the QR factorization of specific TPolyn matrices, i.e., matrices of Taylor polynomials,
 * as well as functions for Taylor arithmetic.
 * 
 */
int testQRTPolynMatrix2()
{
	int i, j, m, n, mX, nX, d, rA;
	double x;
	double **X_f;											// Taylor coeff. X[n][d+1], indep. var.
	double **Y_f;											// Taylor coeff. Y[m][d+1], dep. var.
	double coeffs[ 10 ] = { 1.0, 0.0 };
	TPolyn u( 9 ), v( 9 ), w( 9 ), p( 9 ), q( 9 );
	
	try
	{
		mX = 6;												// nr. of equations in 'asectpol2(...)'
		nX = 1;												// 1 indep. variable, i.e., x
		x = 0.0;											// initialization for x
		d = 10;												// derivative degree
		m = 3;												// matrix dimensions...
		n = 2;
		//
		// Variables declarations
		//
		X_f = new double*[ nX ];
		for( i = 0; i < nX; i++ )
			X_f[ i ] = new double[ d ];
		Y_f = new double*[ mX ];
		for( i = 0; i < mX; i++ )
			Y_f[ i ] = new double[ d ];
		TPolynMatrix A( m, n, d - 1 ),
			R( m, n, d - 1 ),								// = R after QR decomp. of A
			Q( m, m, d - 1 );								// = Q after QR decomp. of A
		int *P = new int[ m ];								// permutations vector

		//
		// Active section: using ADOL-C
		//
		asectpol2( mX, nX, d - 1, x, X_f, Y_f );			// active section for testing

		//
		// Initializations
		//
		A( 0, 0 ).setCoeffs( Y_f[ 0 ] );				
		A( 0, 1 ).setCoeffs( Y_f[ 1 ] );				
		A( 1, 0 ).setCoeffs( Y_f[ 2 ] );				
		A( 1, 1 ).setCoeffs( Y_f[ 3 ] );				
		A( 2, 0 ).setCoeffs( Y_f[ 4 ] );				
		A( 2, 1 ).setCoeffs( Y_f[ 5 ] );				
		A.printtpm( "\nA = \n", black, ioeps );
				
		//
		// QR of A
		//
		hqrfcp( eps, A, P, Q, R, rA, 0, 0, printWhere, pof );
		
		//
		// Other tests
		//
		/*printf( "\n\nu = \n" );
		u.setCoeffs( Y_f[ 1 ] ); //u.setCoeffs( coeffs );
		u.print();
		printf( "\n\nv = \n" );
		v.setCoeffs( Y_f[ 2 ] );
		v.print();
		printf( "\n\np = \n" );
		p.setCoeffs( Y_f[ 5 ] );
		p.print();
		printf( "\n\nw = u*v =\n" );
		w = u * v;
		w.print();
		printf( "\n\nw = - u*v =\n" );
		w = -w;
		w.print();
		printf( "\n\nw = - u*v / p =\n" );
		w /= p;
		w.print();
		printf( "\n\nu*v + p*w == 0???\n" );
		q = u * v + p * w;
		q.print();
		printf( "\n\n" );*/
				
	} // try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:" );
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// others, including rethrowed ones 
	{
		IDException e( "Error when initializing global variables.", 30 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Tests the Matrix class.
 * 
 */
void testMatrix()
{
/*	TPolynMatrix TP1( 2, 3 );
	printf( "\n\nTesting the class tpolynMatrix...\n");
	
	DoubleMatrix M1( 2, 3 ), M2( 2, 3 ), M3, M4( 3, 2 ), M5( 2, 2), M6;
	
	printf( "\n\nTesting the class Matrix...\n");
	
	M1( 0, 0 ) = 1.0;
	M1( 0, 1 ) = 2.0;
	M1( 0, 2 ) = 3.0;
	M1( 1, 0 ) = 4.0;
	M1( 1, 1 ) = 5.0;
	M1( 1, 2 ) = 6.0;
	M1.printm( "\nMatrix M1 = \n" );

	M2( 0, 0 ) = 1.0;
	M2( 0, 1 ) = 2.0;
	M2( 0, 2 ) = 3.0;
	M2( 1, 0 ) = 4.0;
	M2( 1, 1 ) = 5.0;
	M2( 1, 2 ) = 6.0;
	M2.printm( "\nMatrix M2 = \n" );
	
	M4( 0, 0 ) = 1.0;
	M4( 0, 1 ) = 2.0;
	M4( 1, 0 ) = 3.0;
	M4( 1, 1 ) = 4.0;
	M4( 2, 0 ) = 5.0;
	M4( 2, 1 ) = 6.0;
	M4.printm( "\nMatrix M4 = \n" );

	M3 = M1 + M2;
	M3.printm( "\nMatrix M3 = M1 + M2 = \n" );
	M6 = M3.transpm();
	M6.printm( "\nMatrix M6 = M3 transposed = \n" );

	M3 = M1 - M2;
	M3.printm( "\nMatrix M3 = M1 - M2 = \n" );
	
	M3 += M1;
	M3.printm( "\nMatrix M3 += M1 = \n" );

	M3 -= M1;
	M3.printm( "\nMatrix M3 -= M1 = \n" );
	
	M5.mmCaABbC( 1.0, 1.0, M1, M4 );
	M5.printm( "\nMatrix M5 = M1*M4 = \n" );
	
	int *piv = new int[2];
	double *v = new double[2];
	piv[ 0 ] = 1;
	piv[ 1 ] = 0;
	v[ 0 ] = 1.5;
	v[ 1 ] = 0.3;
	M5.cpermutem( piv );
	M5.printm( "\nMatrix M5 columns permuted = \n" );
	M5.rpermutem( piv );
	M5.printm( "\nMatrix M5 rows permuted = \n" );
	M5.set2Id();
	M5.printm( "\nMatrix M5 set to I = \n" );
	M5.set2zero();
	M5.printm( "\nMatrix M5 set to 0 = \n" );
	M5.set2val( 8.5 );
	M5.printm( "\nMatrix M5 set to a value = \n" );
	printdv( 2, v, "v = \n", yellow );
	M5.utsolve( v );
	printdv( 2, v, "M5*x=v  =>  x = \n", yellow );
*/
}

/**
 * Determines the index.
 * 
 */
int daeindex( int argc, char *argv[], IDExample * pobj )
{
	char ch;
	int errcode = 0;
	
	try
	{
		if( argc != 1 )										// Process the command line
			throw IDException( "Wrong command line.", 1 );

//		setInitialTime();									// Start time counter (see time.cpp/.h)
		initialize( pobj );									// Process global data and initialize them

		if( printWhat > 0 )
			if( printWhere == 1 )
				fprintf( pof, "\n*** COMPUTATIONS ***\n\n" );
			else  if( printWhere == 0 )
				printf( "\n*** COMPUTATIONS ***\n\n" );
		asectra( pobj );									// Active section and calc. of x(t)
		asecdyn( pobj ); 									// Active section and calc. of z(t)=d'(x(t),t)
		asecdae( pobj );									// Active section and calc. of f(z(t),x(t),t)
		ABD();												// Construct A, B, and D matrices
		matrixSqc( r );										// Matrix sequence processing

		if( printWhat > 0 )
			if( printWhere == 1 )
				fprintf( pof, "\n*** END OF COMPUTATIONS ***\n" );
			else  if( printWhere == 0 )
				printf( "\n*** END OF COMPUTATIONS ***\n" );
	}
	catch( int e )											// for rethrowed exceptions
	{
		printf( red );
		printf( "\n---Program aborted! Returned code = %d\n", e );
		printf( normal );
/**/		//cleanup();									// Cleaning operations
		return e;
	}
	catch( IDException e )
	{
        e.report();
/**/		//cleanup();									// Cleaning operations
		return e.getErrCode();
	}
	catch(...)
	{
		IDException e( "Unexpected exception type. Program aborted!", 80 );
        e.report();
/**/		//cleanup();									// Cleaning operations
		return e.getErrCode();
	}
/**/		//cleanup();									// Cleaning operations for global var.
//	setFinalTime();											// Stop time counter
//	double elapsedT = getElapsedTime();						// Get elapsed time
//	printf( "\n\nElapsed time: %.3f seconds.\n", elapsedT );
	//printf( "\nType return to continue...\n" );			// Wait some time
	//ch = getchar();
	printf( "\nSuccessful termination! Returned code = %d\n", 0 );
	return 0;
}

/**
 * Determines the index.
 * 
 */
int daeindex( int argc, char *argv[], IDExample * pobj, IDOptions op )
{
	char ch;
	int errcode = 0;
	
	try
	{
		if( argc != 1 )										// Process the command line
			throw IDException( "Wrong command line.", 1 );

//		setInitialTime();									// Start time counter (see time.cpp/.h)
		options = op;										// update values for the global options
		initialize( pobj );									// Process global data and initialize them
		if( printWhat > 0 )
			if( printWhere == 1 )
				fprintf( pof, "\n*** COMPUTATIONS ***\n\n" );
			else  if( printWhere == 0 )
				printf( "\n*** COMPUTATIONS ***\n\n" );
		asectra( pobj );									// Active section and calc. of x(t)
		asecdyn( pobj ); 									// Active section and calc. of z(t)=d'(x(t),t)
		asecdae( pobj );									// Active section and calc. of f(z(t),x(t),t)
		ABD();												// Construct A, B, and D matrices
		matrixSqc( r );										// Matrix sequence processing
		//testMatrix();
		//testTPolyn();
		//testQRTPolynMatrix();
		//testQRTPolynMatrix2();

		if( printWhat > 0 )
			if( printWhere == 1 )
				fprintf( pof, "\n*** END OF COMPUTATIONS ***\n" );
			else  if( printWhere == 0 )
				printf( "\n*** END OF COMPUTATIONS ***\n" );
	}
	catch( int e )											// for rethrowed exceptions
	{
		printf( red );
		printf( "\n---Program aborted! Returned code = %d\n", e );
		printf( normal );
		//cleanup();										// Cleaning operations
		if( printWhere == 1 )
			fclose( pof );									// Close output file
		return e;
	}
	catch( IDException e )
	{
        e.report();
		//cleanup();										// Cleaning operations
		if( printWhere == 1 )
			fclose( pof );									// Close output file
		return e.getErrCode();
	}
	catch(...)
	{
		IDException e( "Unexpected exception type. Program aborted!", 80 );
        e.report();
		//cleanup();										// Cleaning operations
		if( printWhere == 1 )
			fclose( pof );									// Close output file
		return e.getErrCode();
	}
	//cleanup();											// Cleaning operations for global var.
	if( printWhere == 1 )
		fclose( pof );										// Close output file
//	setFinalTime();											// Stop time counter
//	double elapsedT = getElapsedTime();						// Get elapsed time
//	printf( "\n\nElapsed time: %.3f seconds.\n", elapsedT );
	//printf( "\nType return to continue...\n" );			// Wait some time
	//ch = getchar();
	printf( "\nSuccessful termination! Returned code = %d\n", 0 );
	return 0;
}
