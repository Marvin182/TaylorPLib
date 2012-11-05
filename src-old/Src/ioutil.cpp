/*!
 * \file ioutil.cpp
 * \brief Input/output functions for vectors and matrices.
 * \author D. Monett Díaz
 * \date June, 2006
 * 
 * This class implements some input/output functions for vectors and matrices.
 * 
 * Included files:
 * 	\file ioutil.h header file.
 * 
 */

#include "ioutil.h"

/**
 * Prints out a vector of doubles in a given format.
 * 
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] str The character string to be printed out.
 * 
 */
void printdv( int n, double *v, char *str )
{
	printf( str );
	for( int i = 0; i < n; i++ )
		//printf( "%c(%c%.8lg%c)\n", '\t','\t', v[ i ], '\t' );
		printf( "%s%.8lg\n", "  ", v[ i ] );
}

/**
 * Prints out a vector of doubles in a given format.
 * 
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void printdv( int n, double *v, char *str, const char *const color )
{
	printf( color );
	printf( str );
	for( int i = 0; i < n; i++ )
		//printf( "%c(%c%.8lg%c)\n", '\t','\t', v[ i ], '\t' );
		printf( "%s%.8lg\n", "  ", v[ i ] );
	printf( "\n" );
	printf( normal );
}

/**
 * Prints out to a file a vector of doubles in a given format.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprintdv( FILE * fn, int n, double *v, char *str )
{
	fprintf( fn, str );
	for( int i = 0; i < n; i++ )
		//fprintf( fn, "%c(%c%.8lg%c)\n", '\t','\t', v[ i ], '\t' );
		fprintf( fn, "s%.8lg\n", "  ", v[ i ] );
	fprintf( fn, "\n" );
}

/**
 * Prints out a vector of doubles in a given format, with column pivoting.
 * 
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a v.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void printdv( int n, double *v, int *piv, char *str, const char *const color )
{
	printf( color );
	printf( "------\n" );
	printf( str );
	for( int i = 0; i < n; i++ )
		printf( "  %.8lg%c", v[ piv[ i ] ], '\t' );
	printf( "\n" );
	printf( normal );
}

/**
 * Prints out to a file a vector of doubles in a given format, with column pivoting.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a v.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprintdv( FILE * fn, int n, double *v, int *piv, char *str )
{
	fprintf( fn, "------\n" );
	fprintf( fn, str );
	for( int i = 0; i < n; i++ )
		fprintf( fn, "  %.8lg%c", v[ piv[ i ] ], '\t' );
	fprintf( fn, "\n" );
}

/**
 * Prints out a vector of doubles from component a to component b, in a given format.
 * 
 * \param[in] a The starting component in vector \a v.
 * \param[in] b The last component in vector \a v that should be printed out.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 *
 */
void printdv( int a, int b, double *v, char *str, const char *const color )
{
	printf( color );
	printf( "------\n" );
	printf( str );
	for( int i = a; i < b; i++ )
		printf( "  %.8lg%c", v[ i ], '\t' );
	printf( "\n" );
	printf( normal );
}

/**
 * Prints out to a file a vector of doubles from component a to component b, in a given format.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] a The starting component in vector \a v.
 * \param[in] b The last component in vector \a v that should be printed out.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] str The character string to be printed out.
 *
 */
void fprintdv( FILE * fn, int a, int b, double *v, char *str )
{
	fprintf( fn, "------\n" );
	fprintf( fn, str );
	for( int i = a; i < b; i++ )
		fprintf( fn, "  %.8lg%c", v[ i ], '\t' );
	fprintf( fn, "\n" );
}

/**
 * Prints out a vector of doubles from component a to component b, in a given format,
 * with column pivoting.
 * 
 * \param[in] a The starting component in vector \a v.
 * \param[in] b The last component in vector \a v that should be printed out.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a v.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 *
 */
void printdv( int a, int b, double *v, int *piv, char *str, const char *const color )
{
	printf( color );
	printf( "------\n" );
	printf( str );
	for( int i = a; i < b; i++ )
		printf( "  %.8lg%c", v[ piv[ i ] ], '\t' );
	printf( "\n" );
	printf( normal );
}

/**
 * Prints out to a file a vector of doubles from component a to component b, in a given format,
 * with column pivoting.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] a The starting component in vector \a v.
 * \param[in] b The last component in vector \a v that should be printed out.
 * \param[in] v The pointer to \a v, a vector of \type double.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a v.
 * \param[in] str The character string to be printed out.
 *
 */
void fprintdv( FILE * fn, int a, int b, double *v, int *piv, char *str )
{
	fprintf( fn, "------\n" );
	fprintf( fn, str );
	for( int i = a; i < b; i++ )
		fprintf( fn, "  %.8lg%c", v[ piv[ i ] ], '\t' );
	fprintf( fn, "\n" );
}

/**
 * Prints out a vector of integers in a given format.
 * 
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type integer.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void printiv( int n, int *v, char *str, const char *const color )
{
	printf( color );
	printf( "------\n" );
	printf( str );
	for( int i = 0; i < n; i++ )
		printf( "  %d%c", v[ i ], '\t' );
	printf( "\n" );
	printf( normal );
}

/**
 * Prints out to a file a vector of integers in a given format.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type integer.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprintiv( FILE * fn, int n, int *v, char *str )
{
	fprintf( fn, "------\n" );
	fprintf( fn, str );
	for( int i = 0; i < n; i++ )
		fprintf( fn, "  %d%c", v[ i ], '\t' );
	fprintf( fn, "\n" );
}

/**
 * Prints out a vector of integers in a given format, with column pivoting.
 * 
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type integer.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a v.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void printiv( int n, int *v, int *piv, char *str, const char *const color )
{
	printf( color );
	printf( "------\n" );
	printf( str );
	for( int i = 0; i < n; i++ )
		printf( "  %d%c", v[ piv[ i ] ], '\t' );
	printf( "\n" );
	printf( normal );
}

/**
 * Prints out to a file a vector of integers in a given format, with column pivoting.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] n The dimension of the vector \a v.
 * \param[in] v The pointer to \a v, a vector of \type integer.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a v.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprintiv( FILE * fn, int n, int *v, int *piv, char *str )
{
	fprintf( fn, "------\n" );
	fprintf( fn, str );
	for( int i = 0; i < n; i++ )
		fprintf( fn, "  %d%c", v[ piv[ i ] ], '\t' );
	fprintf( fn, "\n" );
}

/**
 * Prints out a matrix in a given format, with the previously set color.
 * 
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * 
 */
void printm( int m, int n, double **M, char *str )
{
	printf( yellow );
	printf( "------\n" );
	printf( str );
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
			printf( "%.8lg%c", M[ i ][ j ], '\t' );
		printf( "\n" );
	}
}

/**
 * Prints out to a file a matrix in a given format, with the previously set color.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprintm( FILE * fn, int m, int n, double **M, char *str )
{
	fprintf( fn, "------\n" );
	fprintf( fn, str );
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
			fprintf( fn, "%.8lg%c", M[ i ][ j ], '\t' );
		fprintf( fn, "\n" );
	}
}

/**
 * Prints out a matrix in a given format.
 * 
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void printm( int m, int n, double **M, char *str, const char *const color )
{
	printf( color );
	printf( "------\n" );
	printf( str );
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
			printf( "%.8lg%c", M[ i ][ j ], '\t' );
		printf( "\n" );
	}
	printf( normal );
}

/**
 * Prints out a matrix in a given format.
 * 
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * 
 */
void printdm( int m, int n, double **M, char *str )
{
	printf( str );
	for( int i = 0; i < m; i++ )
	{
		//printf( "%c(%c", '\t','\t' );
		printf( "%s", "  " );
		for( int j = 0; j < n; j++ )
			printf( "%.8lg%c", M[ i ][ j ], '\t' );
		//printf( ")\n" );
		printf( "\n" );
	}
}

/**
 * Prints out to a file a matrix in a given format.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprintdm( FILE * fn, int m, int n, double **M, char *str )
{
	fprintf( fn, str );
	for( int i = 0; i < m; i++ )
	{
		//fprintf( fn, "%c(%c", '\t','\t' );
		fprintf( fn, "%s", "  " );
		for( int j = 0; j < n; j++ )
			fprintf( fn, "%.8lg%c", M[ i ][ j ], '\t' );
		//fprintf( fn, ")\n" );
		fprintf( fn, "\n" );
	}
}

/**
 * Prints out a matrix in a given format.
 * 
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void printdm( int m, int n, double **M, char *str, const char *const color )
{
	printf( color );
	printf( str );
	for( int i = 0; i < m; i++ )
	{
		//printf( "%c(%c", '\t','\t' );
		printf( "%s", "  " );
		for( int j = 0; j < n; j++ )
			printf( "%.8lg%c", M[ i ][ j ], '\t' );
		//printf( ")\n" );
		printf( "\n" );
	}
	printf( normal );
}

/**
 * Prints out a matrix in a given format.
 * 
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void print3ddm( int m, int n, int p, double ***M, char *str, const char *const color )
{
	printf( color );
	printf( str );
	for( int k = 0; k < p; k++ )
	{
		for( int i = 0; i < m; i++  )
		{
			//printf( "%c(%c", '\t','\t' );
			printf( "%s", "  " );
			for( int j = 0; j < n; j++ )
				printf( "%.8lg%c", M[ i ][ j ][ k ], '\t' );
			//printf( ")\n" );
			printf( "\n" );
		}
		printf( "\n" );
	}
	printf( normal );
}

/**
 * Prints out to a file a matrix in a given format.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprint3ddm( FILE * fn, int m, int n, int p, double ***M, char *str )
{
	fprintf( fn, str );
	for( int k = 0; k < p; k++ )
	{
		for( int i = 0; i < m; i++  )
		{
			//fprintf( fn, "%c(%c", '\t','\t' );
			fprintf( fn, "%s", "  " );
			for( int j = 0; j < n; j++ )
				fprintf( fn, "%.8lg%c", M[ i ][ j ][ k ], '\t' );
			//fprintf( fn, ")\n" );
			fprintf( fn, "\n" );
		}
		fprintf( fn, "\n" );
	}
}

/**
 * Prints out a matrix in a given format, with column pivoting.
 * 
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
 * \param[in] str The character string to be printed out.
 * \param[in] color The desired color for the output.
 * 
 */
void printm( int m, int n, double **M, int *piv, char *str, const char *const color )
{
	printf( color );
	printf( "------\n" );
	printf( str );
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
			printf( "%.8lg%c", M[ i ][ piv[ j ] ], '\t' );
		printf( "\n" );
	}
	printf( normal );
}

/**
 * Prints out a matrix in a given format, with column pivoting.
 * 
 * \param[in] fn The file name where to write the data to.
 * \param[in] m The number of rows.
 * \param[in] n The number of columns.
 * \param[in] M The pointer to \a M, a matrix of \type double.
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
 * \param[in] str The character string to be printed out.
 * 
 */
void fprintm( FILE * fn, int m, int n, double **M, int *piv, char *str )
{
	fprintf( fn, "------\n" );
	fprintf( fn, str );
	for( int i = 0; i < m; i++ )
	{
		for( int j = 0; j < n; j++ )
			fprintf( fn, "%.8lg%c", M[ i ][ piv[ j ] ], '\t' );
		fprintf( fn, "\n" );
	}
}
