#include "Matrix.h"

using namespace std;
using namespace LibMatrix;

//
// C O N S T R U C T O R S
//

/**
 * Default constructor for the class. Creates the object.
 * 
 * It is a 1-by-1 matrix, i.e., it has only one element, which is set to zero:
 * 
 * 		m(0,0) == 0.0
 * 
 */
Matrix::Matrix():
	_rows(1),
	_cols(1),
	_dimT(0)
{
	allocateMemory(true);
}

/**
 * Constructor for the class with both the number of rows and columns as parameters. 
 * Creates the object.
 * 
 * \param[in] r The number of rows.
 * \param[in] c The number of columns.
 * 
 */
Matrix::Matrix(int rows, int cols, bool initialize): 
	_rows(rows), 
	_cols(cols),
	_dimT(0)
{
	allocateMemory(initialize);
}

// /**
//  * Constructor for the class with both the number of rows and columns as parameters,
//  * as well as the dimension of the elements' type, e.g., the Taylor polynomial's grade. 
//  * Creates the object.
//  * 
//  * \param[in] r The number of rows.
//  * \param[in] c The number of columns.
//  * \param[in] dimT The dimension of the type \type T.
//  * 
//  */
Matrix::Matrix(int rows, int cols, int dimT, bool initialize):
	_rows(rows),
	_cols(cols),
	_dimT(dimT)
{
	allocateMemory(initialize);
}

/**
 * Copy constructor.
 * 
 * 		Matrix *newm = new Matrix( (*m) );
 * 
 * \param[in] m A pointer to the \a Matrix object to copy from.
 * 
 */
Matrix::Matrix( const Matrix &m )
{
	copyFrom(m);
}

// /**
//  * A special constructor for the class with both the number of rows and columns as parameters. 
//  * Creates a new object.
//  * 
//  * \param[in] r The number of rows.
//  * \param[in] c The number of columns.
//  * 
//  */
// Matrix Matrix::redim( int r, int c )
// {
// 	int i, j;
	
// 	try
// 	{
// 		if( r <= 0  ||  c <= 0 )							// bounds checking
// 			throw IDException( "Matrix constructor has wrong size.", 35 );
// 		for( i = 0; i < _rows; i++ )						// deallocates the old object
// 			delete [] _data[ i ];
// 		delete [] _data;
// 		setdim( r, c );										// number of rows, columns
// 		setmaxdim( r, c );									// max. number of rows, columns
// 		_data = new double*[ r ];								// data allocation
// 		if( (_data == 0)  ||  (_data == NULL))
// 			throw IDException( "Memory allocation failure.", 3 );
// 		for( i = 0; i < r; i++ )
// 		{
// 			_data[ i ] = new double[ c ];
// 			if( (_data[ i ] == 0)  ||  (_data[ i ] == NULL) )
// 				throw IDException( "Memory allocation failure.", 3 );
// 		}
// 	}
// 	catch( bad_alloc e )
// 	{
// 		printf( red );
// 		printf( "\n***Exception bad_alloc found:");
//         printf( "\n***%s" , e.what() );
// 		printf( normal );
//         throw 4;
// 	}
// 	catch( IDException e )
// 	{
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	catch(...)												// other exceptions
// 	{
// 		IDException e( "Error when allocating a matrix.", 34 );
//         e.report();
// 		throw e.getErrCode();
// 	}
	
// 	return *this;
// }

// /**
//  * A special constructor for the class with both the number of rows and columns as parameters,
//  * as well as the dimension of the elements' type (e.g. the Taylor polynomial's grade). 
//  * Creates a new object.
//  * 
//  * \param[in] r The number of rows.
//  * \param[in] c The number of columns.
//  * \param[in] dimT The dimension of the type \type T.
//  * 
//  */
// Matrix Matrix::redim( int r, int c, int dimT )
// {
// 	int i, j;
	
// 	try
// 	{
// 		if( r <= 0  ||  c <= 0 )							// bounds checking
// 			throw IDException( "Matrix constructor has wrong size.", 35 );
// 		for( i = 0; i < _rows; i++ )						// deallocates the old object
// 			delete [] _data[ i ];
// 		delete [] _data;
// 		setdim( r, c, dimT );								// number of rows, columns, dim. of T
// 		setmaxdim( r, c );									// max. number of rows, columns
// 		_data = new double*[ r ];								// data allocation
// 		if( (_data == 0)  ||  (_data == NULL))
// 			throw IDException( "Memory allocation failure.", 3 );
// 		for( i = 0; i < r; i++ )
// 		{
// 			_data[ i ] = new double[ c ];
// 			if( (_data[ i ] == 0)  ||  (_data[ i ] == NULL) )
// 				throw IDException( "Memory allocation failure.", 3 );
// 			for( int j = 0; j < c; j++ )
// 				_data[ i ][ j ] = T( dimT );
// 		}
// 	}
// 	catch( bad_alloc e )
// 	{
// 		printf( red );
// 		printf( "\n***Exception bad_alloc found:");
//         printf( "\n***%s" , e.what() );
// 		printf( normal );
//         throw 4;
// 	}
// 	catch( IDException e )
// 	{
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	catch(...)												// other exceptions
// 	{
// 		IDException e( "Error when allocating a matrix.", 34 );
//         e.report();
// 		throw e.getErrCode();
// 	}
	
// 	return *this;
// }

// //
// // D E S T R U C T O R
// //

// /**
//  * Destructor. Cleans up the object.
//  * 
//  */
// Matrix::~Matrix()
// {
// /*	for( int i = 0; i < _rows; i++ )
// 		delete [] _data[ i ];
// */	delete [] _data;
// }

// //
// // O V E R L O A D E D   ( )   A N D   =   O P E R A T O R S
// //

/**
 * Implements the () operator.
 * 
 * (T & Matrix::... with '&' allows m(2] = 7, i.e., '( )' also in the left side!!)
 * 
 * \param[in] i The row index.
 * \param[in] j The column index.
 * \return The value of the desired element.
 * 
 */
double &Matrix::operator()(int i, int j)
{ 
	// bounds checking
	if( i >= _rows  ||  i < 0  ||  j >= _cols  ||  j < 0 )
	{
		exception("Wrong matrix indexing.", 36).what();
		throw 36;
	}
	return _data[i][j];
}

/**
 * Implements the assignment operator.
 * 
 * 		newm = m;
 * 
 * \param[in] m The \a Matrix object to assign to.
 * \return A pointer to the \a Matrix new object.
 * 
 */
Matrix Matrix::operator=(const Matrix &m)
{
	if( this != &m )
	  	copyFrom(m);
	
	return *this;
}

// //
// // O V E R L O A D E D   A R I T H M E T I C   O P E R A T O R S
// //

/**
 * Implements the == operator. It compares two matrices.
 * 
 * \param[in] m The matrix to compare with.
 * \return \a true if the matrices are equal. Otherwise it returns \a false.
 * 
 */
bool Matrix::operator==( const Matrix &m )
{
	if( _rows != m._rows || _cols != m._cols ) // TODO _dimT ?
		return false;

	for( int i = 0; i < _rows; i++ )
		for( int j = 0; j < _cols; j++ )
			if( _data[ i ][ j ] != m._data[ i ][ j ] )
				return false;

	return true;
}

/**
 * Implements the != operator. It compares two matrices.
 * 
 * \param[in] m The matrix to compare with.
 * \return \a true if the matrices are not equal. Otherwise it returns \a false.
 * 
 */
bool Matrix::operator!=( const Matrix &m )
{
	return !(*this == m);
}

/**
 * Implements the + operator. It adds up two matrices.
 * 
 * \param[in] m The \a Matrix object to add up.
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator+(const Matrix &m)
{	
	// dimensions checking
	if( _rows != m._rows  ||  _cols != m._cols )
	{
		exception("Cannot substract the matrices. Operation not allowed.", 37).what();
		throw 37;
	}

	// an auxiliary object
	Matrix aux(_rows, _cols, _dimT, false);

	for( int i = 0; i < m._rows; i++ )
		for( int j = 0; j < m._cols; j++ )
			aux._data[i][j] = _data[i][j] + m._data[i][j];

	return aux;
}

/**
 * Implements the += operator. It adds up two matrices.
 * 
 * \param[in] m The \a Matrix object to add up to the current matrix.
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator+=( const Matrix &m )
{	
	// dimensions checking
	if( _rows != m._rows || _cols != m._cols )
	{
		// TODO more information?
		exception( "Cannot add up the matrices. Operation not allowed.", 37 ).what();
		throw 37;
	}
	for( int i = 0; i < m._rows; i++ )
		for( int j = 0; j < m._cols; j++ )
			_data[i][j] += m._data[i][j];

 	return *this;
}

/**
 * Implements the - operator. It substracts two matrices.
 * 
 * \param[in] m The \a Matrix object to substract from.
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator-( const Matrix &m )
{	
	// dimensions checking
	if( _rows != m._rows  ||  _cols != m._cols )
	{
		exception("Cannot substract the matrices. Operation not allowed.", 38 ).what();
		throw 38;
	}

	// an auxiliary object
	Matrix aux(_rows, _cols, _dimT, false);

	for( int i = 0; i < m._rows; i++ )
		for( int j = 0; j < m._cols; j++ )
			aux._data[i][j] = _data[i][j] - m._data[i][j];

	return aux;
}

/**
 * Implements the -= operator. It substracts two matrices.
 * 
 * \param[in] m The \a Matrix object to substract from.
 * \return The resulting matrix.
 * 
 */
Matrix Matrix::operator-=( const Matrix &m )
{		
	// dimensions checking
	if( _rows != m._rows || _cols != m._cols )
	{
		// TODO more information?
		exception( "Cannot substract up the matrices. Operation not allowed.", 38 ).what();
		throw 38;
	}
	for( int i = 0; i < m._rows; i++ )
		for( int j = 0; j < m._cols; j++ )
			_data[i][j] -= m._data[i][j];

 	return *this;
}

/**
 * Implements the unary - operator.
 * 
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator-()
{	
	Matrix aux( _rows, _cols, _dimT);
	
	for( int i = 0; i < _rows; i++ )
		for( int j = 0; j < _cols; j++ )
			aux._data[i][j] = - _data[i][j];
	
	return aux;
}

// /**
//  * Implements the * operator. It multiplies a matrix by a scalar.
//  * 
//  * \param[in] alpha The scalar to multiply by.
//  * \return A pointer to the resulting \a Matrix object.
//  * 
//  */
// Matrix Matrix::operator*( double alpha )
// {
// 	Matrix aux( _rows, _cols, _dimT );					// an auxiliary object

// 	for( int i = 0; i < _rows; i++ )
// 		for( int j = 0; j < _cols; j++ )
// 			aux( i, j ) = _data[ i ][ j ] * alpha;
	
// 	return aux;
// }

// //
// // E S P E C I A L   M A T R I X   M U L T I P L I C A T I O N S
// //

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*B + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		B : p-by-n matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  * 
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaABbC( double alpha, double beta, Matrix &A, Matrix &B )
// {
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < B.ncols(); j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 )// A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( int k = 0; k < B.nrows(); k++ )
// 					_data[ i ][ j ] = A( i, k )*B( k, j )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*B + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		B : p-by-n matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  * 
//  * where the inferior-right block of B is an identity matrix like in:
//  * 
//  * 		( * * * 0 0 )
//  * 		( * * * 0 0 )
//  * 		( 0 0 0 1 0 )
//  * 		( 0 0 0 0 1 )
//  * 
//  * so that a particular block multiplication is needed.
//  * 
//  * \param[in] r The number of rows in \a B that are of interest (2 in the example above).
//  * \param[in] c The number of columns in \a B that are of interest (3 in the example above).
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \return The error code.
//  *
//  */
// int Matrix::bmmCaABbC( int r, int c, double alpha, double beta, 
// 		Matrix &A, Matrix &B )
// {
// 	int i, j, k;
	
// 	try
// 	{
// 		for( i = 0; i < A.nrows(); i++ )
// 			for( j = 0; j < c; j++ )						// only first c columns from B 
// 			{												// are interesting
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 )// A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( k = 0; k < r; k++ )					// only first r rows from B (col. from A)
// 					_data[ i ][ j ] = A( i, k )*B( k, j )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		for( i = 0; i < A.nrows(); i++ )
// 			for( j = r; j < B.ncols(); j++ )					// last columns of A remain the same
// 				_data[ i ][ j ] = A( i, j )*alpha;
// 		//printm( "C = alpha*A*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*B + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		B : p-by-n matrix
//  * 		C : r-by-n matrix (only the last r rows from A are interesting)
//  * 		alpha, beta : real numbers
//  * 
//  * where A ("special" A) is of the form:
//  * 
//  * 		(         )
//  * 		(    X    )
//  * 		( ------- )
//  * 		( * ... * )
//  * 		( * ... * )
//  * 
//  * so that a particular matrix multiplication is needed.
//  * 
//  * \param[in] r The last rows in \a A that are of interest (2 non-zero-rows in the example above).
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \return The error code.
//  *
//  */
// int Matrix::mmCasABbC( int r, double alpha, double beta, Matrix &A, Matrix &B )
// {
// 	int i, j, k;
	
// 	try
// 	{
// 		for( int i = A.nrows() - r; i < A.nrows(); i++ )	// only from row A.nrows()-r on
// 			for( int j = 0; j < B.ncols(); j++ )
// 			{
// 				if( strcmp( _data[ i - A.nrows() + r ][ j ].typeName(), "TPolyn" ) == 0 )
// 															// A zero value indicates 
// 															// that both strings are equal
// 					_data[ i - A.nrows() + r ][ j ].set2zero();// p(x) = 0, initialization
// 				//else _data[ i - A.nrows() + r ][ j ] = 0.0;
// 				for( int k = 0; k < B.nrows(); k++ )
// 					_data[ i - A.nrows() + r ][ j ] = A( i, k )*B( k, j )*alpha + 
// 								_data[ i - A.nrows() + r ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*B + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		B : p-by-n matrix
//  * 		C : m-by-r matrix (only the last r columns from B are interesting).
//  * 		alpha, beta : real numbers
//  * 
//  * where B ("special" B) is of the form:
//  * 
//  * 		(   | * * * )
//  * 		(   | . . . )
//  * 		( X | . . . )
//  * 		(   | . . . )
//  * 		(   | * * * )
//  * 
//  * so that a particular matrix multiplication is needed.
//  * 
//  * \param[in] r The last columns in \a B that are of interest (3 non-zero-columns 
//  * 				in the example above).
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \return The error code.
//  *
//  */
// int Matrix::mmCaAsBbC( int r, double alpha, double beta, Matrix &A, Matrix &B )
// {
// 	int i, j, k;
	
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = B.ncols() - r; j < B.ncols(); j++ )// only from column B.ncols()-r on
// 			{
// 				if( strcmp( _data[ i ][ j - B.ncols() + r ].typeName(), "TPolyn" ) == 0 )
// 															// A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j - B.ncols() + r ].set2zero();// p(x) = 0, initialization
// 				//else _data[ i ][ j - B.ncols() + r ] = 0.0;
// 				for( int k = 0; k < B.nrows(); k++ )
// 					_data[ i ][ j - B.ncols() + r ] = A( i, k )*B( k, j )*alpha + 
// 								_data[ i ][ j - B.ncols() + r ]*beta;
// 			}
// 		//printm( "C = alpha*A*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*UTB + beta*C
//  * 
//  * where UTB means that only the upper triangular part is of interest. Furthermore, a column
//  * pivoting on B is considered.
//  * 
//  * with A : m-by-p matrix
//  * 		B : p-by-n matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  * 
//  * E.g.: Q*R
//  * 
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix. Only its upper triangular 
//  * 		part is interesting.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
//  * \return The error code.
//  *
//  */
// int Matrix::mmCaAUTBPbC( double alpha, double beta, Matrix &A, 
// 		Matrix &B, int *piv )
// {
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < B.ncols(); j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( int k = 0; k < B.nrows(); k++ )
// 					if( k <= j )
// 					_data[ i ][ j ] = A( i, k )*B( k, piv[ j ] )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A*UTB + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*A^T + beta*C
//  * 
//  * with A, C : m-by-m matrix
//  * 		alpha, beta : real numbers
//  *
//  * E.g.: Q * Q^T
//  * 
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is also considered.
//  * \return The error code.
//  *
//  */
// int Matrix::mmCaAATbC( double alpha, double beta, Matrix &A )
// {
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < A.nrows(); j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( int k = 0; k < A.nrows(); k++ )
// 					_data[ i ][ j ] = A( i, k )*A( j, k )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A*A^T + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A^T*A + beta*C
//  * 
//  * with A, C : m-by-m matrix
//  * 		alpha, beta : real numbers
//  *
//  * E.g.: Q^T * Q
//  * 
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is also considered.
//  * \return The error code.
//  *
//  */
// int Matrix::mmCaATAbC( double alpha, double beta, Matrix &A )
// {
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < A.nrows(); j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( int k = 0; k < A.nrows(); k++ )
// 					_data[ i ][ j ] = A( k, i )*A( k, j )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A^T*A + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A^T*B + beta*C
//  * 
//  * with A : p-by-m matrix
//  * 		B : p-by-n matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  *
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is considered.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \return The error code.
//  *
//  */
// int Matrix::mmCaATBbC( double alpha, double beta, Matrix &A, Matrix &B )
// {
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < B.ncols(); j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( int k = 0; k < B.nrows(); k++ )
// 					_data[ i ][ j ] = A( k, i )*B( k, j )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A^T*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A^T*B + beta*C
//  * 
//  * with A : p-by-m matrix
//  * 		B : p-by-n matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  *
//  * and a column pivoting on A's rows.
//  * E.g.: Q^T * A, Q^T * OA
//  * 
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is considered.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaATBPbC( double alpha, double beta, Matrix &A, 
// 		Matrix &B, int *piv )
// {
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < B.ncols(); j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( int k = 0; k < B.nrows(); k++ )
// 					_data[ i ][ j ] = A( k, i )*B( k, piv[ j ] )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A^T*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*B^T + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		B : n-by-p matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  *
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix. Its transpose is considered.
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaABTbC( double alpha, double beta, Matrix &A, Matrix &B )
// {
// 	try
// 	{
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < B.ncols(); j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( int k = 0; k < B.nrows(); k++ )
// 					_data[ i ][ j ] = A( i, k )*B( j, k )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		//printm( "C = alpha*A*B^T + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*B^T + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		B : n-by-p matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  *
//  * After transposing B, either its first or its last rows are considered for multiplication,
//  * 
//  * 		( * ... * )				(    X    )
//  * 		( ------- )		or		( ------- )
//  * 		(         )				( * ... * )
//  * 		(    X    )				( * ... * )
//  * 
//  * according to A dimensions. I.e., the matrix A has less columns than B^T rows has.
//  * 
//  * \param[in] r The number of rows from \a B that should be considered.
//  * \param[in] up The binary parameter to indicate whether the first or the last \a r rows 
//  * 		from \a B should be considered (=1 if first rows; =0 otherwise).
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix. Its transpose is considered.
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaABTbC( int r, bool up, double alpha, double beta, 
// 		Matrix &A, Matrix &B )
// {
// 	try
// 	{
// 		//printf( "\nr=%d", r );
// 		//printf( "\nA.nrows()=%d, A.ncols()=%d", A.nrows(), A.ncols() );
// 		//printf( "\nB.nrows()=%d, B.ncols()=%d", B.nrows(), B.ncols() );
// 		//printf( "\nC.nrows()=%d, C.ncols()=%d", nrows(), ncols() );
// 		for( int i = 0; i < A.nrows(); i++ )
// 			for( int j = 0; j < B.ncols(); j++ )
// 				{
// 					if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) 
// 															// A zero value indicates 
// 															// that both strings are equal
// 						_data[ i ][ j ].set2zero();			// p(x) = 0, initialization
// 					//else _data[ i ][ j ] = 0.0;
// 					if( up )								// first r columns from B are used
// 					{
// 						//printf( "\nr=%d", r );
// 						for( int k = 0; k < r; k++ )
// 							_data[ i ][ j ] = A( i, k )*B( j, k )*alpha + _data[ i ][ j ]*beta;
// 					}
// 					else									// last r columns from B are used
// 						for( int k = B.nrows() - r; k < B.nrows(); k++ )
// 						_data[ i ][ j ] = A( i, k - B.nrows() + r  )*B( j, k )*alpha 
// 							+ _data[ i ][ j ]*beta;
// 		}
// 		//printm( "C = alpha*A*B^T + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*B^T + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		B : n-by-p matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  *
//  * where the inferior-right block of A is an identity matrix like in:
//  * 
//  * 		( * * * 0 0 )
//  * 		( * * * 0 0 )
//  * 		( 0 0 0 1 0 )
//  * 		( 0 0 0 0 1 )
//  * 
//  * so that a particular block multiplication is needed.
//  * 
//  * \param[in] r The number of rows in \a A that are of interest (2 in the example above).
//  * \param[in] c The number of columns in \a A that are of interest (3 in the example above).
//  * \param[in] alpha The scalar value that multiplies \a A*B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix. Its transpose is considered.
//  * \return The error code.
//  * 
//  */
// int Matrix::bmmCaABTbC( int r, int c, double alpha, double beta, 
// 		Matrix &A, Matrix &B )
// {
// 	int i, j, k;

// 	try
// 	{
// 		for( i = 0; i < r; i++ )							// only first r rows from A interesting
// 			for( j = 0; j < c; j++ )						// only first c col. from A (col. from B)
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 				for( k = 0; k < B.nrows(); k++ )
// 					_data[ i ][ j ] = A( i, k )*B( j, k )*alpha + _data[ i ][ j ]*beta;
// 			}
// 		for( i = r; i < A.nrows(); i++ )					// last rows of B remain the same
// 			for( j = 0; j < B.ncols(); j++ )
// 				_data[ i ][ j ] = B( j, i )*alpha;
// 		//printm( "C = alpha*A*B^T + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*I*B + beta*C
//  * 
//  * with I : m-by-p matrix; identity matrix
//  * 		B : p-by-n matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  * 
//  * \param[in] alpha The scalar value that multiplies \a B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaIBbC( double alpha, double beta, Matrix &B )
// {
// 	int i, j;
	
// 	try
// 	{
// 		for( i = 0; i < nrows(); i++ )						// initialize C
// 			for( j = 0; j < B.ncols(); j++ )
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 		if( nrows() >= B.nrows() )							// last rows from I are zeroed
// 		{
// 			for( i = 0; i < B.nrows(); i++ )				// last rows from C are already zeroed!!
// 				for( j = 0; j < B.ncols(); j++ )
// 					_data[ i ][ j ] = B( i, j )*alpha + _data[ i ][ j ]*beta;			
// 		}
// 		else												// last columns from I are zeroed
// 		{
// 			for( i = 0; i < nrows(); i++ )					// C has no more row!!
// 				for( j = 0; j < B.nrows(); j++ )
// 					_data[ i ][ j ] = B( i, j )*alpha + _data[ i ][ j ]*beta;			
// 		}
// 		//printm( "C = alpha*I*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*I*B + beta*C
//  * 
//  * with I : m-by-p matrix; identity matrix permuted according to a vector of permutations, piv
//  * 		B : p-by-n matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  * 
//  * \param[in] alpha The scalar value that multiplies \a B.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on \a I, of type \type int.
//  * \param[in] rows The binary parameter to indicate whether the rows or the columns of I
//  * 		should be permuted (= 0, the rows; = 1, the columns).
//  * \param[in] B The pointer to \a B, an object of type \type Matrix.
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaIBbC( double alpha, double beta, int *piv, bool rows, Matrix &B )
// {
// 	int i, j;
	
// 	try
// 	{
// 		Matrix Id( nrows(), B.nrows(), dimT() );			// Id[m][p][nTcoeff]
// 		Id.set2Id();										// set to the identity
// 		if( rows )											// permute its columns
// 			Id.cpermutem( piv );
// 		else												// permute its rows
// 			Id.rpermutem( piv );
// 		mmCaABbC( alpha, beta, Id, B );
// 		//printm( "C = alpha*I*B + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*I + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		I : p-by-n matrix; identity matrix
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  * 
//  * \param[in] alpha The scalar value that multiplies \a A.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaAIbC( double alpha, double beta, Matrix &A )
// {
// 	int i, j;
	
// 	try
// 	{
// 		for( i = 0; i < A.nrows(); i++ )					// initialize C
// 			for( j = 0; j < ncols(); j++ )
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0, initialization
// 				//else _data[ i ][ j ] = 0.0;
// 		if( A.nrows() <= nrows() )							// last columns from I are zeroed
// 		{
// 			for( i = 0; i < A.nrows(); i++ )				// C has no more row!!
// 				for( j = 0; j < nrows(); j++ )
// 					_data[ i ][ j ] = A( i, j )*alpha + _data[ i ][ j ]*beta;			
// 		}
// 		else												// last rows from I are zeroed
// 		{
// 			for( i = 0; i < A.nrows(); i++ )				// C has no more column!!
// 				for( j = 0; j < ncols(); j++ )
// 					_data[ i ][ j ] = A( i, j )*alpha + _data[ i ][ j ]*beta;			
// 		}
// 		//printm( "C = alpha*A*I + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Matrix multiplication of the form:
//  * 
//  *     C = alpha*A*I + beta*C
//  * 
//  * with A : m-by-p matrix
//  * 		I : p-by-n matrix; identity matrix permuted according to a vector of permutations, piv
//  * 		C : m-by-n matrix
//  * 		alpha, beta : real numbers
//  * 
//  * \param[in] alpha The scalar value that multiplies \a A.
//  * \param[in] beta The scalar value that multiplies \a C.
//  * \param[in] A The pointer to \a A, an object of type \type Matrix.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on \a I, of type \type int.
//  * \param[in] rows The binary parameter to indicate whether the rows or the columns of I
//  * 		should be permuted (= 0, the rows; = 1, the columns).
//  * \return The error code.
//  * 
//  */
// int Matrix::mmCaAIbC( double alpha, double beta, Matrix &A, int *piv, bool rows )
// {
// 	int i, j;
	
// 	try
// 	{
// 		Matrix Id( A.ncols(), ncols(), dimT() );			// Id[p][n][nTcoeff]
// 		Id.set2Id();										// set to the identity
// 		if( rows )											// permute its columns
// 			Id.cpermutem( piv );
// 		else												// permute its rows
// 			Id.rpermutem( piv );
// 		mmCaABbC( alpha, beta, A, Id );
// 		//printm( "C = alpha*A*I + beta*C = \n" );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Error in matrix multiplication. The matrices dimensions are probably wrong.", 10 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// //
// // S O L V I N G   S Y S T E M S   O F   E Q U A T I O N S
// //

// /**
//  * Solves the equation
//  * 
//  * 		U X = B
//  * 
//  * by back-substitution, where:
//  * 
//  * 		U : m-by-m upper triangular matrix, non-singular
//  * 		X : m-by-n matrix
//  * 		B : m-by-n matrix, overwritten with the solution on output.
//  *
//  * The X_ik are calculated making few modifications to the algorithm 3.1.2, p.89 from 
//  * Golub & Van Loan's book:
//  * 
//  * 		x_ik = ( b_ik - sum_{j=i+1}^{m} u_ij * x_jk ) / u_ii     for k=1,...,n
//  * 
//  * \param[in/out] B The pointer to \a B, an object of type \type Matrix that is the independent 
//  * 		term on input. On output it contains the calculated solution \a X.
//  * \return The error code.
//  * 
//  */
// int Matrix::utsolve( Matrix &B )
// {
// 	T sum( dimT() );										// auxiliary variable
	
// 	for( int k = 0; k < B.ncols(); k++ )
// 	{
// 		B( _rows - 1, k ) /= _data[ _rows - 1 ][ _rows - 1 ];
// 		for( int i = _rows - 2; i > -1; i-- )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0, initialization
// 			//else sum = 0.0;
// 			for( int j = i + 1; j < _rows; j++ )
// 				sum += _data[ i ][ j ] * B( j, k );
// 			B( i, k ) = ( B( i, k ) - sum ) / _data[ i ][ i ];
// 		}
// 	}
// 	return 0;
// }

// /**
//  * Solves the equation
//  * 
//  * 		U X = B
//  * 
//  * by back-substitution, where:
//  * 
//  * 		U : m-by-m upper triangular matrix, non-singular
//  * 		X : m-by-n matrix
//  * 		B : m-by-n matrix, overwritten with the solution on output.
//  *
//  * The X_ik are calculated making few modifications to the algorithm 3.1.2, p.89 from 
//  * Golub & Van Loan's book:
//  * 
//  * 		x_ik = ( b_ik - sum_{j=i+1}^{m} u_ij * x_jk ) / u_ii     for k=1,...,n
//  * 
//  * \param[in/out] B The pointer to \a B, an object of type \type Matrix that is the independent 
//  * 		term on input. On output it contains the calculated solution \a X.
//  * \return The error code.
//  * 
//  */
// int Matrix::utsolve( Matrix &B, Matrix &X, int *piv )
// {
// 	T sum( dimT() );
	
// 	for( int k = 0; k < B.ncols(); k++ )
// 	{
// 		X( _rows - 1, k ) = B( piv[ _rows - 1 ], k ) / _data[ piv[ _rows - 1 ] ][ _rows - 1 ];
// 		for( int i = _rows - 2; i > -1; i-- )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0, initialization
// 			//else sum = 0.0;
// 			for( int j = i + 1; j < _rows; j++ )
// 				sum += _data[ piv[ i ] ][ j ] * X( j, k );
// 			X( i , k ) = ( B( piv[ i ], k ) - sum ) / _data[ piv[ i ] ][ i ];
// 		}
// 	}
// 	return 0;
// }

// /**
//  * Solves the equation
//  * 
//  * 		U x = b
//  * 
//  * by back-substitution, where:
//  * 
//  * 		U : n-by-n upper triangular matrix, non-singular
//  * 		x : n vector
//  * 		B : n vector, overwritten with the solution on output.
//  *
//  * The x_i are calculated following the algorithm 3.1.2, p.89 from Golub & Van Loan's book:
//  * 
//  * 		x_i = ( b_i - sum_{j=i+1}^{n} u_ij * x_j ) / u_ii
//  * 
//  * \param[in/out] b The pointer to \a b, an object of type \type Matrix that is the independent 
//  * 		term on input. On output it contains the calculated solution \a x.
//  * \return The error code.
//  * 
//  */
// int Matrix::utsolve( T *b )
// {
// 	T sum( dimT() );
	
// 	//printf( red );
// 	b[ _rows - 1 ] /= _data[ _rows - 1 ][ _rows - 1 ];
// 	//printf( "\nb[ n - 1 ]=%.16g", b[ n - 1 ] );
// 	for( int i = _rows - 2; i > -1; i-- )
// 	{
// 		if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 		// A zero value indicates 
// 															// that both strings are equal
// 			sum.set2zero();									// p(x) = 0, initialization
// 		//else sum = 0.0;
// 		for( int j = i + 1; j < _rows; j++ )
// 			sum += _data[ i ][ j ] * b[ j ];
// 		//printf( "\nsum=%.16g", sum );
// 		b[ i ] = ( b[ i ] - sum ) / _data[ i ][ i ];
// 		//printf( "\ni=%d  -->  ", i );
// 		//printf( "b[ i ]=%.16g", b[ i ] );
// 	}
// 	//printf( "\n" );
// 	//printf( normal );
// 	return 0;
// }

// /**
//  * Solves the equation
//  * 
//  * 		X U = B
//  * 
//  * by back-substitution, where:
//  * 
//  * 		U : m-by-m upper triangular matrix, non-singular
//  * 		X : m-by-n matrix
//  * 		B : m-by-n matrix, overwritten with the solution on output.
//  *
//  * The X_ik are calculated making few modifications to the function 'utsolve' for UX=B.
//  * 
//  * \param[in/out] B The pointer to \a B, an object of type \type Matrix that is the independent 
//  * 		term on input. On output it contains the calculated solution \a X.
//  * \return The error code.
//  * 
//  */
// int Matrix::utxsolve( Matrix &B )
// {
// 	T sum( dimT() );
	
// 	//printm( "A = \n", yellow );
// 	//B.printm( "B = \n", yellow );
// 	for( int k = 0; k < B.nrows(); k++ )
// 	{
// 		B( k, 0 ) /= _data[ 0 ][ 0 ];
// 		for( int i = 1; i < B.ncols(); i++ )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0, initialization
// 			//else sum = 0.0;
// 			for( int j = 0; j < i; j++ )
// 				sum += _data[ j ][ i ] * B( k, j );
// 			B( k, i ) = ( B( k, i ) - sum ) / _data[ i ][ i ];
// 		}
// 	}
// 	return 0;
// }

// /**
//  * Solves the equation
//  * 
//  * 		A X = B
//  * 
//  * by Gaussian elimination using scalling and column pivoting, where:
//  * 
//  * 		A : m-by-m matrix
//  * 		X : m-by-n matrix
//  * 		B : m-by-n matrix, overwritten with the solution on output.
//  *
//  * The algorithm from A. L. Garcia, Numerical Methods for Physics (Prentice Hall, Second edition, 
//  * Englewood Cliffs NJ, 2000) is used.
//  * 
//  * \param[in/out] B The pointer to \a B, an object of type \type Matrix that is the independent 
//  * 		term on input. On output it contains the calculated solution \a X.
//  * \return The error code.
//  * 
//  */
// int Matrix::gsolve( Matrix &B )
// {
// 	int h, i, j, k, jpiv, jindex, *index;
// 	double scalemax = 0.0, ratiomax, ratio, *scale;
// 	T coeff( dimT() );
// 	Matrix X( B.nrows(), B.ncols(), B.dimT() );
  
// 	try
// 	{	
// 		scale = new double[ _rows ];						// scale[m], scale factor
// 		index = new int[ _rows ];							// index[m], row index list
// 		if( (scale == 0)  ||  (scale == NULL)  ||
// 			(index == 0)  ||  (index == NULL) )
// 			throw IDException( "Memory allocation failure.", 3 );
// 		for( i = 0; i < _rows; i++ )						// determine scale factor, for each row
// 		{
// 			index[ i ] = i;									// initialize row index list
// 			for( j = 0; j < _rows; j++ )
// 				if( scalemax < fabs( _data[ i ][ j ].feval() ) )// scale[i] = max(|A[i][j][0]|)
// 					scalemax = fabs( _data[ i ][ j ].feval() );
// 			scale[ i ] = scalemax;
// 		}
// 		for( k = 0; k < _rows - 1; k++ )					// select pivot row...
// 		{													// ...from max(|A[i][k]/scale[i]|)
// 			jpiv = k;
// 			ratiomax = 0.0;
// 			for( i = k; i < _rows; i++ )
// 			{
// 				ratio = fabs( _data[ index[ i ] ][ k ].feval() / scale[ index[ i ] ]);
// 				if( ratio > ratiomax )
// 				{
// 					jpiv = i;
// 					ratiomax = ratio;
// 				}
// 			}
// 			jindex = index[ k ];
// 			if( jpiv != k )									// perform pivoting
// 			{
// 				jindex = index[ jpiv ];
// 				index[ jpiv ] = index[ k ];					// swap index jpiv and k
// 				index[ k ] = jindex;
// 			}
// 			for( i = k + 1; i < _rows; i++ )				// forward elimination
// 			{
// 				coeff = _data[ index[ i ] ][ k ] / _data[ jindex ][ k ];
// 				for( j = k + 1; j < _rows; j++ )
// 					_data[ index[ i ] ][ j ] -= coeff * _data[ jindex ][ j ];
// 				_data[ index[ i ] ][ k ] = coeff;
// 				for( h = 0; h < B.ncols(); h++ )			// for all columns in B
// 				{
// 					B( index[ i ], h ) -= coeff * B( jindex, h );
// 				}
// 			}
// 		}
// 		//printf( "\nGaussian elimination Vi * Wihat = b\n" );
// 		//printiv( _rows, index, "index vector", red );
// 		//double eps = 0.1E-16;
// 		//printtpm( "\nVi overwritten = \n", cyan, eps );
// 		//B.printtpm( "\nb overwritten = \n", cyan, eps );
// 		//printf( "\nBack-substitution..." );
// 		utsolve( B, X, index );								// back-substitution
// 		//B.mmCaABbC( 1.0, 1.0, *this, X );
// 		//B.printtpm( "\nVi * Wihat = \n", cyan, eps );
// 		//X.printtpm( "\nWihat = \n", cyan, eps );
// 		B = X;												// overwrite B: output parameter
// 	} //try
// 	catch( bad_alloc e )
// 	{
// 		delete [] scale;
// 		delete [] index;
// 		printf( red );
// 		printf( "\n***Exception bad_alloc found:");
//         printf( "\n***%s" , e.what() );
// 		printf( normal );
//         throw 4;
// 	}
// 	catch( out_of_range e )
// 	{
// 		delete [] scale;
// 		delete [] index;
// 		printf( red );
// 		printf( "\n***Exception out_of_range found:");
//         printf( "\n***%s" , e.what() );
// 		printf( normal );
//         throw 5;
// 	}
// 	catch( IDException e )
// 	{
// 		delete [] scale;
// 		delete [] index;
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	delete [] scale;
// 	delete [] index;
// 	return 0;        
// }

// /**
//  * Solves the equation
//  * 
//  * 		A x = b
//  * 
//  * by Gaussian elimination using scalling and column pivoting, where:
//  * 
//  * 		A : n-by-n matrix
//  * 		x : n vector
//  * 		b : n vector, overwritten with the solution on output.
//  *
//  * The algorithm from A. L. Garcia, Numerical Methods for Physics (Prentice Hall, Second edition, 
//  * Englewood Cliffs NJ, 2000) is used.
//  * 
//  * \param[in/out] b The pointer to \a b, a vector of type \type T that is the independent 
//  * 		term on input. On output it contains the calculated solution \a x.
//  * \return The error code.
//  * 
//  */
// int Matrix::gsolve( T *b )
// {
// 	int i, j, k, jpiv, jindex, *index;
// 	double scalemax = 0.0, ratiomax = 0.0, ratio, *scale;
// 	T coeff( dimT() );
  
// 	try
// 	{	
// 		scale = new double[ _rows ];						// scale[n], scale factor
// 		index = new int[ _rows ];							// index[n], row index list
// 		if( (scale == 0)  ||  (scale == NULL)  ||
// 			(index == 0)  ||  (index == NULL) )
// 			throw IDException( "Memory allocation failure.", 3 );
// 		for( i = 0; i < _rows; i++ )						// determine scale factor, for each row
// 		{
// 			index[ i ] = i;									// initialize row index list
// 			for( j = 0; j < _rows; j++ )
// 				if( scalemax < fabs( _data[ i ][ j ].feval() ) )// scale[i] = max(|A[i][j][0]|)
// 					scalemax = fabs( _data[ i ][ j ].feval() );
// 			scale[ i ] = scalemax;
// 		}
// 		for( k = 0; k < _rows - 1; k++ )					// select pivot row...
// 		{													// ...from max(|A[i][k]/scale[i]|)
// 			jpiv = k;
// 			for( i = k; i < _rows; i++ )
// 			{
// 				ratio = fabs( _data[ index[ i ] ][ k ].feval() / scale[ index[ i ] ]);
// 				if( ratio > ratiomax )
// 				{
// 					jpiv = i;
// 					ratiomax = ratio;
// 				}
// 			}
// 			jindex = index[ k ];
// 			if( jpiv != k )									// perform pivoting
// 			{
// 				jindex = index[ jpiv ];
// 				index[ jpiv ] = index[ k ];					// swap index jpiv and k
// 				index[ k ] = jindex;
// 			}
// 			for( i = k + 1; i < _rows; i++ )					// forward elimination
// 			{
// 				coeff = _data[ index[ i ] ][ k ] / _data[ jindex ][ k ];
// 				for( j = k + 1; j < _rows; j++ )
// 					_data[ index[ i ] ][ j ] -= coeff * _data[ jindex ][ j ];
// 				_data[ index[ i ] ][ k ] = coeff;
// 				b[ index[ i ] ] -= _data[ index[ i ] ][ k ] * b[ jindex ];
// 			}
// 		}
// 		//printf( "\nGaussian elimination:\n" );
// 		//printm( "A = \n", yellow );
// 		//printdv( _rows, b, "b = \n", yellow );
// 		utsolve( b );										// back-substitution
// 		//printdv( _rows, b, "x = \n", yellow );
// 	} //try
// 	catch( bad_alloc e )
// 	{
// 		delete [] scale;
// 		delete [] index;
// 		printf( red );
// 		printf( "\n***Exception bad_alloc found:");
//         printf( "\n***%s" , e.what() );
// 		printf( normal );
//         throw 4;
// 	}
// 	catch( out_of_range e )
// 	{
// 		delete [] scale;
// 		delete [] index;
// 		printf( red );
// 		printf( "\n***Exception out_of_range found:");
//         printf( "\n***%s" , e.what() );
// 		printf( normal );
//         throw 5;
// 	}
// 	catch( IDException e )
// 	{
// 		delete [] scale;
// 		delete [] index;
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	delete [] scale;
// 	delete [] index;
// 	return 0;        
// }

// //
// // O T H E R   F U N C T I O N S
// //

// /**
//  * Calculates the Frobenius norm of a matrix:
//  * 
//  * 		|| A ||_F = sqrt( sum_i sum_j  a_{ij}^2 )
//  * 
//  * If the matrix is a matrix of Taylor polynomials, then only the first coefficients
//  * of those polynomials are interesting, i.e., the coefficients that correspond to the
//  * evaluation of the function and not to its derivatives (i.e. the resting coefficients
//  * of the Taylor polynomials).
//  * 
//  * \return The resulting value.
//  * 
//  */
// double Matrix::fnorm()
// {
// 	double sum = 0.0;
	
// 	for( int i = 0; i < _rows; i++ )
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 				sum += _data[ i ][ j ].feval() * _data[ i ][ j ].feval();
// 			//else sum += _data[ i ][ j ] * _data[ i ][ j ];	// e.g. a matrix of doubles
// 		}
			
// 	return sqrt( sum );
// }

// /**
//  * Calculates the Frobenius norm of a matrix:
//  * 
//  * 		|| A ||_F = sqrt( sum_i sum_j  a_{ij}^2 )
//  * 
//  * where the coefficients are Taylor polynomials.
//  * 
//  * The Frobenius norm in this case will be a vector, since all Taylor coefficients of 
//  * each polynomial are considered. In other words, the coefficients that correspond to the 
//  * evaluation of the function and the ones of its derivatives (i.e. the resting 
//  * coefficients of the Taylor polynomial) are taken into account.
//  * 
//  * \param[in] nrTC The number of Taylor coefficients to be considered in a polynomial.
//  * \param[out] fn The pointer to \a fn, a vector of type \type double, containing 
//  * 		the resulting values, i.e., Frobenius norms.
//  * \return The error code.
//  * 
//  */
// int Matrix::tcfnorm( int nrTC, double *fn )
// {
// 	double sum;

// 	for( int k = 0; k < nrTC; k++ )
// 	{
// 		sum = 0.0;
// 		for( int i = 0; i < _rows; i++ )
// 			for( int j = 0; j < _cols; j++ )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 																// that both strings are equal
// 					sum += _data[ i ][ j ][ k ] * _data[ i ][ j ][ k ];
// 			}
// 		fn[ k ] = sqrt( sum );
// 	}
			
// 	return 0;
// }

// /**
//  * Calculates the columns norm of a matrix.
//  * 
//  * \param[out] c The pointer to \a c, a vector of type \type double, containing 
//  * 		the resulting values.
//  * \return The error code.
//  * 
//  */
// int Matrix::colnorm( double *c )
// {
// 	try
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			c[ j ] = 0.0;
// 			for( int i = 0; i < _rows; i++ )				// column's norm, sum over the rows
// 				c[ j ] += _data[ i ][ j ].feval() * _data[ i ][ j ].feval();
// 		}
// 		//printdv( _cols, c, "c (column norms) = \n", cyan );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Inconsistent matrix and vector dimensions in parameters.", 11 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Downdates the columns norm.
//  * 
//  * \param[in] pos The starting position.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A (and \a c).
//  * \param[out] c The pointer to \a c, a vector of type \type double, containing the columns norm.
//  * \return The error code.
//  * 
//  */
// int Matrix::colnormdown( int pos, int *piv, double *c )
// {
// 	try
// 	{
// 		for( int i = pos + 1; i < _cols; i++ )
// 			c[ piv[ i ] ] -= _data[ pos ][ piv[ i ] ].feval() * _data[ pos ][ piv[ i ] ].feval();
// 		//printdv( n, c, "c (norm downdate) = \n", cyan );
// 	}
// 	catch(...)
// 	{
// 		IDException e( "Inconsistent matrix and vector dimensions in parameters.", 11 );
//         e.report();
// 		throw e.getErrCode();
// 	}
// 	return 0;
// }

// /**
//  * Permutes the columns of a matrix given a vector of permutations. 
//  * 
//  * For example, in case a matrix A is permuted after a QR decomposition with column pivoting,
//  * then the resulting matrix in the upper triangular matrix R. 
//  * 
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
//  * \return The error code.
//  * 
//  */
// int Matrix::cpermutem( int *piv )
// {
// 	int i, j;
// 	Matrix aux( _rows, _cols, _dimT );					// auxiliary matrix
	
// 	for( i = 0; i < _rows; i++ )
// 		for( j = 0; j < _cols; j++ )
// 			aux( i, j ) = _data[ i ][ piv[ j ] ];
// 	for( i = 0; i < _rows; i++ )
// 		for( j = 0; j < _cols; j++ )
// 			_data[ i ][ j ] = aux( i, j );
// 	return 0;
// }

// /**
//  * Permutes the columns of a matrix given a vector of permutations. 
//  * 
//  * For example, in case a matrix A is permuted after a QR decomposition with column pivoting,
//  * then the resulting matrix in the upper triangular matrix R. 
//  * 
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
//  * \param[in] trans The \a boolean parameter to indicate whether to transpose the vector
//  * 		of permutations \a piv or not (=1, transpose; =0, otherwise).
//  * \return The error code.
//  * 
//  */
// int Matrix::cpermutem( int *piv, bool trans )
// {
// 	int i, j;
// 	Matrix aux( _rows, _cols, _dimT );						// auxiliary matrix
// 	int *pivT = new int[ _cols ];								// auxiliary vector
	
// 	if( trans )
// 		for( i = 0; i < _cols; i++ )							// compute piv^T
// 			pivT[ piv[ i ] ] = i;
// 	for( i = 0; i < _rows; i++ )
// 		for( j = 0; j < _cols; j++ )
// 		{
// 			if( trans )
// 				aux( i, j ) = _data[ i ][ pivT[ j ] ];
// 			else
// 				aux( i, j ) = _data[ i ][ piv[ j ] ];
// 		}
// 	for( i = 0; i < _rows; i++ )
// 		for( j = 0; j < _cols; j++ )
// 			_data[ i ][ j ] = aux( i, j );
// 	return 0;
// }

// /**
//  * Permutes the rows of a matrix given a vector of permutations. 
//  * 
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the rows of \a A.
//  * \return The error code.
//  * 
//  */
// int Matrix::rpermutem( int *piv )
// {
// 	int i, j;
// 	Matrix aux( _rows, _cols, _dimT );						// auxiliary matrix
	
// 	for( i = 0; i < _rows; i++ )
// 		for( j = 0; j < _cols; j++ )
// 			//aux( i, j ) = _data[ piv[ i ] ][ j ];
// 			aux( piv[ i ], j ) = _data[ i ][ j ];
// 	for( i = 0; i < _rows; i++ )
// 		for( j = 0; j < _cols; j++ )
// 			_data[ i ][ j ] = aux( i, j );
// 	return 0;
// }

// /**
//  * Transposes a given matrix. 
//  * 
//  * \return The error code.
//  * 
//  */
// Matrix Matrix::transpm()
// {
// 	Matrix aux( _cols, _rows, _dimT );						// an auxiliary object

// 	for( int i = 0; i < _rows; i++ )
// 		for( int j = 0; j < _cols; j++ )
// 			aux( j, i ) = _data[ i ][ j ];
// 	return aux;
// }

// /**
//  * Implements the shift operator to calculate the derivative of Taylor polynomials
//  * in case the elements of the matrix are such, like in:
//  * 
//  * y(t) = sum_{j=0}^{d} y_j * t^j + O(t^d+1)
//  * 		= y_0 + y_1*t + y_2*t^2 + ... + y_d*t^d
//  * 
//  * y'(t) = y_1 + 2*y_2*t + 3*y_3*t^2 + ... + d*y_d*t^d-1
//  * 
//  * Internally, the coefficients are shifted to the left and the last one is zeroed.
//  * 
//  * \return The error code.
//  * 
//  */
// int Matrix::shift( Matrix &M )
// {	
// 	try
// 	{
// 		if( _rows != M._rows  ||  _cols != M._cols  || _dimT != M._dimT )// dimensions checking
// 			throw IDException( "Cannot shift the matrix.", 39 );
// 		if( strcmp( _data[ 0 ][ 0 ].typeName(), "TPolyn" ) != 0  ||
// 			strcmp( M._data[ 0 ][ 0 ].typeName(), "TPolyn" ) != 0 )// both are matrices of...
// 															// ...Taylor polynomials
// 			throw IDException( "It is not a matrix of Taylor polynomials.", 40 );
// 		M = *this;
// 		for( int i = 0; i < _rows; i++ )
// 			for( int j = 0; j < _cols; j++ )
// 				M._data[ i ][ j ].shift();
// 	}
// 	catch( IDException e )
// 	{
//         e.report();
// 		throw e.getErrCode();
// 	}

// 	return 0;
// }

// /**
//  * Returns \a true in case the given matrix is the identity matrix; \a false otherwise.
//  * 
//  * \return \a true if the matrix is the identity matrix; \a false otherwise.
//  * 
//  */
// bool Matrix::isId()
// {
// 	bool id = true;
	
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _rows; j++ )
// 			if( i == j )
// 				if( !_data[ j ][ i ].isId() )
// 				{
// 					id = false;
// 					break;
// 				} else;
// 			else
// 				if( !_data[ j ][ i ].isZero() )
// 				{
// 					id = false;
// 					break;
// 				}
// 		if( !id )
// 			break;
// 	}
// 	return id;
// }

// /**
//  * Returns \a true in case the given matrix is near the identity matrix; \a false otherwise.
//  * 
//  * \param[in] eps The threshold value for comparing. 
//  * \return \a true if the matrix is the identity matrix; \a false otherwise.
//  * 
//  */
// bool Matrix::isId( double eps )
// {
// 	bool id = true;
	
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _rows; j++ )
// 			if( i == j )
// 				if( !_data[ j ][ i ].isId( eps ) )
// 				{
// 					id = false;
// 					break;
// 				} else;
// 			else
// 				if( !_data[ j ][ i ].isZero( eps ) )
// 				{
// 					id = false;
// 					break;
// 				}
// 		if( !id )
// 			break;
// 	}
// 	return id;
// }

// /**
//  * Returns \a true in case the given submatrix is near the identity matrix; \a false otherwise.
//  * 
//  * \param[in] eps The threshold value for comparing. 
//  * \return \a true if the matrix is the identity matrix; \a false otherwise.
//  * 
//  */
// bool Matrix::isId( int m1, int m2, int n1, int n2, double eps )
// {
// 	bool id = true;
	
// 	for( int i = m1; i < m2; i++ )
// 	{
// 		for( int j = n1; j < n2; j++ )
// 			if( i == j )
// 				if( !_data[ j ][ i ].isId( eps ) )
// 				{
// 					id = false;
// 					break;
// 				} else;
// 			else
// 				if( !_data[ j ][ i ].isZero( eps ) )
// 				{
// 					id = false;
// 					break;
// 				}
// 		if( !id )
// 			break;
// 	}
// 	return id;
// }

// /**
//  * Returns \a true in case the given matrix is the zero matrix; \a false otherwise.
//  * 
//  * \return \a true if the matrix is the zero matrix; \a false otherwise.
//  * 
//  */
// bool Matrix::isZero()
// {
// 	bool z = true;

// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 			if( !_data[ i ][ j ].isZero() )					// all coeff. from p(x) = 0
// 			{
// 				z = false;
// 				break;
// 			}
// 		if( !z )
// 			break;
// 	}
// 	return z;
// }

// /**
//  * Returns \a true in case the given matrix is near the zero matrix; \a false otherwise.
//  * 
//  * \param[in] eps The threshold value for comparing. 
//  * \return \a true if the matrix is the zero matrix; \a false otherwise.
//  * 
//  */
// bool Matrix::isZero( double eps )
// {
// 	bool z = true;

// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 			if( !_data[ i ][ j ].isZero( eps ) )			// at least for one coeff. |p(x)| > eps
// 			{
// 				z = false;
// 				break;
// 			}
// 		if( !z )
// 			break;
// 	}
// 	return z;
// }

// /**
//  * Sets a matrix to the identity one:
//  * 
//  * 		M = I
//  * 
//  * \return The error code.
//  * 
//  */
// int Matrix::set2Id()
// {
// 	for( int i = 0; i < _rows; i++ )
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			if( i == j )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 )// A zero value indicates 
// 															// that both strings are equal
// 				{
// 					_data[ i ][ j ].set2const( 1.0 );		// p(x) = 1
// 				}
// 				//else _data[ i ][ j ] = 1.0;
// 			}
// 			else
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0
// 				//else _data[ i ][ j ] = 0.0;
// 		}
// 	return 0;
// }

// /**
//  * Sets the first submatrix to the identity one:
//  * 
//  * 		e.g. M = ( 1 0 0 |   )
//  * 				 ( 0 1 0 | X )
//  * 				 ( 0 0 1 |   )
//  * 				 (    Y  |   )
//  * 
//  * \param[in] m The number of rows that should be considered.
//  * \param[in] n The number of columns that should be considered.
//  * \return The error code.
//  * 
//  */
// int Matrix::set2Id( int m, int n )
// {
// 	for( int i = 0; i < m; i++ )
// 		for( int j = 0; j < n; j++ )
// 		{
// 			if( i == j )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2const( 1.0 );		// p(x) = 1
// 				//else _data[ i ][ j ] = 1.0;
// 			}
// 			else
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0
// 				//else _data[ i ][ j ] = 0.0;
// 		}
// 	return 0;
// }

// /**
//  * Sets a submatrix to the identity one:
//  * 
//  * 		e.g. M = (    | 1 0 0 |    )
//  * 				 ( M1 | 0 1 0 | M2 )
//  * 				 (    | 0 0 1 |    )
//  * 				 (        M3       )
//  * 
//  * \param[in] m1 The row from which to start on.
//  * \param[in] m2 The last row that should be considered.
//  * \param[in] n1 The column from which to start on.
//  * \param[in] n2 The last column that should be considered.
//  * \return The error code.
//  * 
//  */
// int Matrix::set2Id( int m1, int m2, int n1, int n2 )
// {
// 	for( int i = m1; i < m2; i++ )
// 		for( int j = n1; j < n2; j++ )
// 		{
// 			if( (i - m1) == (j - n1) )
// 			{
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2const( 1.0 );		// p(x) = 1
// 				//else _data[ i ][ j ] = 1.0;
// 			}
// 			else
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0
// 				//else _data[ i ][ j ] = 0.0;
// 		}
// 	return 0;
// }

// /**
//  * Sets a matrix to zero entries.
//  * 
//  * \return The error code.
//  * 
//  */
// int Matrix::set2zero()
// {
// 	set2val( 0.0 );
// 	return 0;
// }

// /**
//  * Sets a matrix to zero entries, for especified rows and columns.
//  * 
//  * \param[in] m The number of rows to be set to zero.
//  * \param[in] n The number of columns to be set to zero.
//  * \return The error code.
//  * 
//  */
// int Matrix::set2zero( int m, int n )
// {
// 	for( int i = 0; i < m; i++ )
// 		for( int j = 0; j < n; j++ )
// 			if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 				_data[ i ][ j ].set2zero();					// p(x) = 0
// 			//else _data[ i ][ j ] = 0.0;
// 	return 0;
// }

// /**
//  * Sets a submatrix to zero:
//  * 
//  * 		e.g. M = (    | 0 0 0 |    )
//  * 				 ( M1 | 0 0 0 | M2 )
//  * 				 (    | 0 0 0 |    )
//  * 				 (        M3       )
//  * 
//  * \param[in] m1 The row from which to start on.
//  * \param[in] m2 The last row that should be considered.
//  * \param[in] n1 The column from which to start on.
//  * \param[in] n2 The last column that should be considered.
//  * \return The error code.
//  * 
//  */
// int Matrix::set2zero( int m1, int m2, int n1, int n2 )
// {
// 	for( int i = m1; i < m2; i++ )
// 		for( int j = n1; j < n2; j++ )
// 			if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 				_data[ i ][ j ].set2zero();					// p(x) = 0
// 			//else _data[ i ][ j ] = 0.0;
// 	return 0;
// }

// /**
//  * Sets a matrix to the value given as parameter.
//  * 
//  * \param[in] v The double value to set the elements to.
//  * \return The error code.
//  * 
//  */
// int Matrix::set2val( double v )
// {
// 	if( v == 0.0 )
// 	{
// 		for( int i = 0; i < _rows; i++ )
// 			for( int j = 0; j < _cols; j++ )
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2zero();				// p(x) = 0
// 				//else _data[ i ][ j ] = 0.0;
// 	}
// 	else
// 		for( int i = 0; i < _rows; i++ )
// 			for( int j = 0; j < _cols; j++ )
// 				if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 					_data[ i ][ j ].set2const( v );			// p(x) = v
// 				//else _data[ i ][ j ] = v;
// 	return 0;
// }

// /**
//  * Sets a matrix element to the value given as parameter.
//  * 
//  * \param[in] i The row position.
//  * \param[in] j The column position.
//  * \param[in] v The type \type T value to set the elements to.
//  * \return The error code.
//  * 
//  */
// int Matrix::set2val( int i, int j, double v )
// {
// 	if( v == 0.0 )
// 	{
// 		if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 			_data[ i ][ j ].set2zero();						// p(x) = 0
// 		//else _data[ i ][ j ] = 0.0;
// 	}
// 	else
// 		if( strcmp( _data[ i ][ j ].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 			_data[ i ][ j ].set2const( v );					// p(x) = v
// 		//else _data[ i ][ j ] = v;
// 	return 0;
// }

// /**
//  * Calculates the inverse (triangular inverse) of a regular, upper triangular matrix.
//  * 
//  * \param[out] Am The pointer to \a Am, an object of type \type Matrix with the inverse of \a A. 
//  * \return The error code.
//  * 
//  */
// int Matrix::trinvm( Matrix &Am )
// {
// 	int i, j, k;
// 	T sum( dimT() ), p( dimT() );
	
// 	Am.set2zero();											// some initialization
// 	p.set2const( 1.0 );
// 	for( i = 0; i < _rows; i++ )
// 	{
// 		Am( i, i ) = p / _data[ i ][ i ];					// elements of the diagonal
// 	}
// 	for( j = _rows - 2; j > -1; j-- )
// 		for( i = j + 1; i < _rows; i++ )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0
// 			//else sum = 0.0;
// 			for( k = j + 1; k < _rows; k++ )
// 				sum += _data[ j ][ k ] * Am( k, i );
// 			Am( j, i ) = - Am( j, j ) * sum;
// 		}
// 	return 0;
// }

// /**
//  * Calculates the inverse (triangular inverse) of a regular, upper triangular matrix,
//  * given the number of rows/columns that are of interest.
//  * 
//  * \param[in] r The rows/columns that are of interest.
//  * \param[out] Am The pointer to \a Am, an object of type \type Matrix with the inverse of \a A. 
//  * \return The error code.
//  * 
//  */
// int Matrix::trinvm( int r, Matrix &Am )
// {
// 	int i, j, k;
// 	T sum( dimT() ), p( dimT() );
	
// 	p.set2const( 1.0 );
// 	for( i = 0; i < r; i++ )
// 		Am( i, i ) = p / _data[ i ][ i ];					// elements of the diagonal
// 	for( j = r - 2; j > -1; j-- )
// 		for( i = j + 1; i < r; i++ )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0
// 			//else sum = 0.0;
// 			for( k = j + 1; k < r; k++ )
// 				sum += _data[ j ][ k ] * Am( k, i );
// 			Am( j, i ) = - Am( j, j ) * sum;
// 		}
// 	return 0;
// }

// /**
//  * Compares two matrices and returns \a true if the matrices are equal.
//  * Otherwise it returns \a false.
//  * 
//  * \param[in] B The pointer to \a B, an object of type \type Matrix to compare to. 
//  * \return \a true if the matrices are equal; \a false otherwise.
//  * 
//  */
// bool Matrix::mcompare( Matrix &B )
// {
// 	bool equal = true;

// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 			if( _data[ i ][ j ] != B( i, j ) )
// 			{
// 				equal = false;
// 				break;
// 			}
// 		if( !equal )
// 			break;
// 	}
// 	return equal;
// }

// /**
//  * Compares two matrices and returns \a true if the matrices are equal.
//  * Otherwise it returns \a false.
//  * 
//  * \param[in] B The pointer to \a B, an object of type \type Matrix to compare to. 
//  * \param[in] eps The threshold value for comparing. 
//  * \return \a true if the matrices are equal; \a false otherwise.
//  * 
//  */
// bool Matrix::mcompare( Matrix &B, double eps )
// {
// 	bool equal = true;
// 	T diff( dimT() );

// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			diff = _data[ i ][ j ] - B( i, j );
// 			if( !diff.isZero( eps ) )
// 			{
// 				equal = false;
// 				break;
// 			}
// 		}
// 		if( !equal )
// 			break;
// 	}
// 	return equal;
// }

// /**
//  * Compares two matrices and returns \a true if the matrices are equal up to the
//  * number of rows of the first matrix. Otherwise it returns \a false.
//  * 
//  * \param[in] B The pointer to \a B, an object of type \type Matrix to compare to. 
//  * \param[in] r The number of rows to be compared to, for all columns. 
//  * \param[in] eps The threshold value for comparing. 
//  * \return \a true if the matrices are equal; \a false otherwise.
//  * 
//  */
// bool Matrix::mcompare( Matrix &B, int r, double eps )
// {
// 	bool equal = true;
// 	T diff( dimT() );

// 	for( int i = 0; i < r; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			diff = _data[ i ][ j ] - B( i, j );
// 			if( !diff.isZero( eps ) )
// 			{
// 				equal = false;
// 				break;
// 			}
// 		}
// 		if( !equal )
// 			break;
// 	}
// 	return equal;
// }

// //
// // P R I N T   O U T   F U N C T I O N S
// //

/**
 * Prints out a matrix in a given format, with the previously set color.
 * 
 * \param[in] str The character string to be printed out.
 * 
*/
void Matrix::printm( char *str )
{
 	printf( str );
 	for( int i = 0; i < _rows; i++ )
 	{
 		for( int j = 0; j < _cols; j++ )
 		{
 			printf("%3.2f", _data[ i ][ j ]);
 			printf( "%c", '\t' );
 		}
 		printf( "\n" );
 	}
}
// /**
//  * Prints out to a file a matrix in a given format, with the previously set color.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] str The character string to be printed out.
//  * 
//  */
// void Matrix::fprintm( FILE * fn, char *str )
// {
// 	fprintf( fn, str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			_data[ i ][ j ].print( fn );
// 			fprintf( fn, "%c", '\t' );
// 		}
// 		fprintf( fn, "\n" );
// 	}
// }

// /**
//  * Prints out a matrix in a given format.
//  * 
//  * \param[in] str The character string to be printed out.
//  * \param[in] color The desired color for the output.
//  * 
//  */
// void Matrix::printm( char *str, const char *const color )
// {
// 	printf( color );
// 	printf( str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			_data[ i ][ j ].print();
// 			printf( "%c", '\t' );
// 		}
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }

// /**
//  * Prints out a matrix in a given format.
//  * 
//  * \param[in] str The character string to be printed out.
//  * 
//  */
// void Matrix::printdm( char *str )
// {
// 	printf( str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		printf( "s", "  " );
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			_data[ i ][ j ].print();
// 			printf( "%c", '\t' );
// 		}
// 		//printf( ")\n" );
// 		printf( "\n" );
// 	}
// }

// /**
//  * Prints out to a file a matrix in a given format.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] str The character string to be printed out.
//  * 
//  */
// void Matrix::fprintdm( FILE * fn, char *str )
// {
// 	fprintf( fn, str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		//fprintf( fn, "%c(%c", '\t','\t' );
// 		fprintf( fn, "%s", "  " );
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			_data[ i ][ j ].print( fn );
// 			fprintf( fn, "%c", '\t' );
// 		}
// 		//fprintf( fn, ")\n" );
// 		fprintf( fn, "\n" );
// 	}
// }

// /**
//  * Prints out a matrix in a given format.
//  * 
//  * \param[in] str The character string to be printed out.
//  * \param[in] color The desired color for the output.
//  * 
//  */
// void Matrix::printdm( char *str, const char *const color )
// {
// 	printf( color );
// 	printf( str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		//printf( "%c(%c", '\t','\t' );
// 		printf( "%s", "  " );
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			_data[ i ][ j ].print();
// 			printf( "%c", '\t' );
// 		}
// 		//printf( ")\n" );
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }

// /**
//  * Prints out a matrix in a given format.
//  * 
//  * \param[in] str The character string to be printed out.
//  * \param[in] color The desired color for the output.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::printdm( char *str, const char *const color, double eps )
// {
// 	printf( color );
// 	printf( str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		//printf( "%c(%c", '\t','\t' );
// 		printf( "%s", "  " );
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			if( _data[ i ][ j ].isZero( eps ) )				//small enough
// 				printf( "%.16g%c", 0.0, '\t' );
// 			else
// 			{
// 				_data[ i ][ j ].print();
// 				printf( "%c", '\t' );
// 			}
// 		}
// 		//printf( ")\n" );
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }
// /**
//  * Prints out a matrix in a given format.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] str The character string to be printed out.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::fprintdm( FILE * fn, char *str, double eps )
// {
// 	fprintf( fn, str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		//fprintf( fn, "%c(%c", '\t', '\t' );
// 		fprintf( fn, "%s", "  " );
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			if( _data[ i ][ j ].isZero( eps ) )				//small enough
// 				fprintf( fn, "%.16g%c", 0.0, '\t' );
// 			else
// 			{
// 				_data[ i ][ j ].print( fn );
// 				fprintf( fn, "%c", '\t' );
// 			}
// 		}
// 		//fprintf( fn, ")\n" );
// 		fprintf( fn, "\n" );
// 	}
// }

// /**
//  * Prints out a matrix in a given format, with column pivoting.
//  * 
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
//  * \param[in] str The character string to be printed out.
//  * \param[in] color The desired color for the output.
//  * 
//  */
// void Matrix::printm( int *piv, char *str, const char *const color )
// {
// 	printf( color );
// 	printf( str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			_data[ i ][ piv[ j ] ].print();
// 			printf( "%c", '\t' );
// 		}
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }

// /**
//  * Prints out a matrix in a given format, with column pivoting.
//  * 
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
//  * \param[in] str The character string to be printed out.
//  * \param[in] color The desired color for the output.
//  * 
//  */
// void Matrix::printdm( int *piv, char *str, 
// 	const char *const color, double eps )
// {
// 	printf( color );
// 	printf( str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		//printf( "%c(%c", '\t','\t' );
// 		printf( "%s", "  " );
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			if( _data[ i ][ piv[ j ] ].isZero( eps ) )		//small enough
// 				printf( "%.16g%c", 0.0, '\t' );
// 			else
// 			{
// 				_data[ i ][ piv[ j ] ].print();
// 				printf( "%c", '\t' );
// 			}
// 		}
// 		//printf( ")\n" );
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }

// /**
//  * Prints out a matrix in a given format, with column pivoting.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
//  * \param[in] str The character string to be printed out.
//  * 
//  */
// void Matrix::fprintm( FILE * fn, int *piv, char *str )
// {
// 	fprintf( fn, str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			_data[ i ][ piv[ j ] ].print( fn );
// 			fprintf( fn, "%c", '\t' );
// 		}
// 		fprintf( fn, "\n" );
// 	}
// }

// /**
//  * Prints out a matrix in a given format, with column pivoting.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
//  * \param[in] str The character string to be printed out.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::fprintdm( FILE * fn, int *piv, char *str, double eps )
// {
// 	fprintf( fn, str );
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		for( int j = 0; j < _cols; j++ )
// 		{
// 			if( _data[ i ][ piv[ j ] ].isZero( eps ) )		//small enough
// 				fprintf( fn, "%.16g%c", 0.0, '\t' );
// 			else
// 			{
// 				_data[ i ][ piv[ j ] ].print( fn );
// 				fprintf( fn, "%c", '\t' );
// 			}
// 		}
// 		fprintf( fn, "\n" );
// 	}
// }

// /**
//  * Prints out a matrix of Taylor polynomials in a given format.
//  * 
//  * \param[in] str The character string to be printed out.
//  * \param[in] color The desired color for the output.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::printtpm( char *str, const char *const color, double eps )
// {
// 	int _coeff = _data[ 0 ][ 0 ].nrcoeff();
// 	printf( color );
// 	printf( str );
// 	for( int k = 0; k < _coeff; k++ )
// 	{
// 		for( int i = 0; i < _rows; i++  )
// 		{
// 			//printf( "%c(%c", '\t','\t' );
// 			printf( "%s", "  " );
// 			for( int j = 0; j < _cols; j++ )
// 			{
// 				if( abs( _data[ i ][ j ][ k ] ) < eps )		//small enough
// 					printf( "%.16lg%c", 0.0, '\t' );
// 				else
// 					printf( "%.16lg%c", _data[ i ][ j ][ k ], '\t' );
// 			}
// 			//printf( ")\n" );
// 			printf( "\n" );
// 		}
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }

// /**
//  * Prints out a submatrix of Taylor polynomials in a given format.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] str The character string to be printed out.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::fprinttpm( FILE * fn, char *str, double eps )
// {
// 	int _coeff = _data[ 0 ][ 0 ].nrcoeff();
// 	fprintf( fn, str );
// 	for( int k = 0; k < _coeff; k++ )
// 	{
// 		for( int i = 0; i < _rows; i++  )
// 		{
// 			//fprintf( fn, "%c(%c", '\t','\t' );
// 			fprintf( fn, "%s", "  " );
// 			for( int j = 0; j < _cols; j++ )
// 			{
// 				if( abs( _data[ i ][ j ][ k ] ) < eps )		//small enough
// 					fprintf( fn, "%.16lg%c", 0.0, '\t' );
// 				else
// 					fprintf( fn, "%.16lg%c", _data[ i ][ j ][ k ], '\t' );
// 			}
// 			//fprintf( fn, ")\n" );
// 			fprintf( fn, "\n" );
// 		}
// 		fprintf( fn, "\n" );
// 	}
// }

// /**
//  * Prints out a submatrix of Taylor polynomials in a given format.
//  * 
//  * \param[in] str The character string to be printed out.
//  * \param[in] color The desired color for the output.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::printtpm( int m1, int m2, int n1, int n2,
// 		char *str, const char *const color, double eps )
// {
// 	int _coeff = _data[ 0 ][ 0 ].nrcoeff();
// 	printf( color );
// 	printf( str );
// 	for( int k = 0; k < _coeff; k++ )
// 	{
// 		for( int i = m1; i < m2; i++  )
// 		{
// 			//printf( "%c(%c", '\t','\t' );
// 			printf( "%s", "  " );
// 			for( int j = n1; j < n2; j++ )
// 			{
// 				if( abs( _data[ i ][ j ][ k ] ) < eps )		//small enough
// 					printf( "%.16lg%c", 0.0, '\t' );
// 				else
// 					printf( "%.16lg%c", _data[ i ][ j ][ k ], '\t' );
// 			}
// 			//printf( ")\n" );
// 			printf( "\n" );
// 		}
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }

// /**
//  * Prints out a submatrix of Taylor polynomials in a given format.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] str The character string to be printed out.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::fprinttpm( FILE * fn, int m1, int m2, int n1, int n2,
// 		char *str, double eps )
// {
// 	int _coeff = _data[ 0 ][ 0 ].nrcoeff();
// 	fprintf( fn, str );
// 	for( int k = 0; k < _coeff; k++ )
// 	{
// 		for( int i = m1; i < m2; i++  )
// 		{
// 			//fprintf( fn, "%c(%c", '\t','\t' );
// 			fprintf( fn, "%s", "  " );
// 			for( int j = n1; j < n2; j++ )
// 			{
// 				if( abs( _data[ i ][ j ][ k ] ) < eps )		//small enough
// 					fprintf( fn, "%.16lg%c", 0.0, '\t' );
// 				else
// 					fprintf( fn, "%.16lg%c", _data[ i ][ j ][ k ], '\t' );
// 			}
// 			//fprintf( fn, ")\n" );
// 			fprintf( fn, "\n" );
// 		}
// 		fprintf( fn, "\n" );
// 	}
// }

// /**
//  * Prints out a matrix of Taylor polynomials in a given format.
//  * 
//  * \param[in] str The character string to be printed out.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
//  * \param[in] color The desired color for the output.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::printtpm( int *piv, char *str, 
// 	const char *const color, double eps )
// {
// 	int _coeff = _data[ 0 ][ 0 ].nrcoeff();
// 	printf( color );
// 	printf( str );
// 	for( int k = 0; k < _coeff; k++ )
// 	{
// 		for( int i = 0; i < _rows; i++  )
// 		{
// 			//printf( "%c(%c", '\t','\t' );
// 			printf( "%s", "  " );
// 			for( int j = 0; j < _cols; j++ )
// 			{
// 				if( abs( _data[ i ][ piv[ j ] ][ k ] ) < eps )//small enough
// 					printf( "%.16lg%c", 0.0, '\t' );
// 				else
// 					printf( "%.16lg%c", _data[ i ][ piv[ j ] ][ k ], '\t' );
// 			}
// 			//printf( ")\n" );
// 			printf( "\n" );
// 		}
// 		printf( "\n" );
// 	}
// 	printf( normal );
// }

// /**
//  * Prints out a matrix of Taylor polynomials in a given format.
//  * 
//  * \param[in] fn The file name where to write the data to.
//  * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a M.
//  * \param[in] str The character string to be printed out.
//  * \param[in] eps The threshold value to print less.
//  * 
//  */
// void Matrix::fprinttpm( FILE * fn, int *piv, char *str, double eps )
// {
// 	int _coeff = _data[ 0 ][ 0 ].nrcoeff();
// 	fprintf( fn, str );
// 	for( int k = 0; k < _coeff; k++ )
// 	{
// 		for( int i = 0; i < _rows; i++  )
// 		{
// 			//fprintf( fn, "%c(%c", '\t','\t' );
// 			fprintf( fn, "%s", "  " );
// 			for( int j = 0; j < _cols; j++ )
// 			{
// 				if( abs( _data[ i ][ piv[ j ] ][ k ] ) < eps )//small enough
// 					fprintf( fn, "%.16lg%c", 0.0, '\t' );
// 				else
// 					fprintf( fn, "%.16lg%c", _data[ i ][ piv[ j ] ][ k ], '\t' );
// 			}
// 			//fprintf( fn, ")\n" );
// 			fprintf( fn, "\n" );
// 		}
// 		fprintf( fn, "\n" );
// 	}
// }






/***************
  P R I V A T E
  **************/

void Matrix::allocateMemory(bool initialize)
{
	try
	{
		// bounds checking
		if( _rows <= 0 || _cols <= 0 )
			throw exception( "Matrix invalid matrix size", 35 ); // TODO what means 35?, add invalid size to message

		_data = new double*[_rows];
		if( _data == 0 || _data == NULL )
			throw exception( "Memory allocation failure.", 3 ); // TODO what meansi 3?
		
		for( int i = 0; i < _rows; i++ )
		{
			_data[i] = new double[_cols];
			if( _data[i] == 0 || _data[i] == NULL )
				throw exception( "Memory allocation failure.", 3 ); // TODO
			
			if(initialize)
			{
				for( int j = 0; j < _cols; j++)
				{
					_data[i][j] = 0.0;
					// _data[i][j].set2zero();
					// _data[i][j] = TPoly(_dimT);
				}
			}
		}
	}
	catch( bad_alloc e )
	{
		// TODO bad!
		printf( "" );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( "" );
        throw 4;
	}
	catch(...)
	{
		exception e( "Error when allocating a matrix.", 34 );
        // e.report();
		// throw e.;
	}
}

void Matrix::deallocateMemory()
{
	for( int i = 0; i < _rows; i++ )
		delete [] _data[ i ];
	delete [] _data;
}

void Matrix::copyFrom(const Matrix &m)
{
	// deallocateMemory();

	_rows = m._rows;
	_cols = m._cols;
	_dimT = m._dimT;

	allocateMemory(false);

	for( int i = 0; i < _rows; i++ )
		for( int j = 0; j < _cols; j++)
			_data[i][j] = m._data[i][j];
}
