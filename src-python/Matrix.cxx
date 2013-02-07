#include "example.h"

using namespace std;

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
	_dimT(0),
	_data(0)
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
	_dimT(0),
	_data(0)
{
	allocateMemory(initialize);
}

/**
 * Constructor for the class with both the number of rows and columns as parameters,
 * as well as the dimension of the elements' type, e.g., the Taylor polynomial's grade. 
 * Creates the object.
 * 
 * \param[in] r The number of rows.
 * \param[in] c The number of columns.
 * \param[in] dimT The dimension of the type \type T.
 * 
 */
Matrix::Matrix(int rows, int cols, int dimT, bool initialize):
	_rows(rows),
	_cols(cols),
	_dimT(dimT),
	_data(0)
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
Matrix::Matrix(const Matrix &m):
	_rows(0),
	_cols(0),
	_dimT(0),
	_data(0)
{
	copyFrom(m);
}

/**
 * Easy constructor for testing.
 *
 */
Matrix::Matrix(int rows, int cols, double *values):
	_rows(rows),
	_cols(cols),
	_dimT(0),
	_data(0)
{
	allocateMemory(false);

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] = values[i * _cols + j];
		}
	}
}

//
// D E S T R U C T O R
//

/**
 * Destructor. Cleans up the object.
 * 
 */
Matrix::~Matrix()
{
	deallocateMemory();
}

//
// G E T T E R
//

/**
 * Gets a single element from the matrix.
 * 
 * \param[in] row The row index.
 * \param[in] col The column index.
 * \return The value of the desired element.
 *
 */
double Matrix::get(int row, int col) const
{ 
	// bounds checking
	if( row >= _rows  ||  row < 0  ||  col >= _cols  ||  col < 0 )
	{
		throw MathException("(%d, %d) is not a valid index of a %d x %d matrix.", row, col, _rows, _cols);
	}
	return _data[row][col];
}


//
// O V E R L O A D E D   ( )   A N D   =   O P E R A T O R S
//

/**
 * Implements the () operator.
 * 
 * (T & Matrix::... with '&' allows m(2, 1)) = 7, i.e., '( )' also in the left side!!)
 * 
 * \param[in] i The row index.
 * \param[in] j The column index.
 * \return The value of the desired element.
 *
 */
double &Matrix::operator()(int row, int col)
{ 
	// bounds checking
	if( row >= _rows  ||  row < 0  ||  col >= _cols  ||  col < 0 )
	{
		throw MathException("(%d, %d) is not a valid index of a %d x %d matrix.", row, col, _rows, _cols);
	}
	return _data[row][col];
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
	if(this != &m)
	{
	 	copyFrom(m);
	}

	return *this;
}

//
// O V E R L O A D E D   A R I T H M E T I C   O P E R A T O R S
//

/**
 * Implements the == operator. It compares two matrices.
 * 
 * \param[in] m The matrix to compare with.
 * \return \a true if the matrices are equal. Otherwise it returns \a false.
 * 
 */
bool Matrix::operator==( const Matrix &m ) const
{
	if(_rows != m._rows || _cols != m._cols) // TODO _dimT ?
	{
		return false;
	}
	
	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			if( _data[i][j] != m._data[i][j] )
			{
				return false;
			}
		}
	}

	return true;
}

/**
 * Implements the != operator. It compares two matrices.
 * 
 * \param[in] m The matrix to compare with.
 * \return \a true if the matrices are not equal. Otherwise it returns \a false.
 * 
 */
bool Matrix::operator!=( const Matrix &m ) const
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
Matrix Matrix::operator+(const Matrix &m) const
{	
	// dimensions checking
	if(_rows != m._rows || _cols != m._cols)
	{
		throw MathException("Cannot add up a %d x %d matrix (left side) and a %d x %d (right side) matrix.", _rows, _cols, m._rows, m._cols);
	}

	// an auxiliary object
	Matrix aux(_rows, _cols, _dimT, false);

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			aux._data[i][j] = _data[i][j] + m._data[i][j];
		}
	}

	return aux;
}

/**
 * Implements the += operator. It adds up two matrices.
 * 
 * \param[in] m The \a Matrix object to add up to the current matrix.
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator+=(const Matrix &m)
{	
	// dimensions checking
	if(_rows != m._rows || _cols != m._cols)
	{
		throw MathException("Cannot add up a %d x %d matrix to a %d x %d matrix.", m._rows, m._cols, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{	
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] += m._data[i][j];
		}
	}

 	return *this;
}

/**
 * Implements the - operator. It substracts two matrices.
 * 
 * \param[in] m The \a Matrix object to substract from.
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator-(const Matrix &m) const
{	
	// dimensions checking
	if(_rows != m._rows || _cols != m._cols)
	{
		throw MathException("Cannot subtract a %d x %d from a %d x %d matrix.", m._rows, m._cols, _rows, _cols);
	}

	// an auxiliary object
	Matrix aux(_rows, _cols, _dimT, false);

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			aux._data[i][j] = _data[i][j] - m._data[i][j];
		}
	}

	return aux;
}

/**
 * Implements the -= operator. It substracts two matrices.
 * 
 * \param[in] m The \a Matrix object to substract from.
 * \return The resulting matrix.
 * 
 */
Matrix Matrix::operator-=(const Matrix &m)
{		
	// dimensions checking
	if(_rows != m._rows || _cols != m._cols)
	{
		// TODO more information?
		throw MathException("Cannot substract a %d x %d matrix from a %d x %d matrix.", m._rows, m._cols, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] -= m._data[i][j];
		}
	}

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
	// an auxiliary object
	Matrix aux( _rows, _cols, _dimT, false);
	
	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			aux._data[i][j] = - _data[i][j];
		}
	}

	return aux;
}

/**
 * Implements the * operator. It multiplies a matrix by a scalar.
 * 
 * \param[in] alpha The scalar to multiply by.
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator*(double alpha) const
{
	// an auxiliary object
	Matrix aux(_rows, _cols, _dimT, false);

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			aux._data[i][j] = _data[i][j] * alpha;
		}
	}

	return aux;
}

/**
 * Implements the *= operator. It multiplies the matrix by a scalar.
 * 
 * \param[in] alpha The scalar to multiply by.
 * \return A pointer to the resulting \a Matrix object.
 * 
 */
Matrix Matrix::operator*=(double alpha)
{
	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] *= alpha;
		}
	}

	return *this;
}

/**
 * Implements the * operator. It multiplies the matrix with another matrix.
 *
 * \param[in] m The matrix to multiply with.
 * \return A pointer ot the resulting \a Matrix object.
 *
 */
Matrix Matrix::operator*(const Matrix &m) const
{
	if (_cols != m._rows) 
	{
		throw MathException("Cannot multiply a %d x %d with a %d x %d matrix. The number of rows in m must match the number of colmuns in this matrix.", _rows, _cols, m._rows, m._cols);
	}

	Matrix aux( _rows, m._cols, _dimT, true);

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < m._cols; j++ )
		{
			// because we initliazed aux all fields are already set to 0
			for( int k = 0; k < _cols; k++ )
			{
				aux._data[i][j] += _data[i][k] * m._data[k][j];
			}
		}
	}

	return aux;	
}

//
// E S P E C I A L   M A T R I X   M U L T I P L I C A T I O N S
//

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * B + beta * C
 * 
 * with A : m-by-p matrix
 * 		B : p-by-n matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 * 
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 * 
 */
void Matrix::mmCaABbC(double alpha, double beta, const Matrix &A, const Matrix &B)
{
	if (A._cols != B._rows)
	{
		throw MathException("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols);
	}
	if (A._rows != _rows || B._cols != _cols)
	{
		throw MathException("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);
			
			for( int k = 0; k < B._rows; k++ )
			{
				h += A._data[i][k] * B._data[k][j];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * B + beta * C
 * 
 * with A : m-by-p matrix
 * 		B : p-by-n matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 * 
 * where the inferior-right block of B is an identity matrix like in:
 * 
 * 		( * * * 0 0 )
 * 		( * * * 0 0 )
 * 		( 0 0 0 1 0 )
 * 		( 0 0 0 0 1 )
 * 
 * so that a particular block multiplication is needed.
 * 
 * \param[in] r The number of rows in \a B that are of interest (2 in the example above).
 * \param[in] c The number of columns in \a B that are of interest (3 in the example above).
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 *
 */
void Matrix::bmmCaABbC(int r, int c, double alpha, double beta, const Matrix &A, const Matrix &B)
{
	if (r > B._rows)
	{
		throw MathException("The parameter r (= %d) must be smaller or equal than the number of rows in B (%d x %d).", r, B._rows, B._cols);
	}
	if (c > B._cols)
	{
		throw MathException("The parameter c (= %d) must be smaller or equal than the number of columns in B (%d x %d).", c, B._rows, B._cols);
	}
	if (A._cols != B._rows)
	{
		throw MathException("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols);
	}
	if (A._rows != _rows || B._cols != _cols)
	{
		throw MathException("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols);
	}

	int cr = c - r;
	// If you wonder why cr = c-r and latter A(i, j-cr)?
	// I don't know, bit take a sheet of paper the the matrices and you will see ;)

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < c; j++ ) 
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < r; k++ )
			{
				h += A._data[i][k] * B._data[k][j];
			}	

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
		
		for( int j = c; j < _cols; j++ )
		{
			_data[i][j] = A._data[i][j-cr] * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * B + beta * C
 * 
 * with A : m-by-p matrix
 * 		B : p-by-n matrix
 * 		C : r-by-n matrix (only the last r rows from A are interesting)
 * 		alpha, beta : real numbers
 * 
 * where A ("special" A) is of the form:
 * 
 * 		(         )
 * 		(    X    )
 * 		( ------- )
 * 		( * ... * )
 * 		( * ... * )
 * 
 * so that a particular matrix multiplication is needed.
 * 
 * \param[in] r The last rows in \a A that are of interest (2 non-zero-rows in the example above).
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 *
 */
void Matrix::mmCasABbC(int r, double alpha, double beta, const Matrix &A, const Matrix &B)
{
	if (r > A._rows)
	{
		throw MathException("The number of rows in A that are of interest (r = %d) must be smaller or equal than the total amount of rows in A (%d x %d).", r, A._rows, A._cols);
	}
	if (A._cols != B._rows)
	{
		throw MathException("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols);
	}
	if (A._rows != _rows || B._cols != _cols)
	{
		throw MathException("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols);
	}

	int n = A._rows - r;
	for( int i = 0; i < n; i++ )
	{
		for( int j = 0; j < B._cols; j++ )
		{
			_data[i][j] *= beta;
		}
	}

	for( int i = n; i < A._rows; i++ )
	{
		for( int j = 0; j < B._cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < B._rows; k++ )
			{
				h += A._data[i][k] * B._data[k][j];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * B + beta * C
 * 
 * with A : m-by-p matrix
 * 		B : p-by-n matrix
 * 		C : m-by-r matrix (only the last r columns from B are interesting).
 * 		alpha, beta : real numbers
 * 
 * where B ("special" B) is of the form:
 * 
 * 		(   | * * * )
 * 		(   | . . . )
 * 		( X | . . . )
 * 		(   | . . . )
 * 		(   | * * * )
 * 
 * so that a particular matrix multiplication is needed.
 * 
 * \param[in] r The last columns in \a B that are of interest (3 non-zero-columns 
 * 				in the example above).
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 *
 */
void Matrix::mmCaAsBbC(int r, double alpha, double beta, const Matrix &A, const Matrix &B)
{
	if (r > B._cols)
	{
		throw MathException("The number of columns in B that are of interest (r = %d) must be smaller or equal than the total amount of columns in B (%d x %d).", r, A._rows, A._cols);
	}
	if (A._cols != B._rows)	
	{
		throw MathException("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols);
	}
	if (A._rows != _rows || B._cols != _cols)
	{
		throw MathException("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols);
	}

	int n = B._cols - r;
	for( int i = 0; i < _rows; i++ )
	{
		// zero-columns of B
		for( int j = 0; j < n; j++ )
		{
			_data[i][j] *= beta;
		}

		// non-zero-columns of B
		for( int j = n; j < _cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < B._rows; k++ )
			{
				h += A._data[i][k] * B._data[k][j];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * UTB + beta * C
 * 
 * where UTB means that only the upper triangular part is of interest. Furthermore, a column
 * pivoting on B is considered.
 * 
 * with A : m-by-p matrix
 * 		B : p-by-n matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 * 
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix. Only its upper triangular 
 * 		part is interesting.
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a B.
 *
 */
void Matrix::mmCaAUTBPbC(double alpha, double beta, const Matrix &A, const Matrix &B, int *piv)
{
	if (A._cols != B._rows)	
	{
		throw MathException("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols);
	}
	if (A._rows != _rows || B._cols != _cols)
	{
		throw MathException("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			double h = 0.0;
			// TPolyn h(_dimT);

			for( int k = 0; k <= j; k++ )
			{
				h += A._data[i][k] * B._data[k][piv[j]];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * A^T + beta * C
 * 
 * with A, C : m-by-m matrix
 * 		alpha, beta : real numbers
 * 
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is also considered.
 *
 */

void Matrix::mmCaAATbC(double alpha, double beta, const Matrix &A)
{
	// if A is a m-by-n matrix A * A^T is always a m-by-m matrix and must have the same dimension as this matrix
	if (A._rows != _rows || A._rows != _cols)
	{
		throw MathException("The size of the matrix A * A^T (%d x %d) must match the current matrix (%d x %d).", A._rows, A._rows, _rows, _cols);
	}
	
	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _rows; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < _rows; k++ )
			{
				h += A._data[i][k] * A._data[j][k];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A^T * A + beta * C
 * 
 * with A, C : m-by-m matrix
 * 		alpha, beta : real numbers
 *
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is also considered.
 *
 */

void Matrix::mmCaATAbC(double alpha, double beta, const Matrix &A)
{
	// if A is a m-by-n matrix A^T * A is always a n-by-n matrix and must have the same dimension as this matrix
	if (A._cols != _rows || _cols != _rows)
	{
		throw MathException("The size of the matrix A^T * A (%d x %d) must match the current matrix (%d x %d).", A._cols, A._cols, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _rows; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < A._rows; k++ )
			{
				h += A._data[k][i] * A._data[k][j];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A^T * B + beta * C
 * 
 * with A : p-by-m matrix
 * 		B : p-by-n matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 *
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is considered.
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 *
 */

void Matrix::mmCaATBbC(double alpha, double beta, const Matrix &A, const Matrix &B)
{
	// A^T und B can only be multiplied if A and B have the same number of rows
	if (A._rows != B._rows)
	{
		throw MathException("Cannot multiply A^T (%d x %d) und B (%d x %d).", A._cols, A._rows, B._rows, B._cols);
	}

	// (A^T * B) must have the same size as this matrix
	if (A._cols != _rows || B._cols != _cols)
	{
		throw MathException("The size of the matrix A^T * B (%d x %d) must match the current matrix (%d x %d).", A._cols, B._cols, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < B._rows; k++ )
			{
				h += A._data[k][i] * B._data[k][j];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A^T * B + beta * C
 * 
 * with A : p-by-m matrix
 * 		B : p-by-n matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 *
 * and a column pivoting on A^Ts rows.
 * 
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix. Its transpose is considered.
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a B.
 *
 */

void Matrix::mmCaATBPbC(double alpha, double beta, const Matrix &A, const Matrix &B, int *piv)
{
	// A^T und B can only be multiplied if A and B have the same number of rows
	if (A._rows != B._rows)
	{
		throw MathException("Cannot multiply A^T (%d x %d) und B (%d x %d).", A._cols, A._rows, B._rows, B._cols);
	}

	// (A^T * B) must have the same size as this matrix
	if (A._cols != _rows || B._cols != _cols)
	{
		throw MathException("The size of the matrix A^T * B (%d x %d) must match the current matrix (%d x %d).", A._cols, B._cols, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < B._rows; k++ )
			{
				h += A._data[k][i] * B._data[k][piv[j]];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * B^T + beta * C
 * 
 * with A : m-by-p matrix
 * 		B : n-by-p matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 *
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix. Its transpose is considered.
 * 
 */

void Matrix::mmCaABTbC(double alpha, double beta, const Matrix &A, const Matrix &B)
{
	// A und B^T can only be multiplied if A and B have the same number of columns
	if (A._cols != B._cols)
	{
		throw MathException("Cannot multiply A (%d x %d) und B^T (%d x %d).", A._rows, A._cols, B._cols, B._rows);
	}

	// (A * B^T) must have the same size as this matrix
	if (A._rows != _rows || B._rows != _cols)
	{
		throw MathException("The size of the matrix A * B^T (%d x %d) must match the current matrix (%d x %d).", A._rows, B._rows, _rows, _cols);
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < A._cols; k++ )
			{
				h += A._data[i][k] * B._data[j][k];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * B^T + beta * C
 * 
 * with A : m-by-p matrix
 * 		B : n-by-p matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 *
 * After transposing B, either its first or its last rows are considered for multiplication,
 * 
 * 		( * ... * )				(    X    )
 * 		( ------- )		or		( ------- )
 * 		(         )				( * ... * )
 * 		(    X    )				( * ... * )
 * 
 * according to A dimensions. I.e., the matrix A has less columns than B^T rows has.
 * 
 * \param[in] r The number of rows from \a B that should be considered.
 * \param[in] up The binary parameter to indicate whether the first or the last \a r rows 
 * 		from \a B should be considered (=1 if first rows; =0 otherwise).
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix. Its transpose is considered.
 *
 */
void Matrix::mmCaABTbC(int r, bool up, double alpha, double beta, const Matrix &A, const Matrix &B)
{
	// A und B^T can only be multiplied if A and B have the same number of columns
	if (A._cols != B._cols)
	{
		throw MathException("Cannot multiply A (%d x %d) und B^T (%d x %d).", A._rows, A._cols, B._cols, B._rows);
	}

	// (A * B^T) must have the same size as this matrix
	if (A._rows != _rows || B._rows != _cols)
	{
		throw MathException("Error in matrix multiplication. The result of A * B must have the same size as C (this matrix).");
	}

	// r <= (colmuns of B^T), which is equal to r <= (rows of B)
	if (r > B._rows)
	{
		throw MathException("Error in matrix multiplication. r must be smaller or equal than the number of rows in B.");
	}

	int bRowsR = B._rows - r;

	for( int i = 0; i < A._rows; i++ )
	{
		for( int j = 0; j < B._cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			if (up)
			{
				// first r columns from B are used
				for( int k = 0; k < r; k++ )
				{
					h += A._data[i][k] * B._data[j][k];
				}
			}
			else
			{
				// last r columns from B are used
				for( int k = bRowsR; k < B._rows; k++ )
				{
					h += A._data[i][k] * B._data[j][k];
				}
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * B^T + beta * C
 * 
 * with A : m-by-p matrix
 * 		B : n-by-p matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 *
 * where the inferior-right block of A is an identity matrix like in:
 * 
 * 		( * * * 0 0 )
 * 		( * * * 0 0 )
 * 		( 0 0 0 1 0 )
 * 		( 0 0 0 0 1 )
 * 
 * so that a particular block multiplication is needed.
 * 
 * \param[in] r The number of rows in \a A that are of interest (2 in the example above).
 * \param[in] c The number of columns in \a A that are of interest (3 in the example above).
 * \param[in] alpha The scalar value that multiplies \a A * B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] B The pointer to \a B, an object of type \type Matrix. Its transpose is considered.
 */
void Matrix::bmmCaABTbC(int r, int c, double alpha, double beta, const Matrix &A, const Matrix &B)
{
	if (r > A._rows)
	{
		throw MathException("Error in matrix multiplication. r cannot be larger than the number of rows of A in total.");
	}
	if (c > A._cols)
	{
		throw MathException("Error in matrix multiplication. c cannot be larger than the number of columns of A in total.");
	}
	if (A._cols != B._cols)
	{
		throw MathException("Error in matrix multiplication. The matrices A und B^T cannot be multiplied because of wrong dimensions.");
	}
	if (A._rows != _rows || B._rows != _cols)
	{
		throw MathException("Error in matrix multiplication. The dimension of the matrix A * B^T must match the current matrix.");
	}

	// only the first r rows from A are interesting
	for( int i = 0; i < r; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			double h = 0.0;
			// TPoly h(_dimT);

			for( int k = 0; k < c; k++ )
			{
				h += A._data[i][k] * B._data[j][k];
			}

			_data[i][j] = h * alpha + _data[i][j] * beta;
		}
	}

	// last rows of A, with exactly one 1 per row
	int d = A._cols - A._rows;
	for ( int i = r; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] = B._data[j][i + d] * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * I * B + beta * C
 * 
 * with I : m-by-p matrix; identity matrix
 * 		B : p-by-n matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 * 
 * \param[in] alpha The scalar value that multiplies \a B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 */
void Matrix::mmCaIBbC(double alpha, double beta, const Matrix &B)
{
	if (B._rows != _rows || B._cols != _cols)
	{
		throw MathException("Error in matrix multiplication. The dimension of the matrix B must match the current matrix.");
	}	

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] = B._data[i][j] * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * I * B + beta * C
 * 
 * with I : m-by-p matrix; identity matrix permuted according to a vector of permutations, piv
 * 		B : p-by-n matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 * 
 * \param[in] alpha The scalar value that multiplies \a B.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a I, of type \type int.
 * \param[in] rows The binary parameter to indicate whether the rows or the columns of I
 * 		should be permuted (= 0, the rows; = 1, the columns).
 * \param[in] B The pointer to \a B, an object of type \type Matrix.
 */
void Matrix::mmCaIBbC(double alpha, double beta, int *piv, bool rows, const Matrix &B)
{
		// Matrix Id( _rows, B._rows, dimT() );			// Id[m][p][nTcoeff]
		// Id.set2Id();										// set to the identity
		// if( rows )											// permute its columns
		// 	Id.cpermutem( piv );
		// else												// permute its rows
		// 	Id.rpermutem( piv );
		// mmCaABbC( alpha, beta, Id, B );
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * I + beta * C
 * 
 * with A : m-by-p matrix
 * 		I : p-by-n matrix; identity matrix
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 * 
 * \param[in] alpha The scalar value that multiplies \a A.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * 
 */
void Matrix::mmCaAIbC(double alpha, double beta, const Matrix &A)
{
	if (A._rows != _rows || A._cols != _cols)
	{
		throw MathException("Error in matrix multiplication. The dimension of the matrix A must match the current matrix.");
	}	

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] = A._data[i][j] * alpha + _data[i][j] * beta;
		}
	}
}

/**
 * Matrix multiplication of the form:
 * 
 *     C = alpha * A * I + beta * C
 * 
 * with A : m-by-p matrix
 * 		I : p-by-n matrix; identity matrix permuted according to a vector of permutations, piv
 * 		C : m-by-n matrix
 * 		alpha, beta : real numbers
 * 
 * \param[in] alpha The scalar value that multiplies \a A.
 * \param[in] beta The scalar value that multiplies \a C.
 * \param[in] A The pointer to \a A, an object of type \type Matrix.
 * \param[in] piv The pointer to \a piv, a vector of permutations on \a I, of type \type int.
 * \param[in] rows The binary parameter to indicate whether the rows or the columns of I
 * 		should be permuted (= 0, the rows; = 1, the columns).
 * 
 */
void Matrix::mmCaAIbC(double alpha, double beta, const Matrix &A, int *piv, bool rows)
{
		// Matrix Id( A._cols, _cols, dimT() );			// Id[p][n][nTcoeff]
		// Id.set2Id();										// set to the identity
		// if( rows )											// permute its columns
		// 	Id.cpermutem( piv );
		// else												// permute its rows
		// 	Id.rpermutem( piv );
		// mmCaABbC( alpha, beta, A, Id );
}

//
// S O L V I N G   S Y S T E M S   O F   E Q U A T I O N S
//

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
	
// 	for( int k = 0; k < B._cols; k++ )
// 	{
// 		B( _rows - 1, k ) /= _data[ _rows - 1 ][ _rows - 1 ];
// 		for( int i = _rows - 2; i > -1; i-- )
// 		{
// 				sum.set2zero();								// p(x) = 0, initialization
// 			for( int j = i + 1; j < _rows; j++ )
// 				sum += _data[i][j] * B._data[j][k];
// 			B._data[i][k] = ( B._data[i][k] - sum ) / _data[ i ][ i ];
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
	
// 	for( int k = 0; k < B._cols; k++ )
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
// 			sum += _data[i][j] * b[ j ];
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
// 	for( int k = 0; k < B._rows; k++ )
// 	{
// 		B._data[k][0] /= _data[ 0 ][ 0 ];
// 		for( int i = 1; i < B._cols; i++ )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0, initialization
// 			//else sum = 0.0;
// 			for( int j = 0; j < i; j++ )
// 				sum += _data[ j ][ i ] * B._data[k][j];
// 			B._data[k][i] = ( B._data[k][i] - sum ) / _data[ i ][ i ];
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
// 	Matrix X( B._rows, B._cols, B.dimT() );
  
// 	try
// 	{	
// 		scale = new double[ _rows ];						// scale[m], scale factor
// 		index = new int[ _rows ];							// index[m], row index list
// 		if( (scale == 0)  ||  (scale == NULL)  ||
// 			(index == 0)  ||  (index == NULL) )
// 			throw IDException( "Memory allocation failure.", 3 );
// 		for( int i = 0; i < _rows; i++ )						// determine scale factor, for each row
// 		{
// 			index[ i ] = i;									// initialize row index list
// 			for( int j = 0; j < _rows; j++ )
// 				if( scalemax < fabs( _data[i][j].feval() ) )// scale[i] = max(|A[i][j][0]|)
// 					scalemax = fabs( _data[i][j].feval() );
// 			scale[ i ] = scalemax;
// 		}
// 		for( int k = 0; k < _rows - 1; k++ )					// select pivot row...
// 		{													// ...from max(|A[i][k]/scale[i]|)
// 			jpiv = k;
// 			ratiomax = 0.0;
// 			for( int i = k; i < _rows; i++ )
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
// 			for( int i = k + 1; i < _rows; i++ )				// forward elimination
// 			{
// 				coeff = _data[ index[ i ] ][ k ] / _data[ jindex ][ k ];
// 				for( int j = k + 1; j < _rows; j++ )
// 					_data[ index[ i ] ][ j ] -= coeff * _data[ jindex ][ j ];
// 				_data[ index[ i ] ][ k ] = coeff;
// 				for( h = 0; h < B._cols; h++ )			// for all columns in B
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
// 		for( int i = 0; i < _rows; i++ )						// determine scale factor, for each row
// 		{
// 			index[ i ] = i;									// initialize row index list
// 			for( int j = 0; j < _rows; j++ )
// 				if( scalemax < fabs( _data[i][j].feval() ) )// scale[i] = max(|A[i][j][0]|)
// 					scalemax = fabs( _data[i][j].feval() );
// 			scale[ i ] = scalemax;
// 		}
// 		for( int k = 0; k < _rows - 1; k++ )					// select pivot row...
// 		{													// ...from max(|A[i][k]/scale[i]|)
// 			jpiv = k;
// 			for( int i = k; i < _rows; i++ )
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
// 			for( int i = k + 1; i < _rows; i++ )					// forward elimination
// 			{
// 				coeff = _data[ index[ i ] ][ k ] / _data[ jindex ][ k ];
// 				for( int j = k + 1; j < _rows; j++ )
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
// 			if( strcmp( _data[i][j].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 															// that both strings are equal
// 				sum += _data[i][j].feval() * _data[i][j].feval();
// 			//else sum += _data[i][j] * _data[i][j];	// e.g. a matrix of doubles
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
// 				if( strcmp( _data[i][j].typeName(), "TPolyn" ) == 0 ) // A zero value indicates 
// 																// that both strings are equal
// 					sum += _data[i][j][ k ] * _data[i][j][ k ];
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
// 				c[ j ] += _data[i][j].feval() * _data[i][j].feval();
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

/**
 * Permutes the columns of a matrix given a vector of permutations. 
 * 
 * For example, in case a matrix A is permuted after a QR decomposition with column pivoting,
 * then the resulting matrix in the upper triangular matrix R. 
 * 
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
 * \param[in] trans The \a boolean parameter to indicate whether to transpose the vector
 * 		of permutations \a piv or not (=1, transpose; =0, otherwise). Default is false.
 * 
 */
void Matrix::cpermutem(int *piv, bool trans)
{
	for( int j = 0; j < _cols; j++ )
	{
		if (piv[j] < 0 || piv[j] >= _cols)
		{
			throw MathException("Invalid row permutation."); // TODO add error code
		}
	}

	double **newData;
	Matrix::allocateMemory(newData, _rows, _cols, false);

	int *pivT = 0;
	if (trans)
	{
		pivT = new int[_cols];
		for( int i = 0; i < _cols; i++ )
		{
			pivT[ piv[i] ] = i;
		}
		piv = pivT;
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			newData[i][j] = _data[i][piv[j]];
		}
	}

	delete pivT;

	Matrix::deallocateMemory(_data, _rows, _cols);
	_data = newData;
}

/**
 * Permutes the rows of a matrix given a vector of permutations. 
 * 
 * \param[in] piv The pointer to \a piv, a vector of permutations on the rows of \a A.
 * 
 */
void Matrix::rpermutem(int *piv)
{
	for( int i = 0; i < _rows; i++ )
	{
		if (piv[i] < 0 || piv[i] >= _rows)
		{
			throw MathException("Invalid row permutation."); // TODO add error code
		}
	}

	double **newData;
	Matrix::allocateMemory(newData, _rows, _cols, false);

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			// mein Weg (bei Monett ausgeklamert
			newData[i][j] = _data[ piv[i] ][j];
			
			// Monetts Weg
			// newData[ piv[i] ][j] = _data[i][j];
		}
	}

	Matrix::deallocateMemory(_data, _rows, _cols);
	_data = newData;
}

/**
 * Transposes this matrix in place.
 *
 */
void Matrix::transpose()
{
	if (_rows == _cols)
	{
		// square matrix can be swapped in place
		for( int i = 1; i < _rows; i++ )
		{
			for( int j = 0; j < i; j++ )
			{
				// std::swap(_data[i][j], _data[j][i]);
				double h = _data[i][j];
				_data[i][j] = _data[j][i];
				_data[j][i] = h;
			}
		}
	}
	else
	{
		// size changes, need to allocate new memomry :(
		double **newData;
		Matrix::allocateMemory(newData, _cols, _rows, false);

		for( int i = 0; i < _rows; i++ )
		{
			for( int j = 0; j < _cols; j++ )
			{
				newData[j][i] = _data[i][j];
			}
		}

		Matrix::deallocateMemory(_data, _rows, _cols);
		_data = newData;
		
		// number of rows and columns swapped
		int r = _rows;
		_rows = _cols;
		_cols = r;
	}
}

/**
 * Creata the transposes of this matrix. This matrix opject remains unchanged. 
 * 
 * \return The transposed matrix object.
 *
 */
Matrix Matrix::asTranspose() const
{
 	Matrix aux(_cols, _rows);

 	for( int i = 0; i < _rows; i++ )
 	{
 		for( int j = 0; j < _cols; j++ )
		{
			aux._data[j][i] = _data[i][j];
		}	
 	}

 	return aux;
}

/**
 * Implements the shift operator to calculate the derivative of Taylor polynomials
 * in case the elements of the matrix are such, like in:
 * 
 * y(t) = sum_{j=0}^{d} y_j * t^j + O(t^d+1)
 * 		= y_0 + y_1*t + y_2*t^2 + ... + y_d*t^d
 * 
 * y'(t) = y_1 + 2*y_2*t + 3*y_3*t^2 + ... + d*y_d*t^d-1
 * 
 * Internally, the coefficients are shifted to the left and the last one is zeroed.
 * 
 * \return The error code.
 * 
 */
void Matrix::shift()
{	
	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			_data[i][j] *= 10.0;
			// _data[i][j].shift();
		}
	}
}

/**
 * Returns \a true in case the given matrix is the identity matrix; \a false otherwise.
 * 
 * \return \a true if the matrix is the identity matrix; \a false otherwise.
 * 
 */
bool Matrix::isId() const
{
	// identity matrices must be also sqare matrices
	if (_rows != _cols)
	{
		return false;
	}

	for( int i = 0; i < _rows; i++ )
	{
		for ( int j = 0; j < _cols; j++ )
		{
			if (i == j)
			{
				if (_data[i][j] != 1.0)
				// if (!_data[i][j].isId())
				{
					return false;
				}
			}
			else
			{
				if (_data[i][j] != 0.0)
				// if (!_data[i][j].isZero())
				{
					return false;
				}
			}
		}
	}

	return true;
}

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

/**
 * Returns \a true in case the given matrix is the zero matrix; \a false otherwise.
 * 
 * \return \a true if the matrix is the zero matrix; \a false otherwise.
 * 
 */
bool Matrix::isZero() const
{
	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			if (_data[i][j] != 0.0)
			// if(!_data[i][j].isZero())
			{
				return false;
			}
		}
	}

	return true;
}

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
// 			if( !_data[i][j].isZero( eps ) )			// at least for one coeff. |p(x)| > eps
// 			{
// 				z = false;
// 				break;
// 			}
// 		if( !z )
// 			break;
// 	}
// 	return z;
// }

/**
 * Sets a matrix to the identity one:
 * 
 * 		M = I
 * 
 */
void Matrix::set2Id()
{
	if (_rows != _cols)
	{
		throw MathException("Only square matrices can be set to an identity matrix."); // TODO add error code
	}

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++ )
		{
			if (i == j)
			{
				_data[i][j] = 1.0;
				// _data[i][j].set2Id();
			}
			else
			{
				_data[i][j] = 0.0;
				// _data[i][j].set2zero();
			}
		}
	}
}


/**
 * Sets a submatrix to the identity one:
 * 
 * 		e.g. M = (    | 1 0 0 |    )
 * 				 ( M1 | 0 1 0 | M2 )
 * 				 (    | 0 0 1 |    )
 * 				 (        M3       )
 * 
 * \param[in] top The number of rows at the top to keep unchanged.
 * \param[in] bottom The number of rows at the bottom to keep unchanged.
 * \param[in] left The number of columns on the left to keep unchanged.
 * \param[in] right The number of columns on the right to keep unchanged.
 * 
 */
void Matrix::set2Id(int top, int bottom, int left, int right)
{
	int lastRow = _rows - bottom - 1;
	int lastCol = _cols - right - 1;
	set2IdFromIndices(top, lastRow, left, lastCol);
}

/**
 * Sets a submatrix to the identity one:
 * 
 * 		e.g. M = (    | 1 0 0 |    )
 * 				 ( M1 | 0 1 0 | M2 )
 * 				 (    | 0 0 1 |    )
 * 				 (        M3       )
 * 
 * \param[in] firstRow The row from which to start on.
 * \param[in] lastRow The last row that should be considered.
 * \param[in] firstCol The column from which to start on.
 * \param[in] lastCol The last column that should be considered.
 * 
 */
void Matrix::set2IdFromIndices(int firstRow, int lastRow, int firstCol, int lastCol)
{
	if (firstRow < 0 || lastRow >= _rows)
	{
		throw MathException("Error in set2Id, row bounds are wrong."); // TODO add error code
	}
	if (firstCol < 0 || lastCol >= _cols)
	{
		throw MathException("Error in set2Id, column bounds are wrong"); // TODO add error code
	}

	int rows = lastRow - firstRow + 1;
	int cols = lastCol - firstCol + 1;

	if (rows != cols)
	{
		throw MathException("Error in set2Id, submatrix must be a square matrix in order to become a identity matrix."); // TODO add error code
	}

	for( int i = firstRow; i <= lastRow; i++ )
	{
		for( int j = firstCol; j <= lastCol; j++ )
		{
			if (i - firstRow == j - firstCol)
			{
				_data[i][j] = 1.0;
				// _data[i][j].set2Id();
			}
			else
			{
				_data[i][j] = 0.0;
				// _data[i][j].set2zero();
			}
		}
	}
}

/**
 * Sets a matrix to zero entries.
 * 
 * \return The error code.
 * 
 */
void Matrix::set2Zero()
{
	set2Val(0.0);
}

/**
 * Sets a submatrix to zero:
 * 
 * 		e.g. M = (    | 0 0 0 |    )
 * 				 ( M1 | 0 0 0 | M2 )
 * 				 (    | 0 0 0 |    )
 * 				 (        M3       )
 * 
 * \param[in] top The number of rows at the top to keep unchanged.
 * \param[in] bottom The number of rows at the bottom to keep unchanged.
 * \param[in] left The number of columns on the left to keep unchanged.
 * \param[in] right The number of columns on the right to keep unchanged.
 * 
 */
void Matrix::set2Zero(int top, int bottom, int left, int right)
{
	int lastRow = _rows - bottom - 1;
	int lastCol = _cols - right - 1;
	set2ZeroFromIndices(top, lastRow, left, lastCol);
}

/**
 * Sets a submatrix to zero:
 * 
 * 		e.g. M = (    | 0 0 0 |    )
 * 				 ( M1 | 0 0 0 | M2 )
 * 				 (    | 0 0 0 |    )
 * 				 (        M3       )
 * 
 * \param[in] firstRow The row from which to start on.
 * \param[in] lastRow The last row that should be considered.
 * \param[in] firstCol The column from which to start on.
 * \param[in] lastCol The last column that should be considered.
 * 
 */
void Matrix::set2ZeroFromIndices(int firstRow, int lastRow, int firstCol, int lastCol)
{
	set2ValFromIndices(firstRow, lastRow, firstCol, lastCol, 0.0);
}

/**
 * Sets a matrix to the value given as parameter.
 * 
 * \param[in] v The double value to set the elements to.
 * 
 */
void Matrix::set2Val(double v)
{
	set2ValFromIndices(0, _rows - 1, 0, _cols - 1, v);
}

/**
 * Sets a submatrix to the value given as parameter:
 * 
 * 		e.g. M = (    | v v v |    )
 * 				 ( M1 | v v v | M2 )
 * 				 (    | v v v |    )
 * 				 (        M3       )
 * 
 * \param[in] top The number of rows at the top to keep unchanged.
 * \param[in] bottom The number of rows at the bottom to keep unchanged.
 * \param[in] left The number of columns on the left to keep unchanged.
 * \param[in] right The number of columns on the right to keep unchanged.
 * \param[in] v The double value to set the elements to.
 * 
 */
void Matrix::set2Val(int top, int bottom, int left, int right, double v)
{
	int lastRow = _rows - bottom - 1;
	int lastCol = _cols - right - 1;
	set2ValFromIndices(top, lastRow, left, lastCol, v);
}

/**
 * Sets a submatrix to the value given as parameter:
 * 
 * 		e.g. M = (    | v v v |    )
 * 				 ( M1 | v v v | M2 )
 * 				 (    | v v v |    )
 * 				 (        M3       )
 * 
 * \param[in] firstRow The row from which to start on.
 * \param[in] lastRow The last row that should be considered.
 * \param[in] firstCol The column from which to start on.
 * \param[in] lastCol The last column that should be considered.
 * \param[in] v The double value to set the elements to.
 * 
 */
void Matrix::set2ValFromIndices(int firstRow, int lastRow, int firstCol, int lastCol, double v)
{
	if (firstRow < 0 || lastRow >= _rows)
	{
		throw MathException("Invalid row bounds."); // TODO add error code
	}
	if (firstCol < 0 || lastCol >= _cols)
	{
		throw MathException("Invalid column bounds."); // TODO add error code
	}

	if (v == 0.0)
	{
		for( int i = firstRow; i <= lastRow; i++ )
		{
			for( int j = firstCol; j <= lastCol; j++ )
			{
				_data[i][j] = 0.0;
				// _data[i][j].set2Zero();
			}
		}
	}
	else
	{
		for( int i = firstRow; i <= lastRow; i++ )
		{
			for( int j = firstCol; j <= lastCol; j++ )
			{
				_data[i][j] = v;
				// _data[i][j].set2Const(v);
			}
		}
	}
}

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
// 	for( int i = 0; i < _rows; i++ )
// 	{
// 		Am( i, i ) = p / _data[ i ][ i ];					// elements of the diagonal
// 	}
// 	for( int j = _rows - 2; j > -1; j-- )
// 		for( int i = j + 1; i < _rows; i++ )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0
// 			//else sum = 0.0;
// 			for( int k = j + 1; k < _rows; k++ )
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
// 	for( int i = 0; i < r; i++ )
// 		Am( i, i ) = p / _data[ i ][ i ];					// elements of the diagonal
// 	for( int j = r - 2; j > -1; j-- )
// 		for( int i = j + 1; i < r; i++ )
// 		{
// 			if( strcmp( sum.typeName(), "TPolyn" ) == 0 ) 	// A zero value indicates 
// 															// that both strings are equal
// 				sum.set2zero();								// p(x) = 0
// 			//else sum = 0.0;
// 			for( int k = j + 1; k < r; k++ )
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
// 			if( _data[i][j] != B._data[i][j] )
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
// 			diff = _data[i][j] - B._data[i][j];
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
// 			diff = _data[i][j] - B._data[i][j];
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
 * Prints out the matrix
 */
void Matrix::print()
{
 	for( int i = 0; i < _rows; i++ )
 	{
 		for( int j = 0; j < _cols; j++ )
 		{
 			printf("%3.2f\t", _data[i][j]);
 		}
 		printf("\n");
 	}
}

/**
 * Prints out the matrix after printing the given string
 *
 * \param[in] before The character string to be printed out before the matrix.
 * 
 */
void Matrix::print(const char *name)
{
 	printf("%s (%dx%d):\n", name, _rows, _cols);
 	print();
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
// 			_data[i][j].print( fn );
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
// 			_data[i][j].print();
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
// 			_data[i][j].print();
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
// 			_data[i][j].print( fn );
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
// 			_data[i][j].print();
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
// 			if( _data[i][j].isZero( eps ) )				//small enough
// 				printf( "%.16g%c", 0.0, '\t' );
// 			else
// 			{
// 				_data[i][j].print();
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
// 			if( _data[i][j].isZero( eps ) )				//small enough
// 				fprintf( fn, "%.16g%c", 0.0, '\t' );
// 			else
// 			{
// 				_data[i][j].print( fn );
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
// 				if( abs( _data[i][j][ k ] ) < eps )		//small enough
// 					printf( "%.16lg%c", 0.0, '\t' );
// 				else
// 					printf( "%.16lg%c", _data[i][j][ k ], '\t' );
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
// 				if( abs( _data[i][j][ k ] ) < eps )		//small enough
// 					fprintf( fn, "%.16lg%c", 0.0, '\t' );
// 				else
// 					fprintf( fn, "%.16lg%c", _data[i][j][ k ], '\t' );
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
// 				if( abs( _data[i][j][ k ] ) < eps )		//small enough
// 					printf( "%.16lg%c", 0.0, '\t' );
// 				else
// 					printf( "%.16lg%c", _data[i][j][ k ], '\t' );
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
// 				if( abs( _data[i][j][ k ] ) < eps )		//small enough
// 					fprintf( fn, "%.16lg%c", 0.0, '\t' );
// 				else
// 					fprintf( fn, "%.16lg%c", _data[i][j][ k ], '\t' );
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

void Matrix::allocateMemory(double **&data, int rows, int cols, bool initialize)
{
	try
	{
		// bounds checking
		if (rows <= 0 || cols <= 0)
		{
			throw MathException("Matrix invalid matrix size");
		}

		data = new double*[rows];
		if (data == 0 || data == NULL)
		{
			throw MathException("Memory allocation failure.");
		}

		for( int i = 0; i < rows; i++ )
		{
			data[i] = new double[cols];
			if (data[i] == 0 || data[i] == NULL)
			{
				throw MathException("Memory allocation failure.");
			}

			if (initialize)
			{
				for( int j = 0; j < cols; j++)
				{
					data[i][j] = 0.0;
					// data[i][j] = TPoly(_dimT);
				}
			}
		}
	}
	catch(bad_alloc e)
	{
		throw MathException(e.what(), 4);
	}
	catch(...)
	{
		throw MathException("Error when allocating a matrix.");
	}
}

void Matrix::deallocateMemory(double **&data, int rows, int cols)
{
	for( int i = 0; i < rows; i++ )
	{
		delete [] data[i];
	}
	delete [] data;
}

void Matrix::copyFrom(const Matrix &m)
{
	deallocateMemory();

	_rows = m._rows;
	_cols = m._cols;
	_dimT = m._dimT;

	allocateMemory(false);

	for( int i = 0; i < _rows; i++ )
	{
		for( int j = 0; j < _cols; j++)
		{
			_data[i][j] = m._data[i][j];
		}
	}
}
