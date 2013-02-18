using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LibMatrix
{

    /// <summary>
    /// Matrix Class
    /// </summary>
    public class Matrix
    {
        #region private vars
        
        /// <summary>
        /// Number of rows
        /// </summary>
        private int _rows;
        
        /// <summary>
        /// Number of columns
        /// </summary>
        private int _cols;

        /// <summary>
        /// The dimension of the Taylor polynomials
        /// </summary>
        private int _dimT;

        /// <summary>
        /// The Taylor Polynomials in an 2 dimensional Matrix Array
        /// </summary>
        private Polynomial[,] _data;

        #endregion

        #region public getter and setter 

        public int nrows() { return _rows; }
        public int ncols() { return _cols; }
        public int dimT() { return _dimT; }
        public Polynomial get(int row, int col)
        {
            if (row >= _rows || row < 0 || col >= _cols || col < 0)
                throw new MathException(String.Format("({0:d}, {1:d}) is not a alid index of {2:d} x {3:d} matrix", row, col, _rows, _cols));

            return _data[row, col];
        }

        #endregion

        #region ctors

        /// <summary>
        /// Default constructor for the class. Creates the object.<para/>
        /// <para/>
        /// It is a 1-by-1 matrix, i.e., it has only one element, which is set to zero:<para/>
        /// <para/>
        /// 		m(0,0) == 0.0<para/>
        /// </summary>
        public Matrix() : this(1,1,0)
        {
        }

        /// <summary>
        /// Constructor for the class with both the number of rows and columns as parameters.<para/> 
        /// Creates the object.
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="cols">The number of columns.</param>
        public Matrix(int rows, int cols) : this(rows,cols,0)
        {
        }

        /// <summary>
        /// Constructor for the class with both the number of rows and columns as parameters,<para/> 
        /// as well as the dimension of the elements' type, e.g., the Taylor polynomial's grade. <para/> 
        /// Creates the object.<para/> 
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="cols">The number of columns.</param>
        /// <param name="dimT">The dimension</param>
		public Matrix(int rows, int cols, int dimT)
        {
            initializeMatrix(rows, cols, dimT, initializePolynoms(rows, cols, dimT));
        }

        /// <summary>
        /// Copy constructor.<para/> 
        /// <para/> 
        /// 		Matrix *newm = new Matrix( (*m) );<para/> 
        /// </summary>
        /// <param name="matrix"> A Matrix object to copy from.</param>
		public Matrix(Matrix matrix)
        {
            initializeMatrix(matrix._rows, matrix._cols, matrix._dimT, matrix._data);
        }

        /// <summary>
        /// Easy constructor for testing.
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="cols">The number of columns.</param>
        /// <param name="values">an initialised Polynom</param>
        public Matrix(int rows, int cols, Polynomial[,] values)
        {
            int dimt = 0;

            if (values[0, 0] != null)
                dimt = values[0, 0].order;

            initializeMatrix(rows, cols, dimt, values);
        }

        /// <summary>
        /// Destructor. Cleans up the object.
        /// </summary>
        ~Matrix()
        {
            _rows = 1;
            _cols = 1;
            _dimT = 0;
            _data = new Polynomial[1,1];
        }
        
        #endregion

        #region public functions

        #region operators

        /// <summary>
        /// Implements the [] operator.
        /// </summary>
        /// <param name="row">The row index</param>
        /// <param name="col">The column index</param>
        /// <returns>The value of the desired element.</returns>
        public Polynomial this[int row, int col]
        {
            get
            {
                if (row >= _rows || row < 0 || col >= _cols || col < 0)
                    throw new MathException(String.Format("({0:d}, {1:d}) is not a alid index of {2:d} x {3:d} matrix", row, col, _rows, _cols));

                return _data[row, col];
            }
            set
            {
                if (row >= _rows || row < 0 || col >= _cols || col < 0)
                    throw new MathException(String.Format("({0:d}, {1:d}) is not a alid index of {2:d} x {3:d} matrix", row, col, _rows, _cols));

                _data[row, col] = value;
            }
        }

        /// <summary>
        /// Implements the == operator. It compares two matrices.
        /// </summary>
        /// <param name="a">Matrix on the left side</param>
        /// <param name="b">Matrix on the right side</param>
        /// <returns>true if the matrices are equal. Otherwise it returns false.</returns>
        public static bool operator ==(Matrix a, Matrix b)
        {
            if (a._rows != b._rows || a._cols != b._cols || a._dimT != b._dimT) // TODO _dimT ?
            {
                return false;
            }

            for (int i = 0; i < a._rows; i++)
            {
                for (int j = 0; j < a._cols; j++)
                {
                    if (a._data[i,j] != b._data[i,j])
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// Implements the != operator.<para/>
        /// Compares two matrices.<para/>
        /// </summary>
        /// <param name="a">Matrix on the left side</param>
        /// <param name="b">Matrix on the right side</param>
        /// <returns>true if not equal</returns>
        public static bool operator !=(Matrix a, Matrix b) { return !(a == b); }

        /// <summary>
        /// Implements the + operator. It adds up two matrices.
        /// </summary>
        /// <param name="a">Matrix on the left side</param>
        /// <param name="b">Matrix on the right side</param>
        /// <returns>the resulting Matrix object.</returns>
        public static Matrix operator +(Matrix a, Matrix b)
        {
	        // dimensions checking
	        if(a._rows != b._rows || a._cols != b._cols)
	        {
                throw new MathException(String.Format("Cannot add up a {0:d} x {1:d} matrix (left side) and a {2:d} x {3:d} (right side) matrix.", a._rows, a._cols, b._rows, b._cols));
	        }

	        // an auxiliary object
	        Matrix aux = new Matrix(a._rows, a._cols, a._dimT);

	        for( int i = 0; i < a._rows; i++ )
	        {
		        for( int j = 0; j < a._cols; j++ )
		        {
			        aux._data[i,j] = a._data[i,j] + b._data[i,j];
		        }
	        }

	        return aux;
        }

        /// <summary>
        /// Implements the - operator. It substracts two matrices.
        /// </summary>
        /// <param name="a">Matrix on the left side</param>
        /// <param name="b">Matrix on the right side</param>
        /// <returns>the resulting Matrix object.</returns>
        public static Matrix operator -(Matrix a, Matrix b)
        {
            // dimensions checking
            if (a._rows != b._rows || a._cols != b._cols)
            {
                throw new MathException(String.Format("Cannot add up a {0:d} x {1:d} matrix (left side) and a {2:d} x {3:d} (right side) matrix.", a._rows, a._cols, b._rows, b._cols));
            }

            // an auxiliary object
            Matrix aux = new Matrix(a._rows, a._cols, a._dimT);

            for (int i = 0; i < a._rows; i++)
            {
                for (int j = 0; j < a._cols; j++)
                {
                    aux._data[i, j] = a._data[i, j] - b._data[i, j];
                }
            }

            return aux;
        }

        /// <summary>
        /// Implements the unary - operator.
        /// </summary>
        /// <param name="a">Matrix to operate on</param>
        /// <returns>the resulting Matrix object.</returns>
        public static Matrix operator -(Matrix a)
        {
            Matrix retval = new Matrix(a._rows, a._cols, a._dimT);
            for (int i = 0; i < a._rows; i++)
                for (int j = 0; j < a._cols; j++)
                    retval[i, j] = -a[i, j];

            return retval;
        }

        /// <summary>
        /// Implements the * operator. It multiplies a matrix by a scalar.
        /// </summary>
        /// <param name="a">The Matrix to multiply with.</param>
        /// <param name="alpha">The scalar to multiply by.</param>
        /// <returns>the resulting Matrix object.</returns>
        public static Matrix operator *(Matrix a, double alpha)
        {
        	Matrix aux = new Matrix(a._rows, a._cols, a._dimT);

	        for( int i = 0; i < a._rows; i++ )
	        {
                for (int j = 0; j < a._cols; j++)
		        {
                    aux._data[i, j] = a._data[i,j] * alpha;
		        }
	        }

	        return aux;
        }

        /// <summary>
        /// Implements the * operator. It multiplies the matrix with another matrix.
        /// </summary>
        /// <param name="a">Matrix on the left side</param>
        /// <param name="b">Matrix on the right side</param>
        /// <returns>the resulting Matrix object.</returns>
        public static Matrix operator *(Matrix a, Matrix b)
        {
            if (a._cols != b._rows)
                throw new MathException(String.Format("Cannot multiply a {0:d} x {1:d} with a {2:d} x {3:d} matrix. The number of rows in m must match the number of colmuns in this matrix.", a._rows, a._cols, b._rows, b._cols));

            Matrix aux = new Matrix(a._rows, b._cols, a._dimT);

            for (int i = 0; i < a._rows; i++)
                for (int j = 0; j < b._cols; j++)
                    for (int k = 0; k < a._cols; k++)
                    {
                        aux._data[i, j] += a._data[i, k] * b._data[k, j];
                    }

            return aux;
        }

        #endregion

        #region functions

        #region E S P E C I A L   M A T R I X   M U L T I P L I C A T I O N S

        /// <summary>
        /// Matrix multiplication of the form:<para/>
        /// <para/>
        ///     C = alpha * A * B + beta * C<para/>
        /// <para/>
        /// with A : m-by-p matrix<para/>
        /// 	 B : p-by-n matrix<para/>
        /// 	 C : m-by-n matrix<para/>
        /// 	 alpha, beta : real numbers<para/>
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * B</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix</param>
        /// <param name="B">an object of type Matrix</param>
        public void mmCaABbC(double alpha, double beta, Matrix A, Matrix B)
        {
            if (A._cols != B._rows)
            {
                throw new MathException(String.Format("Cannot multiply A ({0:d} x {1:d}) und B ({2:d} x {3:d}).", A._rows, A._cols, B._rows, B._cols));
            }
            if (A._rows != _rows || B._cols != _cols)
            {
                throw new MathException(String.Format("The size of the matrix A * B ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._rows, B._cols, _rows, _cols));
            }

            for( int i = 0; i < _rows; i++ )
	        {
		        for( int j = 0; j < _cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);
			
			        for( int k = 0; k < B._rows; k++ )
			        {
				        h += A._data[i,k] * B._data[k,j];
			        }

			        _data[i,j] = h * alpha + _data[i,j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:<para/>
        /// <para/>
        ///     C = alpha * A * B + beta * C <para/>
        /// <para/>
        /// with A : m-by-p matrix<para/>
        /// 	 B : p-by-n matrix<para/>
        ///      C : m-by-n matrix<para/>
        ///      alpha, beta : real numbers<para/>
        /// <para/>
        /// where the inferior-right block of B is an identity matrix like in:<para/>
        ///     ( * * * 0 0 )<para/>
        ///     ( * * * 0 0 )<para/>
        ///     ( 0 0 0 1 0 )<para/>
        ///     ( 0 0 0 0 1 )<para/>
        /// <para/>
        /// so that a particular block multiplication is needed.<para/>
        /// </summary>
        /// <param name="r">The number of rows in B that are of interest (2 in the example above)</param>
        /// <param name="c">The number of columns in B that are of interest (3 in the example above)</param>
        /// <param name="alpha">The scalar value that multiplies A * B</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix</param>
        /// <param name="B">an object of type Matrix</param>
        public void bmmCaABbC(int r, int c, double alpha, double beta, Matrix A, Matrix B)
        {
            if (r > B._rows)
            {
                throw new MathException(String.Format("The parameter r (= {0:d}) must be smaller or equal than the number of rows in B ({1:d} x {2:d}).", r, B._rows, B._cols));
            }
            if (c > B._cols)
            {
                throw new MathException(String.Format("The parameter c (= {0:d}) must be smaller or equal than the number of columns in B ({1:d} x {2:d}).", c, B._rows, B._cols));
            }
            if (A._cols != B._rows)
            {
                throw new MathException(String.Format("Cannot multiply A ({0:d} x {1:d}) und B ({2:d} x {3:d}).", A._rows, A._cols, B._rows, B._cols));
            }
            if (A._rows != _rows || B._cols != _cols)
            {
                throw new MathException(String.Format("The size of the matrix A * B ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._rows, B._cols, _rows, _cols));
            }

            int cr = c - r;
            // If you wonder why cr = c-r and latter A(i, j-cr)?
            // I don't know, bit take a sheet of paper the the matrices and you will see ;)

            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    Polynomial h = new Polynomial(_dimT);

                    for (int k = 0; k < r; k++)
                    {
                        h += A._data[i, k] * B._data[k, j];
                    }
                    _data[i, j] = h * alpha + _data[i, j] * beta;
                }
            }

            for (int i = r; i < _rows; i++)
            {
                for (int j = c; j < _cols; j++)
                {
                    _data[i, j] = A._data[i, j - cr] * B._data[j - cr, i] * alpha + _data[i, j] * beta;
                }
            }
        }

        /// <summary>
        /// Matrix multiplication of the form:<para/>
        /// <para/>
        ///     C = alpha * A * B + beta * C<para/>
        /// <para/>
        /// with A : m-by-p matrix<para/>
        /// 	 B : p-by-n matrix<para/>
        ///      C : r-by-n matrix (only the last r rows from A are interesting)<para/>
        ///      alpha, beta : real numbers<para/>
        /// <para/>
        /// where A ("special" A) is of the form:<para/>
        /// <para/>
        ///     (         )<para/>
        ///     (    X    )<para/>
        ///     ( ------- )<para/>
        ///     ( * ... * )<para/>
        ///     ( * ... * )<para/>
        /// <para/>
        /// so that a particular matrix multiplication is needed.<para/>
        /// </summary>
        /// <param name="r">The last rows in A that are of interest (2 non-zero-rows in the example above)</param>
        /// <param name="alpha">The scalar value that multiplies A * B</param>
        /// <param name="beta">beta The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix</param>
        /// <param name="B">an object of type Matrix</param>
        public void mmCasABbC(int r, double alpha, double beta, Matrix A, Matrix B)
        {
	        if (r > A._rows)
	        {
                throw new MathException(String.Format("The number of rows in A that are of interest (r = {0:d}) must be smaller or equal than the total amount of rows in A ({1:d} x {2:d}).", r, A._rows, A._cols));
	        }
	        if (A._cols != B._rows)
	        {
                throw new MathException(String.Format("Cannot multiply A ({0:d} x {1:d}) und B ({2:d} x {3:d}).", A._rows, A._cols, B._rows, B._cols));
	        }
	        if (A._rows != _rows || B._cols != _cols)
	        {
                throw new MathException(String.Format("The size of the matrix A * B ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._rows, B._cols, _rows, _cols));
	        }

	        int n = A._rows - r;
	        for( int i = 0; i < n; i++ )
	        {
		        for( int j = 0; j < B._cols; j++ )
		        {
			        _data[i,j] *= beta;
		        }
	        }

	        for( int i = n; i < A._rows; i++ )
	        {
		        for( int j = 0; j < B._cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < B._rows; k++ )
			        {
				        h += A._data[i,k] * B._data[k,j];
			        }

			        _data[i,j] = h * alpha + _data[i,j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:<para/>
        ///<para/>
        ///     C = alpha * A * B + beta * C<para/>
        /// <para/>
        /// with A : m-by-p matrix<para/>
        /// 		B : p-by-n matrix<para/>
        /// 		C : m-by-r matrix (only the last r columns from B are interesting).<para/>
        /// 		alpha, beta : real numbers<para/>
        /// <para/>
        /// where B ("special" B) is of the form:<para/>
        /// <para/>
        /// 		(   | * * * )<para/>
        /// 		(   | . . . )<para/>
        /// 		( X | . . . )<para/>
        /// 		(   | . . . )<para/>
        /// 		(   | * * * )<para/>
        /// <para/>
        /// so that a particular matrix multiplication is needed.<para/>
        /// </summary>
        /// <param name="r">The last columns in B that are of interest (3 non-zero-columns - * in the example above)</param>
        /// <param name="alpha">The scalar value that multiplies A * B</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix</param>
        /// <param name="B">an object of type Matrix</param>
        public void mmCaAsBbC(int r, double alpha, double beta, Matrix A, Matrix B)
        {
            if (r > B._cols)
	        {
                throw new MathException(String.Format("The number of columns in B that are of interest (r = {0:d}) must be smaller or equal than the total amount of columns in B ({1:d} x {2:d}).", r, A._rows, A._cols));
	        }
	        if (A._cols != B._rows)	
	        {
                throw new MathException(String.Format("Cannot multiply A ({0:d} x {1:d}) und B ({2:d} x {3:d}).", A._rows, A._cols, B._rows, B._cols));
	        }
	        if (A._rows != _rows || B._cols != _cols)
	        {
                throw new MathException(String.Format("The size of the matrix A * B ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._rows, B._cols, _rows, _cols));
	        }

	        int n = B._cols - r;
	        for( int i = 0; i < _rows; i++ )
	        {
		        // zero-columns of B
		        for( int j = 0; j < n; j++ )
		        {
			        _data[i,j] *= beta;
		        }

		        // non-zero-columns of B
		        for( int j = n; j < _cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < B._rows; k++ )
			        {
				        h += A._data[i,k] * B._data[k,j];
			        }

			        _data[i,j] = h * alpha + _data[i,j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:<para/>
        ///<para/>
        ///     C = alpha * A * UTB + beta * C<para/>
        /// <para/>
        /// where UTB means that only the upper triangular part is of interest. Furthermore, a column<para/>
        /// pivoting on B is considered.<para/>
        /// <para/>
        /// with A : m-by-p matrix<para/>
        /// 		B : p-by-n matrix<para/>
        /// 		C : m-by-r matrix (only the last r columns from B are interesting).<para/>
        /// 		alpha, beta : real numbers<para/>
        /// <para/>
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * B</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix</param>
        /// <param name="B">an object of type Matrix</param>
        /// <param name="piv">a vector of permutations on the columns of B</param>
        public void mmCaAUTBPbC(double alpha, double beta, Matrix A, Matrix B, int[] piv)
        {
	        if (A._cols != B._rows)	
	        {
		        throw new MathException(String.Format("Cannot multiply A ({0:d} x {1:d}) und B ({2:d} x {3:d}).", A._rows, A._cols, B._rows, B._cols));
	        }
	        if (A._rows != _rows || B._cols != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A * B ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._rows, B._cols, _rows, _cols));
	        }

	        for( int i = 0; i < _rows; i++ )
	        {
		        for( int j = 0; j < _cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

                    for (int k = 0; k < _cols; k++)
			        {
				        h += A._data[i,k] * B._data[k,piv[j]];
			        }

			        _data[i,j] += h * alpha + _data[i,j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A * A^T + beta * C
        /// 
        /// with A, C : m-by-m matrix
        ///      alpha, beta : real numbers
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * A^T</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix. Its transpose is also considered</param>
        public void mmCaAATbC(double alpha, double beta, Matrix A) 
        { 
            // if A is a m-by-n matrix A * A^T is always a m-by-m matrix and must have the same dimension as this matrix
	        if (A._rows != _rows || A._rows != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A * A^T ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._rows, A._rows, _rows, _cols));
	        }
	
	        for( int i = 0; i < _rows; i++ )
	        {
		        for( int j = 0; j < _rows; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < _rows; k++ )
			        {
                        h += A._data[i, k] * A._data[j, k];
			        }

                    _data[i, j] = h * alpha + _data[i, j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A^T * A + beta * C
        /// 
        /// with A, C : m-by-m matrix
        ///      alpha, beta : real numbers
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * A^T</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix. Its transpose is also considered</param>
        public void mmCaATAbC(double alpha, double beta, Matrix A) 
        { 
           	// if A is a m-by-n matrix A^T * A is always a n-by-n matrix and must have the same dimension as this matrix
	        if (A._cols != _rows || _cols != _rows)
	        {
		        throw new MathException(String.Format("The size of the matrix A^T * A ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._cols, A._cols, _rows, _cols));
	        }

	        for( int i = 0; i < _rows; i++ )
	        {
		        for( int j = 0; j < _rows; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < A._rows; k++ )
			        {
                        h += A._data[k, i] * A._data[k, j];
			        }

                    _data[i, j] = h * alpha + _data[i, j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A^T * B + beta * C
        /// 
        /// with A : p-by-m matrix
        ///      B : p-by-n matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * A^T</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix. Its transpose is also considered</param>
        /// /// <param name="B">an object of type Matrix. Its transpose is also considered</param>
        public void mmCaATBbC(double alpha, double beta, Matrix A, Matrix B) 
        { 
	        // A^T und B can only be multiplied if A and B have the same number of rows
	        if (A._rows != B._rows)
	        {
		        throw new MathException(String.Format("Cannot multiply A^T ({0:d} x {1:d}) und B ({2:d} x {3:d}).", A._cols, A._rows, B._rows, B._cols));
	        }

	        // (A^T * B) must have the same size as this matrix
	        if (A._cols != _rows || B._cols != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A^T * B ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._cols, B._cols, _rows, _cols));
	        }

	        for( int i = 0; i < _rows; i++ )
	        {
		        for( int j = 0; j < _cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < B._rows; k++ )
			        {
                        h += A._data[k, i] * B._data[k, j];
			        }

                    _data[i, j] = h * alpha + _data[i, j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A^T * B + beta * C
        /// 
        /// with A : p-by-m matrix
        ///      B : p-by-n matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// and a column pivoting on A^Ts rows
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * A^T</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix. Its transpose is also considered</param>
        /// <param name="B">an object of type Matrix</param>
        /// <param name="piv">a vector of permutations on the columns of B</param>
        public void mmCaATBPbC(double alpha, double beta, Matrix A, Matrix B, int[] piv) 
        { 
            // A^T und B can only be multiplied if A and B have the same number of rows
	        if (A._rows != B._rows)
	        {
		        throw new MathException(String.Format("Cannot multiply A^T ({0:d} x {1:d}) und B ({2:d} x {3:d}).", A._cols, A._rows, B._rows, B._cols));
	        }

	        // (A^T * B) must have the same size as this matrix
	        if (A._cols != _rows || B._cols != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A^T * B ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._cols, B._cols, _rows, _cols));
	        }

	        for( int i = 0; i < _rows; i++ )
	        {
		        for( int j = 0; j < _cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < B._rows; k++ )
			        {
                        h += A._data[k, i] * B._data[k, piv[j]];
			        }

                    _data[i, j] = h * alpha + _data[i, j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A * B^T + beta * C
        /// 
        /// with A : m-by-p matrix
        ///      B : n-by-p matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * B^T</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix.</param>
        /// <param name="B">an object of type Matrix. Its transpose is also considered</param>
        public void mmCaABTbC(double alpha, double beta, Matrix A, Matrix B) 
        { 
            // A und B^T can only be multiplied if A and B have the same number of columns
	        if (A._cols != B._cols)
	        {
		        throw new MathException(String.Format("Cannot multiply A ({0:d} x {1:d}) und B^T ({2:d} x {3:d}).", A._rows, A._cols, B._cols, B._rows));
	        }

	        // (A * B^T) must have the same size as this matrix
	        if (A._rows != _rows || B._rows != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A * B^T ({0:d} x {1:d}) must match the current matrix ({2:d} x {3:d}).", A._rows, B._rows, _rows, _cols));
	        }

	        for( int i = 0; i < _rows; i++ )
	        {
		        for( int j = 0; j < _cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < A._cols; k++ )
			        {
                        h += A._data[i, k] * B._data[j, k];
			        }

                    _data[i, j] = h * alpha + _data[i, j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A * B^T + beta * C
        /// 
        /// with A : m-by-p matrix
        ///      B : n-by-p matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// After transposing B, either its first or its last rows are considered for multiplication,
        /// 
        ///     ( * ... * )				(    X    )
        ///     ( ------- )		or		( ------- )
        ///     (         )				( * ... * )
        ///     (    X    )				( * ... * )
        ///     
        /// according to A dimensions. I.e., the matrix A has less columns than B^T rows has.
        /// </summary>
        /// <param name="r">The number of rows from B that should be considered.</param>
        /// <param name="up">The binary parameter to indicate whether the first or the last r rows</param>
        /// <param name="alpha">The scalar value that multiplies A * B^T</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix.</param>
        /// <param name="B">an object of type Matrix. Its transpose is also considered</param>
        public void mmCaABTbC(int r, bool up, double alpha, double beta, Matrix A, Matrix B) 
        { 
	        // A und B^T can only be multiplied if A and B have the same number of columns
	        if (A._cols != B._cols)
	        {
		        throw new MathException(String.Format("Cannot multiply A ({0:d} x {1:d}) und B^T ({2:d} x {3:d}).", A._rows, A._cols, B._cols, B._rows));
	        }

	        // (A * B^T) must have the same size as this matrix
	        if (A._rows != _rows || B._rows != _cols)
	        {
		        throw new MathException("Error in matrix multiplication. The result of A * B must have the same size as C (this matrix).");
	        }

	        // r <= (colmuns of B^T), which is equal to r <= (rows of B)
	        if (r > B._rows)
	        {
		        throw new MathException("Error in matrix multiplication. r must be smaller or equal than the number of rows in B.");
	        }

	        int bRowsR = B._rows - r;

	        for( int i = 0; i < A._rows; i++ )
	        {
		        for( int j = 0; j < B._cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        if (up)
			        {
				        // first r columns from B are used
				        for( int k = 0; k < r; k++ )
				        {
                            h += A._data[i, k] * B._data[j, k];
				        }
			        }
			        else
			        {
				        // last r columns from B are used
				        for( int k = bRowsR; k < B._rows; k++ )
				        {
                            h += A._data[i, k] * B._data[j, k];
				        }
			        }

                    _data[i, j] = h * alpha + _data[i, j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A * B^T + beta * C
        /// 
        /// with A : m-by-p matrix
        ///      B : n-by-p matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// where the inferior-right block of A is an identity matrix like in:
        ///     
        ///     ( * * * 0 0 )
        ///     ( * * * 0 0 )
        ///     ( 0 0 0 1 0 )
        ///     ( 0 0 0 0 1 )
        ///     
        /// so that a particular block multiplication is needed.
        /// </summary>
        /// <param name="r">The number of rows in A that are of interest (2 in the example above).</param>
        /// <param name="c">The number of columns in A that are of interest (3 in the example above).</param>
        /// <param name="alpha">The scalar value that multiplies A * B^T</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix.</param>
        /// <param name="B">an object of type Matrix. Its transpose is also considered</param>
        public void bmmCaABTbC(int r, int c, double alpha, double beta, Matrix A, Matrix B) 
        { 
            if (r > A._rows)
	        {
		        throw new MathException("Error in matrix multiplication. r cannot be larger than the number of rows of A in total.");
	        }
	        if (c > A._cols)
	        {
		        throw new MathException("Error in matrix multiplication. c cannot be larger than the number of columns of A in total.");
	        }
	        if (A._cols != B._cols)
	        {
		        throw new MathException("Error in matrix multiplication. The matrices A und B^T cannot be multiplied because of wrong dimensions.");
	        }
	        if (A._rows != _rows || B._rows != _cols)
	        {
		        throw new MathException("Error in matrix multiplication. The dimension of the matrix A * B^T must match the current matrix.");
	        }

	        // only the first r rows from A are interesting
	        for( int i = 0; i < r; i++ )
	        {
		        for( int j = 0; j < _cols; j++ )
		        {
			        Polynomial h = new Polynomial(_dimT);

			        for( int k = 0; k < c; k++ )
			        {
                        h += A._data[i, k] * B._data[j, k];
			        }

                    _data[i, j] = h * alpha + _data[i, j] * beta;
		        }
	        }

	        // last rows of A, with exactly one 1 per row
	        int d = A._cols - A._rows;
	        for ( int i = r; i < _rows; i++ )
	        {
		        for( int j = 0; j < _cols; j++ )
		        {
                    _data[i, j] = B._data[j, i + d] * alpha + _data[i, j] * beta;
		        }
	        }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * I * B + beta * C
        /// 
        /// with I : m-by-p matrix; identity matrix
        ///      B : p-by-n matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies I * B</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="B">an object of type Matrix</param>
        public void mmCaIBbC(double alpha, double beta, Matrix B) 
        {
            if (B._rows != _rows || B._cols != _cols)
            {
                throw new MathException("Error in matrix multiplication. The dimension of the matrix B must match the current matrix.");
            }

            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    _data[i, j] = B._data[i, j] * alpha + _data[i, j] * beta;
                }
            }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * I * B + beta * C
        /// 
        /// with I : m-by-p matrix; identity matrix permuted according to a vector of permutations, piv
        ///      B : p-by-n matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies I * B</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="piv">a vector of permutations on I, of type int</param>
        /// <param name="rows">The binary parameter to indicate whether the rows or the columns of I should be permuted (true for the rows; false for the columns)</param>
        /// <param name="B">an object of type Matrix</param>
        public void mmCaIBbC(double alpha, double beta, int[] piv, bool rows, Matrix B) 
        { 
           	//TODO: implement a faster version without identity matrix
	        Matrix Id = new Matrix(_rows, B._rows, _dimT);
	        Id.set2Id();				

	        if (rows)
	        {
		        // permute its columns
		        Id.cpermutem(piv);
	        }
	        else
	        {
		        // permute its rows
		        Id.rpermutem(piv);
	        }
	        mmCaABbC(alpha, beta, Id, B);
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A * I + beta * C
        /// 
        /// with A : m-by-p matrix
        ///      I : p-by-n matrix; identity matrix
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * I</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix</param>
        public void mmCaAIbC(double alpha, double beta, Matrix A) 
        {
            if (A._rows != _rows || A._cols != _cols)
            {
                throw new MathException("Error in matrix multiplication. The dimension of the matrix A must match the current matrix.");
            }

            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    _data[i, j] = A._data[i, j] * alpha + _data[i, j] * beta;
                }
            }
        }

        /// <summary>
        /// Matrix multiplication of the form:
        /// 
        ///     C = alpha * A * I + beta * C
        /// 
        /// with A : m-by-p matrix
        ///      I : p-by-n matrix; identity matrix permuted according to a vector of permutations, piv
        ///      C : m-by-n matrix
        ///      alpha, beta : real numbers
        ///      
        /// </summary>
        /// <param name="alpha">The scalar value that multiplies A * I</param>
        /// <param name="beta">The scalar value that multiplies C</param>
        /// <param name="A">an object of type Matrix</param>
        /// <param name="piv">a vector of permutations on I, of type int</param>
        /// <param name="rows">The binary parameter to indicate whether the rows or the columns of I should be permuted (true for the rows; false for the columns)</param>
        public void mmCaAIbC(double alpha, double beta, Matrix A, int[] piv, bool rows) 
        { 
           	// TODO: implement a faster version without identity matrix
	        Matrix Id = new Matrix( A._cols, _cols, _dimT);
	        Id.set2Id();

	        if(rows)
	        {
		        // permute its columns
		        Id.cpermutem(piv);
	        }
	        else
	        {
		        // permute its rows
		        Id.rpermutem(piv);
	        }

	        mmCaABbC(alpha, beta, A, Id);
        }

        #endregion

        #region S O L V I N G   S Y S T E M S   O F   E Q U A T I O N S

        /// <summary>
        /// Permutes the columns of a matrix given a vector of permutations. <para/>
        /// <para/>
        /// For example, in case a matrix A is permuted after a QR decomposition with column pivoting,<para/>
        /// then the resulting matrix in the upper triangular matrix R. <para/>
        /// <para/>
        /// </summary>
        /// <param name="piv">a vector of permutations on the columns of A</param>
        /// <param name="trans">The boolean parameter to indicate whether to transpose the vector of permutations piv or not (=1, transpose; =0, otherwise). Default is false</param>
        public void cpermutem(int[] piv, bool trans = false)
        {
            for (int i = 0; i < _cols; i++)
            {
                if (piv[i] < 0 || piv[i] >= _cols)
                    throw new MathException("Invalid row permutation.");
            }

            Polynomial[,] newValues = initializePolynoms(_rows, _cols, _dimT);

            int[] pivT;

            if (trans)
            {
                pivT = new int[_cols];
                for (int i = 0; i < _cols; i++)
                {
                    pivT[piv[i]] = i;
                }
                piv = pivT;
            }

            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    newValues[i, j] = _data[i, piv[j]];
                }
            }

            _data = null;
            // GC.Collect();
            _data = newValues;
        }

        /// <summary>
        /// Permutes the rows of a matrix given a vector of permutations. 
        /// </summary>
        /// <param name="piv">a vector of permutations on the rows of A</param>
        public void rpermutem(int[] piv)
        {
            for (int i = 0; i < _rows; i++)
            {
                if (piv[i] < 0 || piv[i] >= _rows)
                {
                    throw new MathException("Invalid row permutation.");
                }
            }

            Polynomial[,] newValues = initializePolynoms(_rows, _cols, _dimT);
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    newValues[i, j] = _data[piv[i], j];
                }
            }

            _data = null;
            // GC.Collect();
            _data = newValues;
        }

        /// <summary>
        /// Transposes this matrix in place.
        /// </summary>
        public void transpose()
        {
            if (_rows == _cols)
            {
                // square matrix can be swapped in place
                for (int i = 1; i < _rows; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        Polynomial temp = _data[i, j];
                        _data[i, j] = _data[j, i];
                        _data[j, i] = temp;
                    }
                }
            }
            else
            {
                Polynomial[,] newValues = initializePolynoms(_cols, _rows, _dimT);
                for (int i = 0; i < _rows; i++)
                {
                    for (int j = 0; j < _cols; j++)
                    {
                        newValues[j, i] = _data[i, j];
                    }
                }
                _data = newValues;
                // number of rows and columns swapped
                int r = _rows;
                _rows = _cols;
                _cols = r;
            }
        }

        /// <summary>
        /// Creata the transposes of this matrix. This matrix opject remains unchanged. <para/>
        /// </summary>
        /// <returns>The transposed matrix object</returns>
        public Matrix asTranspose()
        {
            Matrix retval = new Matrix(this);
            retval.transpose();
            return retval;
        }

        /// <summary>
        /// Implements the shift operator to calculate the derivative of Taylor polynomials<para/>
        /// in case the elements of the matrix are such, like in:<para/>
        /// <para/>
        /// y(t) = sum_{j=0}^{d} y_j * t^j + O(t^d+1)<para/>
        /// 	 = y_0 + y_1*t + y_2*t^2 + ... + y_d*t^d<para/>
        /// <para/>
        /// y'(t) = y_1 + 2*y_2*t + 3*y_3*t^2 + ... + d*y_d*t^d-1<para/>
        /// <para/>
        /// Internally, the coefficients are shifted to the left and the last one is zeroed.<para/>
        /// </summary>
        public void shift()
        {
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    _data[i,j].shift();
                }
            }
        }

        /// <summary>
        /// Returns true in case the given matrix is the identity matrix; false otherwise.
        /// </summary>
        /// <returns>true if the matrix is the identity matrix; false otherwise.</returns>
        public bool isId()
        {
            // identity matrices must be also sqare matrices
            if (_rows != _cols)
            {
                return false;
            }

            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    if (i == j)
                    {
                        if (!_data[i,j].isId())
                        {
                            return false;
                        }
                    }
                    else
                    {
                        if (!_data[i,j].isZero())
                        {
                            return false;
                        }
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Returns true in case the given matrix is the zero matrix; false otherwise.
        /// </summary>
        /// <returns>true if the matrix is the zero matrix; false otherwise.</returns>
        public bool isZero()
        {
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    if (!_data[i, j].isZero())
                        return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Sets a matrix to the identity one:<para/>
        /// <para/>
        ///     M = I<para/>
        /// <para/>
        /// </summary>
        public void set2Id()
        {
            if (_rows != _cols)
            {
                throw new MathException("Only square matrices can be set to an identity matrix.");
            }
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                {
                    if (i == j)
                    {
                        _data[i,j].set2Id();
                    }
                    else
                    {
                        _data[i,j].set2Zero();
                    }
                }
            }
        }

        /// <summary>
        /// Sets a submatrix to the identity one:<para/>
        /// <para/>
        ///     e.g. M = (    | 1 0 0 |    )<para/>
        ///              ( M1 | 0 1 0 | M2 )<para/>
        ///              (    | 0 0 1 |    )<para/>
        ///              (        M3       )<para/>
        /// <para/>
        /// </summary>
        /// <param name="top">The number of rows at the top to keep unchanged</param>
        /// <param name="bottom">The number of rows at the bottom to keep unchanged</param>
        /// <param name="left">The number of columns on the left to keep unchanged</param>
        /// <param name="right">The number of columns on the right to keep unchanged</param>
        public void set2Id(int top, int bottom, int left, int right)
        {
            int lastRow = _rows - bottom - 1;
            int lastCol = _cols - right - 1;
            set2IdFromIndices(top, lastRow, left, lastCol);
        }

        /// <summary>
        /// Sets a submatrix to the identity one:<para/>
        /// <para/>
        ///     e.g. M = (    | 1 0 0 |    )<para/>
        ///              ( M1 | 0 1 0 | M2 )<para/>
        ///              (    | 0 0 1 |    )<para/>
        ///              (        M3       )<para/>
        /// </summary>
        /// <param name="firstRow">The row from which to start on</param>
        /// <param name="lastRow">The last row that should be considered</param>
        /// <param name="firstCol">The column from which to start on</param>
        /// <param name="lastCol">The last column that should be considered</param>
        public void set2IdFromIndices(int firstRow, int lastRow, int firstCol, int lastCol)
        {
            if (firstRow < 0 || lastRow >= _rows)
            {
                throw new MathException("Error in set2Id, row bounds are wrong.");
            }
            if (firstCol < 0 || lastCol >= _cols)
            {
                throw new MathException("Error in set2Id, column bounds are wrong");
            }

            int rows = lastRow - firstRow + 1;
            int cols = lastCol - firstCol + 1;

            if (rows != cols)
            {
                throw new MathException("Error in set2Id, submatrix must be a square matrix in order to become a identity matrix.");
            }

            for (int i = firstRow; i <= lastRow; i++)
            {
                for (int j = firstCol; j <= lastCol; j++)
                {
                    if (i - firstRow == j - firstCol)
                    {
                        _data[i,j].set2Id();
                    }
                    else
                    {
                        _data[i,j].set2Zero();
                    }
                }
            }

        }

        /// <summary>
        /// Sets a matrix to zero entries.
        /// </summary>
        public void set2Zero()
        {
	        set2Val(0.0);
        }

        /// <summary>
        /// Sets a submatrix to zero:<para/>
        /// <para/>
        ///     e.g. M = (    | 0 0 0 |    )<para/>
        ///              ( M1 | 0 0 0 | M2 )<para/>
        ///              (    | 0 0 0 |    )<para/>
        ///              (        M3       )<para/>
        /// </summary>
        /// <param name="top">The number of rows at the top to keep unchanged</param>
        /// <param name="bottom">The number of rows at the bottom to keep unchanged</param>
        /// <param name="left">The number of columns on the left to keep unchanged</param>
        /// <param name="right">The number of columns on the right to keep unchanged</param>
        public void set2Zero(int top, int bottom, int left, int right) 
        {
            int lastRow = _rows - bottom - 1;
            int lastCol = _cols - right - 1;
            set2ZeroFromIndices(top, lastRow, left, lastCol);
        }

        /// <summary>
        /// Sets a submatrix to zero:<para/>
        /// <para/>
        ///     e.g. M = (    | 0 0 0 |    )<para/>
        ///              ( M1 | 0 0 0 | M2 )<para/>
        ///              (    | 0 0 0 |    )<para/>
        ///              (        M3       )<para/>
        /// </summary>
        /// <param name="firstRow">The row from which to start on</param>
        /// <param name="lastRow">The last row that should be considered</param>
        /// <param name="firstCol">The column from which to start on</param>
        /// <param name="lastCol">The last column that should be considered</param>
        public void set2ZeroFromIndices(int firstRow, int lastRow, int firstCol, int lastCol) 
        {
            set2ValFromIndices(firstRow, lastRow, firstCol, lastCol, 0.0);
        }

        /// <summary>
        /// Sets a matrix to the value given as parameter.<para/>
        /// </summary>
        /// <param name="v">The double value to set the elements to</param>
        public void set2Val(double v)
        {
	        set2ValFromIndices(0, _rows - 1, 0, _cols - 1, v);
        }

        /// <summary>
        /// Sets a submatrix to the value given as parameter:<para/>
        /// <para/>
        ///     e.g. M = (    | v v v |    )<para/>
        ///              ( M1 | v v v | M2 )<para/>
        ///              (    | v v v |    )<para/>
        ///              (        M3       )<para/>
        /// </summary>
        /// <param name="top">The number of rows at the top to keep unchanged</param>
        /// <param name="bottom">The number of rows at the bottom to keep unchanged</param>
        /// <param name="left">The number of columns on the left to keep unchanged</param>
        /// <param name="right">The number of columns on the right to keep unchanged</param>
        /// <param name="v">The double value to set the elements to</param>
        public void set2Val(int top, int bottom, int left, int right, double v)
        {
            int lastRow = _rows - bottom - 1;
            int lastCol = _cols - right - 1;
            set2ValFromIndices(top, lastRow, left, lastCol, v);
        }

        /// <summary>
        /// Sets a submatrix to the value given as parameter:<para/>
        /// <para/>
        ///     e.g. M = (    | v v v |    )<para/>
        ///              ( M1 | v v v | M2 )<para/>
        ///              (    | v v v |    )<para/>
        ///              (        M3       )<para/>
        /// </summary>
        /// <param name="firstRow">The row from which to start on</param>
        /// <param name="lastRow">The last row that should be considered</param>
        /// <param name="firstCol">The column from which to start on</param>
        /// <param name="lastCol">The last column that should be considered</param>
        /// <param name="v">The double value to set the elements to</param>
        public void set2ValFromIndices(int firstRow, int lastRow, int firstCol, int lastCol, double v)
        {
            if (firstRow < 0 || lastRow >= _rows)
            {
                throw new MathException("Invalid row bounds.");
            }
            if (firstCol < 0 || lastCol >= _cols)
            {
                throw new MathException("Invalid column bounds.");
            }

            if (v == 0.0)
            {
                for (int i = firstRow; i <= lastRow; i++)
                {
                    for (int j = firstCol; j <= lastCol; j++)
                    {
                        _data[i,j].set2Zero();
                    }
                }
            }
            else
            {
                for (int i = firstRow; i <= lastRow; i++)
                {
                    for (int j = firstCol; j <= lastCol; j++)
                    {
                        _data[i,j].set2const(v);
                    }
                }
            }
        }

        #endregion

        /// <summary>
        /// Returns a String of the Matrix
        /// </summary>
        /// <returns>String of matrix</returns>
        override public String ToString()
        {
            String retval = String.Empty;
            for (int i = 0; i < _rows; i++)
            {
                for (int j = 0; j < _cols; j++)
                    retval += _data[i, j].ToString() + "\t";

                retval += System.Environment.NewLine;
            }

            return retval;
        }

        #endregion

        #endregion

        #region private functions

        private void initializeMatrix(int rows, int cols, int dimT, Polynomial[,] p)
        {
            _rows = rows;
            _cols = cols;
            _dimT = dimT;
            // creates a real copy, not only a pointer to the object
            _data = new Polynomial[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    _data[i,j] = new Polynomial(p[i,j]);
        }

        private Polynomial[,] initializePolynoms(int rows, int cols, int dimT)
        {
            Polynomial[,] p = new Polynomial[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    p[i, j] = new Polynomial(dimT);

            return p;
        }

        #endregion
    }
}
