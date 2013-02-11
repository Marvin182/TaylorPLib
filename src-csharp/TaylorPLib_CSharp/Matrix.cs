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
                throw new MathException(String.Format("(%d, %d) is not a alid index of %d x %d matrix", row, col, _rows, _cols));

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
                    throw new MathException(String.Format("(%d, %d) is not a alid index of %d x %d matrix", row, col, _rows, _cols));

                return _data[row, col];
            }
            set
            {
                if (row >= _rows || row < 0 || col >= _cols || col < 0)
                    throw new MathException(String.Format("(%d, %d) is not a alid index of %d x %d matrix", row, col, _rows, _cols));

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
		        throw new MathException(String.Format("Cannot add up a %d x %d matrix (left side) and a %d x %d (right side) matrix.", a._rows, a._cols, b._rows, b._cols));
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
                throw new MathException(String.Format("Cannot add up a %d x %d matrix (left side) and a %d x %d (right side) matrix.", a._rows, a._cols, b._rows, b._cols));
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
                throw new MathException(String.Format("Cannot multiply a %d x %d with a %d x %d matrix. The number of rows in m must match the number of colmuns in this matrix.", a._rows, a._cols, b._rows, b._cols));

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
                throw new MathException(String.Format("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols));
            }
            if (A._rows != _rows || B._cols != _cols)
            {
                throw new MathException(String.Format("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols));
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
                throw new MathException(String.Format("The parameter r (= %d) must be smaller or equal than the number of rows in B (%d x %d).", r, B._rows, B._cols));
            }
            if (c > B._cols)
            {
                throw new MathException(String.Format("The parameter c (= %d) must be smaller or equal than the number of columns in B (%d x %d).", c, B._rows, B._cols));
            }
            if (A._cols != B._rows)
            {
                throw new MathException(String.Format("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols));
            }
            if (A._rows != _rows || B._cols != _cols)
            {
                throw new MathException(String.Format("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols));
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
		        throw new MathException(String.Format("The number of rows in A that are of interest (r = %d) must be smaller or equal than the total amount of rows in A (%d x %d).", r, A._rows, A._cols));
	        }
	        if (A._cols != B._rows)
	        {
		        throw new MathException(String.Format("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols));
	        }
	        if (A._rows != _rows || B._cols != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols));
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
		        throw new MathException(String.Format("The number of columns in B that are of interest (r = %d) must be smaller or equal than the total amount of columns in B (%d x %d).", r, A._rows, A._cols));
	        }
	        if (A._cols != B._rows)	
	        {
		        throw new MathException(String.Format("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols));
	        }
	        if (A._rows != _rows || B._cols != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols));
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
		        throw new MathException(String.Format("Cannot multiply A (%d x %d) und B (%d x %d).", A._rows, A._cols, B._rows, B._cols));
	        }
	        if (A._rows != _rows || B._cols != _cols)
	        {
		        throw new MathException(String.Format("The size of the matrix A * B (%d x %d) must match the current matrix (%d x %d).", A._rows, B._cols, _rows, _cols));
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
