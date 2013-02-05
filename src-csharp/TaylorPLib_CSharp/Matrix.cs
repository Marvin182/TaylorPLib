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
