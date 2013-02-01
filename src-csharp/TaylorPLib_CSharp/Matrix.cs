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
        
        private int _rows,
                    _cols,
                    _dimT;
        private Polynomial[,] _data;

        #endregion

        #region public getter and setter 

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
            initializeMatrix(rows, cols, dimT, initializePolynoms(rows, cols));
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
            initializeMatrix(rows, cols, 0, values);
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

        private Polynomial[,] initializePolynoms(int rows, int cols)
        {
            Polynomial[,] p = new Polynomial[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    p[i, j] = new Polynomial();

            return p;
        }

        #endregion
    }
}
