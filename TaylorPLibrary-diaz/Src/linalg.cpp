/*!
 * \file linalg.cpp
 * \brief Linear algebra, matrix computations.
 * \author D. Monett Diaz
 * \date June, 2006
 * 
 * This class implements some linear algebra algorithms and related functions.
 * 
 * Included files:
 * 	\file linalg.h header file.
 * 
 * \todo
 * - release all temporary used memory
 * 
 */

#include "linalg.h"

using namespace std;

/**
 * Returns the index of the maximum value in a vector, from position a to b.
 * 
 * \param[in] v The pointer to \a v, a vector of type \type T.
 * \param[in] selected The pointer to \a selected, a vector of type \type int.
 * \param[in] P The pointer to \a P, a vector of permutations on \a v.
 * \param[in] a The starting position.
 * \param[in] b The ending position.
 * \return The index of the maximum value that was found, as an type \type int.
 */
int max( double *v, bool *selected, int *P, int a, int b )
{
	int i, j = -1;											// j = index for the max. value
	double m = -1;											// only usefull in case all v[i] >= 0
	/*printf( " a in max = %d\n", a );
	printf( " b in max = %d\n", b );*/
	if( a == b )
		j = a;
	for( i = a; i < b; i++ )
		if( !selected[ P[ i ] ] )
			if( v[ P[ i ] ] > m )							// values in v are not swapped but in P!!
			{
				j = i;
				m = v[ P[ i ] ];
			}
	selected[ P[ j ] ] = true;
	/*printf( " found j in max = %d\n", j );
	printf( " P[ j ] in max = %d\n", P[ j ] );
	printf( " selected[ P[ j ] ] in max = %d\n", selected[ P[ j ] ] );
	printf( "selected (in max)= \n" );
	for( i = 0; i < b; i++ )
		printf( "  %d%c", selected[ i ], '\t' );
	printf( "\n\n" );*/
	return j;
}

/**
 * Returns the index i such that v[ i ] == t.
 * 
 * It is assumed that such a value exists.
 * 
 * \param[in] v The pointer to \a v, a vector of type \type int.
 * \param[in] t The type \type int value to compare with.
 * \return The index that is found, as an \type int.
 * 
 */
int findIndex( int *v, int t )
{
	int i = -1;
	do {
		i++;
	}
	while( v[ i ] != t ); 
	return i;
}

/**
 * Dot product of two vectors from components r to m-1.
 * 
 *     c = x^T y
 * 
 * \param[in] r The starting component.
 * \param[in] m The ending component.
 * \param[in] x The pointer to \a x, a vector of type \type T.
 * \param[in] y The pointer to \a y, a vector of type \type T.
 * \return The result of the calculation, as a type \type T.
 * 
 */
template <typename T> T dpc( int r, int m, T *x, T *y )
{
	T c( x[ 0 ].order() );
	
	c.set2zero();
	for( int i = r; i < m; i++ )
		c += x[ i ] * y[ i ];
	//c.print();
	return c;
}

/**
 * Transposes a given vector of permutations. 
 * 
 * \param[in] n The dimension of the vector.
 * \param[in] x The pointer to \a x, a vector of \type int.
 * \param[in] y The pointer to \a y, the resulting transposed vector, of \type int.
 * \return The error code.
 * 
 */
int transpv( int n, int *x, int *y )
{
	try
	{
		for( int i = 0; i < n; i++ )
			y[ x[ i ] ] = i;
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Returns false in case at least one element from a given array is greater that a threshold.
 * 
 * \param[in] n The dimension of the array \a v.
 * \param[in] v The pointer to \a v, with elements of \type double.
 * \param[in] eps The threshold value to compare with.
 * \param[out] elem The first element (its index) of the array \a v not meeting the condition.
 * \return \a true if all elements are lesser or equal that the threshold. 
 * 		Otherwise it returns \a false.
 * 
 */
bool vleqeps( int n, double *v, double eps, int & elem )
{
	bool leq = true;
	
	elem = -1;												// some initialization
	for( int i = 0; i < n; i++ )
		if( v[ i ] > eps )
		{
			leq = false;
			elem = i;
			break;
		}

	return leq;
}

/**
 * Constructs the matrix R from a given matrix A.
 * 
 * The upper triangular part of A contains R, after A is pivoted. Then, in case an upper diagonal 
 * matrix R is desired, then the call should look like this:
 * 
 * 			constructR( A, piv, R, 1 )
 * 
 * If the resulting R should have non-permuted columns, then call
 * 
 * 			constructR( A, piv, R, 0 )
 * 
 * \param[in] A The pointer to \a A, a matrix of \type T.
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
 * \param[out] R The pointer to \a R, the resulting matrix.
 * \param[in] permuted The binary parameter to indicate whether R should be 
 * 		permuted or not (= 1, yes).
 * \return The error code.
 * 
 */
template <typename T> int constructR( Matrix<T> &A, int *piv, Matrix<T> &R, bool permuted )
{
	int i, j, k;
	
	try
	{
		if( permuted )										// for a permuted R, permuted=1
		{
			for( i = 0; i < A.nrows(); i++ )
				for( j = 0; j < A.ncols(); j++ )
				{
					if( j >= i )
						R( i, j ) = A( i, piv[ j ] );
					else
						R( i, j ).set2zero();				// set to the null polynomial
				}
		} 
		else 
		{													// for a non-permuted R, permuted=0
			for( j = 0; j < A.ncols(); j++ )
			{
				k = findIndex( piv, j );
				for( i = 0; i < k+1; i++ )
					R( i, j ) = A( i, j );
				for( i = k+1; i < A.nrows(); i++ )
					R( i, j ).set2zero();					// set to the null polynomial
			}
		}
	}
	catch(...)
	{
		IDException e( "Error when constructing the matrix R from a given matrix A.", 12 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Applies a Householder reflection to a matrix from the left hand side.
 * 
 * See Section 5.4.1 "Applying Householder Matrices", p.211 from Golub & Van Loan's book.
 * 
 * 		P*A = (beta*v*v^T)*A 
 * 			= A - v*w^T
 * 
 * with  w = beta*A^T*v
 * 
 * where:
 * 		A : m-by-n matrix
 * 		P : m-by-m matrix
 * 		v : Householder m-vector
 * 		beta: calculated by Householder
 * 	
 * Observation: matrix-vector multiplication and outer product updates are CHEAPER 
 * 				than matrix-matrix multiplication!!!
 * 
 * \param[in] pos The starting position.
 * \param[in] beta The value of \a beta after calculating the Householder vector.
 * \param[in/out] A The pointer to \a A, a matrix of \type T. On output, its lower triangular
 * 			part contains the Householder vectors.
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
 * \param[out] v The pointer to \a v, a vector of \type T that is the Householder vector.
 * \return The error code.
 * 
 */
template <typename T> int hhlh( int pos, T beta, Matrix<T> &A, int *piv, T *v )
{
	int i, j;
	T *w;
	
	try
	{
		w = new T[ A.ncols() ];
		if( (w == 0)  ||  (w == NULL) )
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < A.ncols(); i++ )
			w[ i ] = T( beta.order() );						// grade of Taylor polynomials
															// (same as of beta)
		for( j = pos; j < A.ncols(); j++ )					// calculate w = beta A^T v  with...
		{													// ...A(r:m,r:n) and v m-r vector
			w[ j ] = A( pos, piv[ j ] );					// initialize w_j, when v[r]=1 !!
			for( i = pos + 1; i < A.nrows(); i++ )
				w[ j ] += A( i, piv[ j ] ) * v[ i ];
			w[ j ] *= beta;
		}
		for( i = pos; i < A.nrows(); i++ )					// calculate A = A - v w^T
			for( j = pos; j < A.ncols(); j++ )
				A( i, piv[ j ] ) -= v[ i ] * w[ j ];
		for( i = pos + 1; i < A.nrows(); i++ )
			A( i, piv[ pos ] ) = v[ i ];
	
		delete [] w;
	}
	catch( IDException e )
	{
		delete [] w;
        e.report();
		throw e.getErrCode();
	}
	catch(...)
	{
		delete [] w;
		IDException e( "Error when applying a Householder reflection to a matrix.", 16 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Applies a Householder reflection to a matrix from the left hand side.
 * 
 * See Section 5.4.1 "Applying Householder Matrices", p.211 from Golub & Van Loan's book.
 * 
 * 		P*A = (beta*v*v^T)*A 
 * 			= A - v*w^T
 * 
 * with  w = beta*A^T*v
 * 
 * where:
 * 		A : m-by-n matrix
 * 		P : m-by-m matrix
 * 		v : Householder m-vector
 * 		beta: calculated by Householder
 * 	
 * Observation: matrix-vector multiplication and outer product updates are CHEAPER 
 * 				than matrix-matrix multiplication!!!
 * 
 * \param[in] pos The starting position.
 * \param[in] beta The value of \a beta after calculating the Householder vector.
 * \param[in/out] A The pointer to \a A, a matrix of \type T. On output, its lower triangular
 * 			part contains zeros!!
 * \param[in] piv The pointer to \a piv, a vector of permutations on the columns of \a A.
 * \param[in] v The pointer to \a v, a vector of \type T, containing the Householder vector.
 * \return The error code.
 * 
 */
template <typename T> int hhlr( int pos, T beta, Matrix<T> &A, int *piv, T *v )
{
	int i, j;
	T *w;
	
	try
	{
		w = new T[ A.ncols() ];
		if( (w == 0)  ||  (w == NULL) )
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < A.ncols(); i++ )
			w[ i ] = T( beta.order() );						// grade of Taylor polynomials
															// (same as of beta)
		for( j = pos; j < A.ncols(); j++ )					// calculate w = beta A^T v  with...
		{													// ...A(r:m,r:n) and v m-r vector
			w[ j ] = A( pos, piv[ j ] );					// initializew_j, when v[r]=1 !!
			for( i = pos + 1; i < A.nrows(); i++ )
				w[ j ] += A( i, piv[ j ] ) * v[ i ];
			w[ j ] *= beta;
		}
		for( i = pos; i < A.nrows(); i++ )					// calculate A = A - v w^T
			for( j = pos; j < A.ncols(); j++ )
				A( i, piv[ j ] ) -= v[ i ] * w[ j ];
		for( i = pos + 1; i < A.nrows(); i++ )				// make lower triangular A zero!!
															// => R is constructed
			A( i, piv[ pos ] ).set2zero();					// set to p(x) = 0
	
		delete [] w;
	}
	catch( IDException e )
	{
		delete [] w;
        e.report();
		throw e.getErrCode();
	}
	catch(...)
	{
		delete [] w;
		IDException e( "Error when applying a Householder reflection to a matrix.", 16 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Applies a Householder reflection to a matrix from the right hand side (see /211 Golub's book):
 * 
 * See Section 5.4.1 "Applying Householder Matrices", p.211 from Golub & Van Loan's book.
 * 
 * 		A*P = A*(beta*v*v^T)
 * 			= A - w*v^T
 * 
 * with  w = beta*A*v
 * 
 * where:
 * 		A : m-by-n matrix
 * 		P : m-by-m matrix
 * 		v : Householder m-vector
 * 		beta: calculated by Householder
 * 		
 * Observation: matrix-vector multiplication and outer product updates are CHEAPER 
 * 				than matrix-matrix multiplication!!!
 * 
 * \param[in] rowpos The starting row.
 * \param[in] colpos The starting column.
 * \param[in] beta The value of \a beta after calculating the Householder vector.
 * \param[in/out] A The pointer to \a A, a matrix of \type T.
 * \param[in] v The pointer to \a v, a vector of \type T, containing the Householder vector.
 * \return The error code.
 * 
 */
template <typename T> int hhr( int rowpos, int colpos, T beta, Matrix<T> &A, T *v )
{
	int i, j;
	T *w;
	
	try
	{
		w = new T[ A.nrows() ];
		if( (w == 0)  ||  (w == NULL) )
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < A.nrows(); i++ )
			w[ i ] = T( beta.order() );						// grade of Taylor polynomials
															// (same as of beta)
		/*printf( "beta =" );
		beta.print();
		printf( "\nv =" );
		for( i = 0; i < A.nrows(); i++ )
		{
			printf( "\n  [ %d ] : ", i );
			v[ i ].print();
		}
		printf( "\n" );*/
		for( i = rowpos; i < A.nrows(); i++ )				// calculate, e.g., w = beta Q v
		{
			w[ i ] = A( i, colpos );						// initialize w_i, when v[r]=1 !!
			for( j = colpos + 1; j < A.nrows(); j++ )
				w[ i ] += A( i, j ) * v[ j ] ;
			w[ i ] *= beta;
		}
		for( j = colpos; j < A.nrows(); j++ )				// calculate, e.g., Q = Q - w v^T
			for( i = rowpos; i < A.nrows(); i++ )
				A( i, j ) -= w[ i ] * v[ j ];
	
		delete [] w;
	}
	catch( IDException e )
	{
		delete [] w;
        e.report();
		throw e.getErrCode();
	}
	catch(...)
	{
		delete [] w;
		IDException e( "Error when applying a Householder reflection to a matrix.", 16 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the Householder matrix:
 * 
 *     H = 2 (v v^T) / (v^T v)
 * 
 * given a Householder vector v, and verifies some of its properties.
 * 
 * Observation: Only components r to m of v are interesting.
 * 
 * \param[in] r The starting component.
 * \param[in] eps The threshold value used in stop criterion.
 * \param[in] A The pointer to \a A, a matrix of \type T.
 * \param[in] P The pointer to \a P, a vector of permutations on the columns of \a A.
 * \param[in] v The pointer to \a v, a vector of \type T, containing the Householder vector.
 * \return The error code.
 * 
 */
template <typename T> int hm( int r, double eps, Matrix<T> &A, int *P, T *v, T beta )
{
	try
	{
		//
		// When TESTING hm(...) !!
		//
		// For computing the Householder matrix given, e.g.,  v = [1.4773,0.4444,1.0]^T
		// in the case of v being a vector of doubles
		//
		// v1 = 2*v*v^T/(v^T*v)
		// v1 =
		//		1.2914 0.3885 0.8742
		//		0.3885 0.1169 0.2630
		//		0.8742 0.2630 0.5917
		//
		// H =
		//		-0.2914 -0.3885 -0.8742
		//		-0.3885 0.8831 -0.2630
		//		-0.8742 -0.2630 0.4083
		//
		// Then do:
		//
		/*T *vv = new T[ 3 ];
		vv[ 0 ] = 1.4773;
		vv[ 1 ] = 0.4444;
		vv[ 2 ] = 1.0;
		hm( 0, eps, vv, P, A ); 							// for this example
		hm( r, eps, v, P, A );								// for the general case
		// In case H should be declared in the calling program and passed as parameter, do:
		hm( 0, eps, vv, P, A, H ); 							// for this example
		hm( r, eps, v, P, A, H );							// for the general case
		// and eliminate the internal declaration of H in hm
		//
		// See examples like D:\AD\doc\matrix\householderQR-espannol.pdf and others
		*/
	
		int i, j, k, m = A.nrows();
		T d( A.dimT() ), cTP( A.dimT() ), den( A.dimT() );	// aux. Taylor polynomials
		Matrix<T> H( m, m, A.dimT() ), I( m, m, A.dimT() ),	// aux.matrices
			C1( m, m, A.dimT() ), C2( m, m, A.dimT() );
			
		printf( "\nComputing the Householder matrix and testing...\n" );
		I.set2Id();
		printf( "v =");
		for( i = r; i < m; i++ )
		{
			printf( "\n  [ %d ] : ", i );
			v[ i ].print();
		}
		printf( "\n\n" );
		//den = dpc( r, m, v, v );							// the dot product v^T*v is calculated
		//printf( "v^T*v = " );
		//den.print();
		//printf( "\n" );
		//if( den.feval() == 0.0 )
		//	throw IDException( "Division by zero.", 20 );
		//cTP.set2const( 2.0 );
		//d = cTP / den;
		for( i = r; i < m; i++ )
			for( j = r; j < m; j++ )
			{
				H( i, j ) = v[ j ] * v[ i ];				// this is an outer product v*v^T
				H( i, j ) *= beta;							// equivalent to H( i, j ) *= d;
			}
		H.printtpm( r, m, r, m, "beta * (v v^T) = \n", red, eps );
		//H.printtpm( r, m, r, m, "2 * (v v^T) / (v^T v) = \n", red, eps );
		cTP.set2const( 1.0 );
		for( i = r; i < m; i++ )
			for( j = r; j < m; j++ )
			{
				if( i == j ) {
					H( i, j ) = cTP - H( i, j );			// diagonal from the Id matrix is used
				} 
				else 
				{
					H( i, j ) = - H( i, j );				// rest of the elements of I are zero
				}
			}
		H.printtpm( r, m, r, m, "H = I - beta*(v v^T) = \n", red, eps );
		//H.printtpm( r, m, r, m, "H = I - 2*(v v^T)/(v^T v) = \n", red, eps );
		//
		// Verify: H * H = I   !!
		// The complete matrix multiplication is done, but only for testing!!!
		//
		C1.set2zero( r, m, r, m );
		C2.set2zero( r, m, r, m );
		for( j = r; j < m; j++ )
			for( k = r; k < m; k++ )
				for( i = r; i < m; i++ )
				{
					C1( i, j ) += H( i, k ) * H( j, k );
					C2( i, j ) += H( i, k ) * H( k, j );
				}
		C1.printtpm( r, m, r, m, "H * H^T = (should be I) \n", red, eps );
		if( C1.isId( r, m, r, m, eps ) )					// H * H^T == I ?
			printf( "H * H^T = I\n\n" );
		else
		{
			printf( red );					
			printf( "-->  H * H^T <> I (for eps=%.8lg)\n\n", eps );
			printf( normal );
		}
		C2.printtpm( r, m, r, m, "H * H = (should also be I) \n", red, eps );
		if( C2.isId( r, m, r, m, eps ) )					// H * H == I ?
			printf( "H * H = I\n\n" );
		else
		{
			printf( red );					
			printf( "-->  H * H <> I (for eps=%.8lg)\n\n", eps );
			printf( normal );
		}
		//
		// Verify: H * A zeroes the P[r] column from A(r+1,P[r]) on   !!
		//
		printf( "H * A = ( P[r] column from A zeroes from A[r+1, P[r]] on) \n" );
		for( i = r; i < m; i++ )
		{
			d.set2zero();
			printf( "\nd = " );
			d.print();
			for( j = r; j < m; j++ )
			{
				printf( "\n\nH( %d, %d ) = ", i, j );
				H( i, j ).print();
				printf( "\n\nA( %d, %d ) = ", j, P[ r ] );
				A( j, P[ r ] ).print();
				d += H( i, j ) * A( j, P[ r ] );
				printf( "\n\nd += H( %d, %d ) * A( %d, %d ) = ", i, j, j, P[ r ] );
				d.print();
			}
			printf( green );
			printf( "\n\n   ==> updated A( %d, %d ) = ", i, P[ r ] );
			d.print();
			printf( normal );
			printf( "\n" );
		}
		printf( "\n" );
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)
	{
		IDException e( "Error when computing a Householder reflection.", 17 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes a Householder vector.
 * 
 * This function implements the Algorithm 5.1.1 p. 210, from Golub & Van Loan's
 * "Matrix<T> Computations", 3rd Edition, Johns Hopkins University Press, 1996.
 * 
 * Observation: The interesting part of v is the one from r to m (and not as in 
 * the book where length of v = m - r because A(r:m,r) is to be considered).
 * 
 * \param[in] r The starting position.
 * \param[in] A The pointer to \a A, a matrix of \type T.
 * \param[in] P The pointer to \a P, a vector of permutations on the columns of \a A.
 * \param[out] v The pointer to \a v, a vector of \type T that is the Householder vector.
 * \param[in] eps The threshold value used in stop criterion.
 * \param[out] beta The value of \a beta after calculating the Householder vector.
 * \return The error code.
 * 
 */
template <typename T> int hhv( int r, Matrix<T> &A, int *P, T *v, double eps, T &beta )
{
	try
	{
		int i, 
			ord = beta.order();								// grade of Taylor polynomials
		T sigma( ord ), aux( ord ), miu( ord ), cTP( ord ), tempTP( ord );
		
		/*printf( "\nComputing the Householder vector...\n" );
		A.printtpm( "Matrix A (at start) = \n", yellow, eps );
		printf( "\nr = %d\n", r );
		printf( "\nP[ r ] = %d\n", P[ r ] );
		printf( "\nElement A[0][0] =\n  ");
		A( 0, 0 ).print();
		printf( "\ncolumn %d from A =", P[ r ] );
		for( i = 0; i < A.nrows(); i++ )
		{
			printf( "\n  [ %d ] : ", i );
			A( i, P[ r ] ).print();
		}
		printf( "\n" );*/
		sigma.set2zero();
		//printf( "\nsigma (init) =\n  ");
		//sigma.print();
		for( i = r + 1; i < A.nrows(); i++ )				// calculate A(r+1:m,r) norm
		{
			//printf( "\n  %d, sigma += A( i, P[ r ] ) * A( i, P[ r ] ): ", i );
			sigma += A( i, P[ r ] ) * A( i, P[ r ] );		// rows r+1 to m at col. P[ r ]
			//sigma.print();
		}
		//printf( green );
		//printf( "\nv (initialization) = ");
		v[ r ].set2const( 1.0 );							// v(1) = 1
		//printf( "\n  [ %d ] : ", r );
		//v[ r ].print();
		for( i = r + 1; i < A.nrows(); i++ )				// v(2:length) = x(2:length)
		{
			v[ i ] = A( i, P[ r ] );
			//printf( "\n  [ %d ] : ", i );
			//v[ i ].print();
		}
		//printf( normal );
		//printf( "\n" );
		aux = A( r, P[ r ] );								// first element
		if( fabs( sigma.feval() ) <= eps )					// feval() = first Taylor coefficient!
			beta.set2zero();
		else 
		{
			tempTP = aux * aux;
			tempTP += sigma;
			miu = tempTP.sqrt();							// miu = sqrt( x(1)^2 + sigma )
			//printf( "\n\nmiu =\n  ");
			//miu.print();
			cTP.set2zero();
			if( aux <= cTP )								// if x(1) <= 0
			{ 												// implicit: if aux.feval() <= cTP.feval()
				v[ r ] = aux - miu;
				/*printf( "\naux" );
				aux.print();
				printf( "\nmiu" );
				miu.print();
				printf( "\nA( %d, %d ) <= 0  then  ", r, P[ r ] );
				printf( green );
				printf( "v[ %d ] = ", r );
				v[ r ].print();*/
			} 
			else 
			{
				/*printf( "\naux before v[ r ] = - sigma / ( aux + miu ) =\n  " );
				aux.print();
				printf( "\nmiu before v[ r ] = - sigma / ( aux + miu ) =\n  " );
				miu.print();
				printf( "\nsigma before v[ r ] = - sigma / ( aux + miu ) =\n  " );
				sigma.print();*/
				v[ r ] = - sigma / ( aux + miu );
				/*printf( "\naux after =\n  " );
				aux.print();
				printf( "\nmiu after =\n  " );
				miu.print();
				printf( "\nsigma after =\n  " );
				sigma.print();
				printf( "\nA( %d, %d ) > 0  then v[ %d ] = ", r, P[ r ], r );
				printf( "- sigma / ( A( %d, %d ) + miu ) =\n  ", r, P[ r ] );
				printf( green );
				v[ r ].print();*/
			}
			//printf( normal );
			tempTP = v[ r ];								// overrides old value of tempTP
			tempTP.sqr();
			//printf( "\n\ntempTP = sqr of v[ %d ] =\n  ", r );
			//tempTP.print();
			tempTP = v[ r ] * v[ r ];						// overrides old value of tempTP
			//printf( "\n\nv[ %d ]*v[ %d ] =\n  ", r, r );
			//tempTP.print();
			tempTP = sigma / tempTP;
			//printf( "\n\nsigma/(v[ %d ]*v[ %d ]) =\n  ", r, r );
			//tempTP.print();
			cTP.set2const( 1.0 );
			tempTP += cTP;
			//printf( "\n\n[sigma/(v[ %d ]*v[ %d ])]+1 =\n  ", r, r );
			//tempTP.print();
			beta.set2const( 2.0 );
			beta /= tempTP;									// same as 2*v(1)^2 / (sigma+v(1)^2)
			/*printf( magenta );
			printf( "\n\nbeta = 2/(sigma/(v[ %d ]*v[ %d ])+1) =\n  ", r, r );
			beta.print();
			//old version: beta = 2 / ( 1 + sigma / ( v[ r ] * v[ r ] ) );;
			//beta.print();
			printf( green );
			printf( "\n\nv /= v[ 0 ]");*/
			aux = v[ r ];
			for( i = r; i < A.nrows(); i++ )				// update resting values
			{
				v[ i ] /= aux;
				//printf( "\n  [ %d ] : ", i );
				//v[ i ].print();
			}
			//printf( normal );
			//printf( "\n" );
		}
	}
	catch(...)
	{
		IDException e( "Error when computing a Householder vector.", 18 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Householder QR Factorization with Column Pivoting.
 * 
 * This function implements the Algorithm 5.4.1 p. 249, from Golub & Van Loan's
 * "Matrix<T> Computations", 3rd Edition, Johns Hopkins University Press, 1996.
 * 
 * In the case of rank deficiency, the QR factorization of a matrix A can be computed 
 * as follows:
 * 
 *     A P = Q R
 * 
 * where P is a column permutation.
 * 
 * A is a m-by-n matrix, with m >= n (m rows, n columns). Its rank r is calculated by the algorithm.
 * Matrix<T> Q = H_1 ... H_r is a product of Householder matrices H_i, with 1<=i<=r.
 * Matrix<T> P = P_1 ... P_r is a product of interchange matrices P_j, with 1<=j<=r.
 * 
 * In particular:
 * 
 *     Q = Q_k ... Q_2 Q_1 
 * 
 * where 
 * 
 *     Q_i = tau_i v_i v_i^T 
 * 
 * and v_i is the Householder vector 
 * 
 *     v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i))
 * 
 * and tau_i is related to the permutations that are done.
 * 
 * \param[in] eps The threshold value used in stop criterion.
 * \param[in/out] A The pointer to \a A, a matrix of \type T. On input it contains the matrix 
 * 			to be factorized. On output, the matrix R (diagonal and upper triangular part).
 * \param[out] P The pointer to \a P, a vector of permutations on the columns of \a A. 
 * \param[out] Q The pointer to \a Q, a matrix of \type T, containing the 
 * 		resulting \a Q matrix.
 * \param[out] rank The calculated rank.
 * \param[in] permuted The binary parameter to indicate whether A should be permuted or 
 * 		not on output (= 1, yes).
 * \param[in] printWhat The \type int parameter for printing out.
 * \param[in] printWhere The \type int parameter for printing out to the screen, to a file, or nowhere.
 * \param[in] fn The name of the file to print out to.
 * \return The error code.
 * 
 */
template <typename T> int hqrfcp( double eps, Matrix<T> &A, int *P, Matrix<T> &Q, 
	int& rank, bool permuted, int printWhat, int printWhere, FILE * fn )
{	
	int i, j, k, r, t, e, 
		ord = A.dimT();										// order of Taylor polynomials
	T *v;
	T beta( ord );
	bool *selected, cond = false;
	double *c, pivNotCond, pivCond;
	
	try
	{
		//printf( "\nEntering hqrfcp...\n" );
		
		//////
		//  S E T U P 
		//////
		
		c = new double[ A.ncols() ];						// column norms
		selected = new bool[ A.ncols() ];					// whether c[j] considered or not
		v = new T[ A.nrows() ];								// for the Householder vectors
		if( (c == 0)  ||  (c == NULL)  ||  (selected == 0)  ||  (selected == NULL)  ||
			(v == 0)  ||  (v == NULL) )
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < A.nrows(); i++ )
			v[ i ] = T( ord );								// grade of Taylor polynomials
		Matrix<T> OA( A.nrows(), A.ncols(), A.dimT() );		// copy of A, for debugging purposes
		Matrix<T> C( A.nrows(), A.nrows(), A.dimT() );		// used only for debugging purposes
		Matrix<T> C2( A.nrows(), A.ncols(), A.dimT() );		// used only for debugging purposes
		
		Q.set2Id();											// Q = I_m at start
		for( i = 0; i < A.nrows(); i++ )
			for( j = 0; j < A.ncols(); j++ )
				OA( i, j ) = A( i, j );						// original matrix, for testing purposes
		for( j = 0; j < A.ncols(); j++ )
		{
			selected[ j ] = false;							// at start: no element is selected
			P[ j ] = j ;									// initialize P
		}
		//A.printm( "A (at start) = \n", normal );
	
		A.colnorm( c );										// column norms
		//printdv( A.ncols(), c, "c (column norms at start) = \n", normal );
		//printiv( A.ncols(), selected, "selected (at start)= ", normal );
	
		//////
		//  F I R S T   P I V O T I N G 
		//////
		
		k = max( c, selected, P, 0, A.ncols() );			// max in c is searched for
		r = -1;												// initialize index in P
/**///		printf( "\nk = %d", k );
/**///		printf( "\nP[ k ] = %d", P[ k ] );

		//////
		//  P R I N T I N G   P I V O T   E L E M E N T S
		//////
		
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintf( fn, "\nFirst Householder pivot element = %.16lg", c[ P[ k ] ] );
			else if( printWhere == 0 )
			{
				printf( green );
				printf( "\nFirst Householder pivot element = %.16lg", c[ P[ k ] ] );
				printf( normal );
			}
		}
		if( c[ P[ k ] ] > eps )								// initializations
		{
			cond = true;
			pivCond = c[ P[ k ] ];
		}
		else pivNotCond = c[ P[ k ] ];

		//////
		//  M A I N   L O O P  hqrfcp
		//////
		
		while( c[ P[ k ] ] > eps )							// epsilon is a parameter, e.g., 1E-13
		{
			pivCond = c[ P[ k ] ];							// update pivot holding condition
			r++;
			e = P[ r ];										// swap columns r and k in P
			P[ r ] = P[ k ];
			P[ k ] = e;
			//printiv( A.ncols(), P, "P = \n", cyan );
			printf( normal );
			
			//////
			//  R E D U C T I O N  
			//////
			
			hhv( r, A, P, v, eps, beta );					// compute Householder vector v
			
			//////
			//  M A T R I X   U P D A T E
			//////
			
			hhlr( r, beta, A, P, v );						// (beta v v^T) A
			//A.printm( "A (after update, with zeros in upper triangular part, 
						//non-permuted) = \n", yellow );
			//A.printm( P, "A (after update, with zeros in upper triangular part, 
						//w.col.piv.) = \n", yellow );
	
			//////
			//  N O R M   D O W N D A T E 
			//////
			
			A.colnormdown( r, P, c );
			//printdv( A.ncols(), c, P, "c (downdate) = \n", normal );
	
			//////
			//  N E X T   P I V O T I N G
			//////
		
			if( r < A.ncols() - 1 )							// last column of A is not considered
			{
				k = max( c, selected, P, r + 1, A.ncols() );// max in c is searched for
				if( c[ P[ k ] ] > eps )						// new pivot element
				{
					pivCond = c[ P[ k ] ];
					cond = true;
				}
				else
				{
					pivNotCond = c[ P[ k ] ];
					cond = false;
				}
			}
			else 
			{
				if( c[ P[ k ] ] > eps )						// last one before setting to zero
				{
					pivCond = c[ P[ k ] ];
					cond = true;
				}
				else
				{
					pivNotCond = c[ P[ k ] ];
					cond = false;
				}
				c[ P[ k ] ] = 0.0;
			}
	
			//////
			//  Q   C O M P U T A T I O N 
			//////
	
			//
			// Forward accumulation for computing Q is used 
			// (see p. 213 from Golub's book): Q = Q*Q_j.
			// 
			// Then:
			//
			//		Q+ = Q (beta v v^T) = Q - w v^T 
			//
			// with  w = beta Q v
			//
			// The Householder vector v remains the same as for A, as well as beta does.
			// Components from 0 to r-1 of v should be zero; they are anyway unneeded 
			// for the calculations.
			// 
			
			hhr( 0, r, beta, Q, v );						// Q (beta v v^T)
	
			//////
			//  T E S T S   I N S I D E   T H E   L O O P  :  Q, Q^T*origA, Q^T*Q, Q*Q^T
			//////
	
			//Q.printm( "Q (hqrfcp) = \n", magenta );		// Q
	
			C2.mmCaATBPbC( 1.0, 1.0, Q, OA, P );			// Q^T * OA (with col. piv.)
			//C2.printm( "Q^T * A (w.col.piv.) = \n", yellow );
			
			C.mmCaATAbC( 1.0, 1.0, Q );						// Q^T * Q
			//C.printm( "Q^T * Q = (should be I)\n", yellow );
			
			C.mmCaAATbC( 1.0, 1.0, Q );						// Q * Q^T
			//C.printm( "Q * Q^T = (should be also I)\n", yellow );
			
		}													// END OF big while
	
		//////
		//  P R I N T I N G   P I V O T   E L E M E N T S
		//////
		
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
			{
				if( cond )
					fprintf( fn, "\nLast Householder pivot = %.16lg > eps (it holds the condition for eps=%.16lg)", pivCond, eps );
				fprintf( fn, "\nLast Householder pivot = %.16lg <= eps (it does not hold the condition for eps=%.16lg)", pivNotCond, eps );
			}
			else if( printWhere == 0 )
			{
				printf( green );
				if( cond )
					printf( "\nLast Householder pivot = %.16lg > eps (it holds the condition for eps=%.16lg)", pivCond, eps );
				printf( "\nLast Householder pivot = %.16lg <= eps (it does not hold the condition for eps=%.16lg)", pivNotCond, eps );
				printf( normal );
			}
		}
		
		//////
		//  F I N A L   T E S T S  :  Q*R, Q*R-A*P, A, R
		//////
	
		C2.mmCaAUTBPbC( 1.0, 1.0, Q, A, P );				// Q * R, upper triangular A contains R
		/*if( printWhat > 21 )
		{
			if( printWhere == 1 )
			{
				C2.fprinttpm( fn, "Q*R = \n", eps );
				OA.fprinttpm( fn, P, "(original)A*P = \n", eps );// A * P, with the original A!!!
			}
			else if( printWhere == 0 )
			{
				C2.printtpm( "Q*R = \n", yellow, eps );
				OA.printtpm( P, "(original)A*P = \n", yellow, eps );// A * P, with the original A!!!
			}
		}*/
		
		for( i = 0; i < A.nrows(); i++ )
			for( j = 0; j < A.ncols(); j++ )
				C2( i, j ) -= OA( i, P[ j ] );				// Q*R - A*P
		/*if( printWhat > 21 )
		{
			if( printWhere == 1 )
				C2.fprinttpm( fn, "Q*R - A*P = (should be zeroed)\n", eps );
			else if( printWhere == 0 )
				C2.printtpm( "Q*R - A*P = (should be zeroed)\n", yellow, eps );
		}*/
		
		if( permuted )										// permuted = 1, 
		{													// A really upper triangular on output!
			A.cpermutem( P );
			/*if( printWhat > 21 )
			{
				if( printWhere == 1 )
					A.fprinttpm( fn, "A (on output = R, A permuted) = \n", eps );
				else if( printWhere == 0 )
					A.printtpm( "A (on output = R, A permuted) = \n", yellow, eps );
			}*/
		} 
		else 
		{													// permuted=0, col. order not considered
			/*if( printWhat > 21 )
			{
				if( printWhere == 1 )
					A.fprinttpm( fn, "A (on output = R, A non-permuted) = \n", eps );
				else if( printWhere == 0 )
					A.printtpm( "A (on output = R, A non-permuted) = \n", yellow, eps );
			}*/
		}
	
		rank = r + 1;
		//printf( "\nRank = %d\n", rank );
		printf( normal );
	
		//////
		//  C L E A N I N G   O P E R A T I O N S  (last-in-first-out)
		//////
		delete [] v;
		delete [] selected;
		delete [] c;		
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)
	{
		IDException e( "Error in the QR factorization with column pivoting.", 19 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Householder QR Factorization with Column Pivoting.
 * 
 * This function implements the Algorithm 5.4.1 p. 249, from Golub & Van Loan's
 * "Matrix<T> Computations", 3rd Edition, Johns Hopkins University Press, 1996.
 * 
 * In the case of rank deficiency, the QR factorization of a matrix A can be computed 
 * as follows:
 * 
 *     A P = Q R
 * 
 * where P is a column permutation.
 * 
 * A is a m-by-n matrix, with m >= n (m rows, n columns). 
 * Its rank r is calculated by the algorithm.
 * Matrix<T> Q = H_1 ... H_r is a product of Householder matrices H_i, with 1<=i<=r.
 * Matrix<T> P = P_1 ... P_r is a product of interchange matrices P_j, with 1<=j<=r.
 * 
 * In particular:
 * 
 *     Q = Q_k ... Q_2 Q_1 
 * 
 * where 
 * 
 *     Q_i = tau_i v_i v_i^T 
 * 
 * and v_i is the Householder vector 
 * 
 *     v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i))
 * 
 * and tau_i is related to the permutations that are done.
 * 
 * \param[in] eps The threshold value used in stop criterion.
 * \param[in/out] A The pointer to \a A, a matrix of \type T. On input it contains the matrix 
 * 			to be factorized. On output, its diagonal and upper triangular part contains 
 * 			the matrix R; its lower triangular part contains the Householder vectors (which 
 * 			encode the orthogonal matrix Q).
 * \param[out] P The pointer to \a P, a vector of permutations on the columns of \a A. 
 * \param[out] Q The pointer to \a Q, a matrix of \type T, containing the resulting 
 * 			\a Q matrix.
 * \param[out] R The pointer to \a R, a matrix of \type T, containing the resulting 
 * 			\a R matrix.
 * \param[out] rank The calculated rank.
 * \param[in] permuted The binary parameter to indicate whether R should be permuted or not 
 * 			on output (= 1, yes).
 * \param[in] printWhat The \type int parameter for printing out.
 * \param[in] printWhere The \type int parameter for printing out to the screen, to a file, or nowhere.
 * \param[in] fn The name of the file to print out to.
 * \return The error code.
 * 
 */
template <typename T> int hqrfcp( double eps, Matrix<T> &A, int *P, Matrix<T> &Q, 
		Matrix<T> &R, int& rank, bool permuted, int printWhat, int printWhere, FILE * fn )
{	
	int i, j, k, r, t, e, 
		ord = A.dimT();										// order of Taylor polynomials
	T *v;													// Householder vectors
	T beta( ord );
	bool *selected, cond = false;
	double *c, pivNotCond, pivCond;

	try
	{
		//printf( "\nEntering hqrfcp...\n" );
		
		//////
		//  S E T U P 
		//////
		
		c = new double[ A.ncols() ];						// column norms
		selected = new bool[ A.ncols() ];					// whether c[j] considered or not
		v = new T[ A.nrows() ];
		if( (c == 0)  ||  (c == NULL)  ||  (selected == 0)  ||  (selected == NULL)  ||
			(v == 0)  ||  (v == NULL) )
			throw IDException( "Memory allocation failure.", 3 );
		for( i = 0; i < A.nrows(); i++ )
			v[ i ] = T( ord );								// with the grade of Taylor polynomials
		Matrix<T> OA( A.nrows(), A.ncols(), A.dimT() );		// copy of original A
		Matrix<T> C( A.nrows(), A.nrows(), A.dimT() );		// used only for debugging purposes
		Matrix<T> C2( A.nrows(), A.ncols(), A.dimT() );		// used only for debugging purposes

		Q.set2Id();											// Q = I_m at start
		for( i = 0; i < A.nrows(); i++ )
			for( j = 0; j < A.ncols(); j++ )
				OA( i, j ) = A( i, j );						// original matrix, for testing purposes
		for( j = 0; j < A.ncols(); j++ )
		{
			selected[ j ] = false;							// at start: no element is selected
			P[ j ] = j ;									// initialize P
		}
/**///		A.printtpm( "Matrix (at start) = \n", yellow, eps );
		A.colnorm( c );										// column norms
/**///		printdv( A.ncols(), c, "c (column norms at start) = \n", normal );
	
		//////
		//  F I R S T   P I V O T I N G 
		//////
		
		k = max( c, selected, P, 0, A.ncols() );			// max in c is searched for
		r = -1;												// initialize index in P
		
		//////
		//  P R I N T I N G   P I V O T   E L E M E N T S
		//////
		
/**///		printf( "\nk = %d", k );
/**///		printf( "\nP[ k ] = %d", P[ k ] );
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
				fprintf( fn, "\nFirst Householder pivot element = %.16lg", c[ P[ k ] ] );
			else if( printWhere == 0 )
			{
				printf( green );
				printf( "\nFirst Householder pivot element = %.16lg", c[ P[ k ] ] );
				printf( normal );
			}
		}
		if( c[ P[ k ] ] > eps )								// initializations
		{
			cond = true;
			pivCond = c[ P[ k ] ];
		}
		else pivNotCond = c[ P[ k ] ];

		//////
		//  M A I N   L O O P  hqrfcp
		//////
		
		while( c[ P[ k ] ] > eps )							// epsilon is a parameter, e.g., 1E-13
		{
			pivCond = c[ P[ k ] ];							// update pivot holding condition
			r++;
			e = P[ r ];										// swap columns r and k in P
			P[ r ] = P[ k ];
			P[ k ] = e;
/**///		printiv( A.ncols(), P, "P = \n", cyan );
			
			//////
			//  R E D U C T I O N  
			//////
			
			hhv( r, A, P, v, eps, beta );					// compute Householder vector v
/**///			for( i = 0; i < A.nrows(); i++ )
/**///			{
/**///				printf( "Householder v[ %d ] = ", i );
/**///				v[ i ].print();
/**///				printf( "\n" );
/**///			}
/**///			printf( "beta = " );
/**///			beta.print();
/**///			printf( "\n" );
			
			//////
			// T E S T
			//////
/**///			hm( r, eps, A, P, v, beta );				// compute Householder matrix and test
			
			//////
			//  M A T R I X   U P D A T E
			//////
			
			hhlh( r, beta, A, P, v );						// (beta v v^T) A
/**///			A.printtpm( "A (after update, with Householder vectors) = \n", green, eps );
	
			//////
			//  N O R M   D O W N D A T E 
			//////
			
			A.colnormdown( r, P, c );
/**///			printdv( A.ncols(), c, P, "c (downdate) = \n", normal );
	
			//////
			//  N E X T   P I V O T I N G
			//////
		
			if( r < A.ncols() - 1 )							// last column of A is not considered
			{
				k = max( c, selected, P, r + 1, A.ncols() );// max in c is searched for
				if( c[ P[ k ] ] > eps )						// new pivot element
				{
					pivCond = c[ P[ k ] ];
					cond = true;
				}
				else
				{
					pivNotCond = c[ P[ k ] ];
					cond = false;
				}
			}
			else 
			{
				if( c[ P[ k ] ] > eps )						// last one before setting to zero
				{
					pivCond = c[ P[ k ] ];
					cond = true;
				}
				else
				{
					pivNotCond = c[ P[ k ] ];
					cond = false;
				}
				c[ P[ k ] ] = 0.0;
			}
/**///			printf( magenta );
/**///			printf( "\nNext pivot Householder= %.8lg\n", c[ P[ k ] ] );
/**///			printf( normal );
				
			//////
			//  Q   C O M P U T A T I O N 
			//////
	
			//
			// Forward accumulation for computing Q is used
			// (see p. 213 from Golub's book): Q = Q*Q_j.
			// 
			// Then:
			//
			//		Q+ = Q (beta v v^T) = Q - w v^T 
			//
			// with  w = beta Q v
			//
			// The Householder vector v remains the same as for A, as well as beta does.
			// Components from 0 to r-1 of v should be zero; they are anyway unneeded 
			// for the calculations.
			// 
			
/**///		Q.printtpm( "Q (before forward accumulation) = \n", magenta, eps );	// Q
/**///		printf( "\nr = %d\n", r );
			hhr( 0, r, beta, Q, v );						// Q (beta v v^T)
/**///		Q.printtpm( "Q (after forward accumulation) = \n", magenta, eps );	// Q
	
			//////
			//  T E S T S   I N S I D E   T H E   L O O P  :  Q, R, P, Q^T*Q, Q*Q^T
			//////
	
/**///			C2.mmCaATBPbC( 1.0, 1.0, Q, OA, P );			// Q^T * OA (with col. piv.)
/**///			C2.printtpm( "Q^T * A (w.col.piv., should be R) = \n", yellow, eps );
			
/**///			C.mmCaATAbC( 1.0, 1.0, Q );						// Q^T * Q
/**///			C.printtpm( "Q^T * Q = (should be I)\n", yellow, eps );
			
/**///			C.mmCaAATbC( 1.0, 1.0, Q );						// Q * Q^T
/**///			C.printtpm( "Q * Q^T = (should be also I)\n", yellow, eps );
			
		}													// END OF big while
		
		//////
		//  P R I N T I N G   P I V O T   E L E M E N T S
		//////
		
		if( printWhat > 1 )
		{
			if( printWhere == 1 )
			{
				if( cond )
					fprintf( fn, "\nLast Householder pivot = %.16lg > eps (it holds the condition for eps=%.16lg)", pivCond, eps );
				fprintf( fn, "\nLast Householder pivot = %.16lg <= eps (it does not hold the condition for eps=%.16lg)", pivNotCond, eps );
			}
			else if( printWhere == 0 )
			{
				printf( green );
				if( cond )
					printf( "\nLast Householder pivot > eps (it holds the condition for eps=%.16lg) = %.16lg", eps, pivCond );
				else printf( "\nLast Householder pivot <= eps (it does not hold the condition for eps=%.16lg) = %.16lg", eps, pivNotCond );
				printf( normal );
			}
		}
		
		//////
		//  F I N A L   T E S T S  :  Q*R, Q*R-A*P
		//////
	
		if( permuted )										// permuted R, upper triangular on output
		{
			constructR( A, P, R, 1 );
			/*if( printWhat > 21 )
			{
				if( printWhere == 1 )
					R.fprinttpm( fn, "R (w.col.piv.) = \n", eps );
				else if( printWhere == 0 )
					R.printtpm( "R (w.col.piv.) = \n", yellow, eps );
			}*/
			C2.mmCaABbC( 1.0, 1.0, Q, R );					// Q * R
			/*if( printWhat > 21 )
			{
				if( printWhere == 1 )
					C2.fprinttpm( fn, "Q*R (R w.col.piv.) = \n", eps );
				else if( printWhere == 0 )
					C2.printtpm( "Q*R (R w.col.piv.) = \n", yellow, eps );
			}*/
		} 
		else 
		{													// non-permuted R, altered columns
			//A.printm( "A (before constructR) = \n", yellow );
			constructR( A, P, R, 0 );
/**///			R.printtpm( "R = \n", cyan, eps );
/**///			Q.printtpm( "Q = \n", magenta, eps );
			/*if( printWhat > 21 )
			{
				if( printWhere == 1 )
				{
					R.fprinttpm( fn, "R = \n", eps );
					Q.fprinttpm( fn, "Q = \n", eps );
				}
				else if( printWhere == 0 )
				{
					R.printtpm( "R = \n", yellow, eps );
					Q.printtpm( "Q = \n", yellow, eps );
				}
			}*/
			C2.mmCaAUTBPbC( 1.0, 1.0, Q, A, P );			// Q * R, upper triangular A after 
															// permuted contains R
			/*if( printWhat > 21 )
			{
				if( printWhere == 1 )
					C2.fprinttpm( fn, "Q*R = \n", eps );
				else if( printWhere == 0 )
					C2.printtpm( "Q*R = \n", yellow, eps );
			}*/
/**///			C2.printtpm( "Q*R = \n", blue, eps );
		}
		/*if( printWhat > 21 )
			if( printWhere == 1 )
				OA.fprinttpm( fn, P, "(original)A*P = \n", eps );// A * P, with the original A!!!
			else if( printWhere == 0 )
				OA.printtpm( P, "(original)A*P = \n", yellow, eps );// A * P, with the original A!!!
		*/
/**///			OA.printtpm( P, "(original)A*P = \n", red, eps );
		for( i = 0; i < A.nrows(); i++ )
			for( j = 0; j < A.ncols(); j++ )
				C2( i, j ) -= OA( i, P[ j ] );				// Q*R - A*P
		/*if( printWhat > 21 )
			if( printWhere == 1 )
				C2.fprinttpm( fn, "Q*R - A*P = (should be zeroed)\n", eps );
			else if( printWhere == 0 )
				C2.printtpm( "Q*R - A*P = (should be zeroed)\n", yellow, eps );
		*/
/**///			C2.printtpm( "Q*R - A*P = (should be zeroed)\n", normal, eps );
		rank = r + 1;
		//printf( "\nRank = %d\n", rank );
		printf( normal );
	
		//////
		//  C L E A N I N G   O P E R A T I O N S   (last-in-first-out)
		//////
		delete [] v;
		delete [] selected;
		delete [] c;		
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)
	{
		IDException e( "Error in the QR factorization with column pivoting.", 19 );
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Calculates P, W, and Q for the matrix sequence, where
 * 
 * 		P = G^- * G
 * 		W = G * G^-
 * 		Q = P
 * 
 * \param[in] m The number of rows/columns in \a G.
 * \param[out] P The pointer to \a P, a matrix of \type T with the calculated \a P. 
 * \param[out] W The pointer to \a W, a matrix of \type T with the calculated \a W. 
 * \param[out] Q The pointer to \a Q, a matrix of \type T with the calculated \a Q. 
 * \return The error code.
 * 
 */
template <typename T> int pwq( Matrix<T> &G, Matrix<T> &Gm, Matrix<T> &P, 
		Matrix<T> &W, Matrix<T> &Q )
{
	int m = G.nrows();
	Matrix<T> I( m, m, G.dimT() );

	I.set2Id();												// set to the identity matrix
	P.mmCaABbC( 1.0, 1.0, Gm, G );
	W.mmCaABbC( 1.0, 1.0, G, Gm );
	W = I - W;
	Q = I - P;
	
	return 0;
}

/**
 * Forms the matrix SL:
 * 
 * 		SL = ( SA1  SA2 )
 * 			 (  0    I  )
 * 
 * \param[in] r The row number that delimits the size of \a SA1.
 * \param[in] SA1 The pointer to \a SA1, a submatrix of \type T. 
 * \param[in] SA2 The pointer to \a SA2, a submatrix of \type T. 
 * \param[out] SL The pointer to \a SL, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int formSL( int r, Matrix<T> &SA1, Matrix<T> &SA2, Matrix<T> &SL )
{
	int i, j, n = SL.ncols();								// SL[n][n]
	
	if( r == 0 )											// SA1 does not exist but SA2
	{														// then
		for( i = 0; i < r; i++ )							// SL = ( SA2 )
			for( j = 0; j < n; j++ )						//      (  I  )
				SL( i, j ) = SA2( i, j );
		SL.set2Id( r, n, 0, n );
	}
	else if( n - r == 0 )									// SA2 does not exist but SA1
	{														// then
		for( i = 0; i < r; i++ )							// SL = ( SA1 )
			for( j = 0; j < n; j++ )						//      (  0  )
				SL( i, j ) = SA1( i, j );
		SL.set2zero( r, n, 0, n );
	}
	else													// both SA1 and SA2 exist
	{														// then
		for( i = 0; i < r; i++ )							// SL = ( SA1  SA2 )
			for( j = 0; j < r; j++ )						//      (  0    I  )
				SL( i, j ) = SA1( i, j );
		for( i = 0; i < r; i++ )
			for( j = r; j < n; j++ )
				SL( i, j - r ) = SA2( i, j );
		SL.set2zero( r, n, 0, r );
		SL.set2Id( n - r, n, n - r, n );
	}
	return 0;
}

/**
 * Forms the matrix SR:
 * 
 * 		SR = ( SD1^T  0 )
 * 			 ( SD2^T  I )
 * 
 * \param[in] r The row number that delimits the size of \a SD1.
 * \param[in] SD1 The pointer to \a SD1, a submatrix of \type T. 
 * \param[in] SD2 The pointer to \a SD2, a submatrix of \type T. 
 * \param[out] SR The pointer to \a SR, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int formSR( int r, Matrix<T> &SD1, Matrix<T> &SD2, Matrix<T> &SR )
{
	int i, j, n = SR.ncols();								// SR[n][n]
	
	try
	{	
		Matrix<T> SD1T( SD1.ncols(), SD1.nrows(), SD1.dimT() ),	// SD1^T
			SD2T( SD2.ncols(), SD2.nrows(), SD2.dimT() ); 	// SD2^T
			
		if( r == 0 )										// SD1 does not exist but SD2
		{													// then
			SD2T = SD2.transpm();							// SR = ( SD2^T  I )
			for( i = 0; i < r; i++ )
				for( j = 0; j < n; j++ )
					SR( i, j ) = SD2T( i, j );
			SR.set2Id( 0, r, r, n );
		}
		else if( n - r == 0 )								// SD2 does not exist but SD1
		{													// then
			SD1T = SD1.transpm();							// SR = ( SD1^T  0 )
			for( i = 0; i < r; i++ )
				for( j = 0; j < n; j++ )
					SR( i, j ) = SD1T( i, j );
			SR.set2zero( 0, r, r, n );
		}
		else												// both SD1 and SD2 exist
		{													// then
			SD1T = SD1.transpm();							// SR = ( SD1^T  0 )
			SD2T = SD2.transpm();							//      ( SD2^T  I )
			for( i = 0; i < r; i++ )
				for( j = 0; j < r; j++ )
					SR( i, j ) = SD1T( i, j );
			for( i = r; i < n; i++ )
				for( j = 0; j < r; j++ )
					SR( i, j ) = SD2T( i - r, j );
			SR.set2zero( 0, r, r, n );
			SR.set2Id( n - r, n, n - r, n );
		}
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Decomposes a matrix into submatrices, where the submatrices M11 and M12 are of interest:
 * 
 * 		M = [ M11  M12 ]
 *          [ ...  ... ]
 * 
 * \param[in] r The row number that delimits the size of \a M11.
 * \param[in] M The pointer to \a M, a matrix of \type T to be decomposed into two submatrices. 
 * \param[out] M11 The pointer to \a M11, the first submatrix of \a M of \type T. 
 * \param[out] M12 The pointer to \a M12, the second submatrix of \a M of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int decompM( int r, Matrix<T> &M, Matrix<T> &M11, Matrix<T> &M12 )
{
	int i, j;
	
	for( i = 0; i < r; i++ )
	{
		for( j = 0; j < r; j++ )
			M11( i, j ) = M( i, j );
		for( j = r; j < M.ncols(); j++ )
			M11( i, j - r ) = M( i, j );
	}
	return 0;
}

/**
 * Decomposes a matrix into two submatrices R1 and R2:
 * 
 * 		    ( R1  R2 )
 * 		R = (      0 )
 * 
 * The columns of R may be permuted. The vector of permutations is given and the column
 * number where to make the cut in R, too.
 * 
 * Special cases:
 * 
 *     r = 0  =>  R1 does not exist but R2
 * 0 < r < n  =>  both R1 and R2 exist
 * 0 < r = n  =>  R2 does not exist but R1
 * 
 * \param[in] r The column number that delimits the size of \a R1.
 * \param[in] R The pointer to \a R, a matrix of \type T to be decomposed. 
 * \param[in] P The pointer to \a P, the vector of column permutations in \a R. 
 * \param[out] R1 The pointer to \a R1, a sub matrix of \a R of \type T. 
 * \param[out] R2 The pointer to \a R2, a sub matrix of \a R of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int decompS( int r, Matrix<T> &R, int *P, Matrix<T> &R1, Matrix<T> &R2 )
{
	int i, j;
	
	if( r == 0 )											// R1 does not exist but R2
	{
		for( i = 0; i < r; i++ )
			for( j = 0; j < R.ncols(); j++ )
				R2( i, j ) = R( i, P[ j ] );
	}
	else if( R.ncols() - r == 0 )							// R2 does not exist but R1
	{
		for( i = 0; i < r; i++ )
			for( j = 0; j < R.ncols(); j++ )
				R1( i, j ) = R( i, P[ j ] );
	}
	else													// both R1 and R2 exist
	{
		for( i = 0; i < r; i++ )
			for( j = 0; j < r; j++ )
				R1( i, j ) = R( i, P[ j ] );
		for( i = 0; i < r; i++ )
			for( j = r; j < R.ncols(); j++ )
				R2( i, j - r ) = R( i, P[ j ] );
	}
	return 0;
}

/**
 * Decomposes a matrix into two submatrices V1 and V2:
 * 
 * 		V = [ V1 | V2 ]
 * 
 * Special cases:
 * 
 *     r = 0  =>  V1 does not exist but V2
 * 0 < r < n  =>  both V1 and V2 exist
 * 0 < r = n  =>  V2 does not exist but V1
 * 
 * \param[in] r The column number that delimit the size of \a V1.
 * \param[in] V The pointer to \a V, a matrix of \type T to be decomposed into two submatrices. 
 * \param[out] V1 The pointer to \a V1, the first submatrix of \a V of \type T. 
 * \param[out] V2 The pointer to \a V2, the second submatrix of \a V of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int decompV( int r, Matrix<T> &V, Matrix<T> &V1, Matrix<T> &V2 )
{
	int i, j;
	
	if( r == 0 )											// V1 does not exist but V2
		V2 = V;
	else if( V.ncols() - r == 0 )							// V2 does not exist but V1
		V1 = V;
	else													// both V1 and V2 exist
	{
		for( i = 0; i < V.nrows(); i++ )
		{
			for( j = 0; j < r; j++ )
				V1( i, j ) = V( i, j );
			for( j = r; j < V.ncols(); j++ )
				V2( i, j - r ) = V( i, j );
		}
	}
	return 0;
}

/**
 * Decomposes a matrix into two submatrices W1 and W2:
 * 
 * 		W = [ W1 ]
 *          [ W2 ]
 * 
 * Special cases:
 * 
 *     r = 0  =>  W1 does not exist but W2
 * 0 < r < n  =>  both W1 and W2 exist
 * 0 < r = n  =>  W2 does not exist but W1
 * 
 * \param[in] r The row number that delimit the size of \a W1.
 * \param[in] W The pointer to \a W, a matrix of \type T to be decomposed into two submatrices. 
 * \param[out] W1 The pointer to \a W1, the first submatrix of \a W of \type T. 
 * \param[out] W2 The pointer to \a W2, the second submatrix of \a W of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int decompW( int r, Matrix<T> &W, Matrix<T> &W1, Matrix<T> &W2 )
{
	int i, j;
	
	if( r == 0 )											// W1 does not exist but W2
		W2 = W;
	else if( W.nrows() - r == 0 )							// W2 does not exist but W1
		W1 = W;
	else													// both W1 and W2 exist
	{
		for( j = 0; j < W.ncols(); j++ )
		{
			for( i = 0; i < r; i++ )
				W1( i, j ) = W( i, j );
			for( i = r; i < W.nrows(); i++ )
				W2( i - r, j ) = W( i, j );
		}
	}
	return 0;
}

/**
 * Forms the matrix BM:
 * 
 * 		BM (or "big matrix") =  ( BM11  BM12 )
 * 								( BM21  BM22 )
 * 
 * \param[in] d The row number that delimits the size of \a BM11.
 * \param[in] BM11 The pointer to \a BM11, a submatrix of \type T. 
 * \param[in] BM12 The pointer to \a BM12, a submatrix of \type T. 
 * \param[in] BM21 The pointer to \a BM21, a submatrix of \type T. 
 * \param[in] BM22 The pointer to \a BM22, a submatrix of \type T. 
 * \param[out] BM The pointer to \a BM, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int formBM( int d, Matrix<T> &BM11, Matrix<T> &BM12, Matrix<T> &BM21, 
		Matrix<T> &BM22, Matrix<T> &BM )
{
	int i, j;
	
	for( i = 0; i < d; i++ )
	{
		for( j = 0; j < d; j++ )							// BM11
			BM( i, j ) = BM11( i, j );
		for( j = d; j < BM.ncols(); j++ )					// BM12
			BM( i, j ) = BM12( i, j - d );
	}
	for( i = d; i < BM.nrows(); i++ )
	{
		for( j = 0; j < d; j++ )							// BM21
			BM( i, j ) = BM21( i - d, j );
		for( j = d; j < BM.ncols(); j++ )					// BM22
			BM( i, j ) = BM22( i - d, j - d );
	}
	return 0;
}

/**
 * Forms the matrix BM:
 * 
 * 		BM (or "big matrix") = ( BM11  BM12 )  OR  ( BM11 )
 * 												   ( BM21 )
 * 
 * depending on a boolean variable that indicates the orientation.
 * 
 * \param[in] d The row number that delimits the size of \a BM11.
 * \param[in] vert The boolean variable indicating whether to form \a BM vertically or
 * 			horizontally (=1, vertically; horizintally otherwise).
 * \param[in] BM11 The pointer to \a BM11, a submatrix of \type T. 
 * \param[in] BM12 The pointer to \a BM12, a submatrix of \type T. 
 * \param[out] BM The pointer to \a BM, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int formBM( int d, bool vert, Matrix<T> &BM1, Matrix<T> &BM2, Matrix<T> &BM )
{
	int i, j;
	
	if( vert )												// BM11 and BM21 needed
	{
		for( j = 0; j < BM.ncols(); j++ )
		{
			for( i = 0; i < d; i++ )						// BM11
				BM( i, j ) = BM1( i, j );
			for( i = d; i < BM.nrows(); i++ )				// BM21
				BM( i, j ) = BM2( i - d, j );
		}
	}
	else													// BM11 and BM12 needed
	{
		for( i = 0; i < BM.nrows(); i++ )
		{
			for( j = 0; j < d; j++ )						// BM11
				BM( i, j ) = BM1( i, j );
			for( j = d; j < BM.ncols(); j++ )				// BM12
				BM( i, j ) = BM2( i, j - d );
		}
	}
	return 0;
}

/**
 * Updates the matrix S as follows:
 * 
 * 		S = ( S     )
 *			(    S1 )
 * 
 * with:
 * 		S : m-by-m matrix (before update)
 * 			and (m+r)-by-(m+r) matrix (after update)
 * 		S1 : r-by-r matrix
 * 
 * The right dimensions of updated S should be guaranteed in the calling program.
 * 
 * \param[in] r The number of rows/columns in \a S1.
 * \param[in] S1 The pointer to \a S1, a submatrix of \type T. 
 * \param[out] S The pointer to \a S, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int filloutS( int r, Matrix<T> &S1, Matrix<T> &S )
{
	int i, j, m = S.nrows(), n = S.ncols();
	Matrix<T> aux( m, n, S.dimT() );						// auxiliary matrix
	
	if( r == 0 )											// only the first time
		S = S1;
	else
	{
		aux = S;											// save old S
		S.redim( m + r, n + r, S.dimT() );					// redim. S
		S.set2zero( 0, m, n, n + r );						// 0
		S.set2zero( m, m + r, 0, n );						// 0
		for( i = 0; i < m; i++ )
			for( j = 0; j < n; j++ )
				S( i, j ) = aux( i, j );					// copy old S
		for( i = m; i < m + r; i++ )
			for( j = n; j < n + r; j++ )
				S( i, j ) = S1( i - m, j - n );				// S1
	}
	//S.printtpm( "\nS = \n", yellow );
	return 0;
}

/**
 * Fills out a matrix M as follows:
 * 
 * 		M = ( M , V )
 * 
 * with:
 * 		M : m-by-d1 matrix (before update)
 * 			m-by-d1+d2 matrix (after update)
 * 		V : m-by-d2 matrix
 * 
 * \param[in] d1 The number of the column from \a M from which to add the new columns on.
 * \param[in] d2 The number of columns that are added to, from \a V.
 * \param[in] V The pointer to \a V, a matrix of \type T. 
 * \param[out] M The pointer to \a M, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int filloutM( int d1, int d2, Matrix<T> &V, Matrix<T> &M )
{
	int i, j;
	
	if( d2 != 0 )											// there are columns to add to
	{
		for( i = 0; i < M.nrows(); i++ )
			for( j = d1; j < d1 + d2; j++ )					// from column d1 on
				M( i, j ) = V( i, j - d1 );
		//M.printm( "M = ( M , V ) = \n", yellow );
	}
	return 0;
}

/**
 * Computes the matrix Tlim:
 * 
 * 		Tlim =  (   I     	) * ( I          ) = (   I        0    )
 * 				( mi1*Si  I )   (    Vitilde )	 ( mi1*Si  Vitilde )
 * 
 * with:
 * 		mi1 : (m-r)-by-r matrix
 * 		Si : r-by-r matrix
 * 		Tlim : m-by-m matrix
 * 		Vitilde : (m-r)-by-(m-r) matrix; it might be permuted
 * 
 * \param[in] mi1 The pointer to \a mi1, a matrix of \type T. 
 * \param[in] Si The pointer to \a Si, a matrix of \type T. 
 * \param[in] Vitilde The pointer to \a Vitilde, a matrix of \type T. 
 * \param[out] Tlim The pointer to \a Tli, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computeTlim( Matrix<T> &mi1, Matrix<T> &Si, 
		Matrix<T> &Vitilde, Matrix<T> &Tlim )
{
	int i, j, r = Si.nrows(), m = Tlim.nrows();
	
	try
	{	
		Matrix<T> aux( m - r, r, Si.dimT() );				// aux[m-r][r]
		
		//mi1.printm( "\nmi1", yellow );
		//Si.printm( "\nSi", yellow );
		//Vitilde.printm( "\nVitilde", yellow );
		aux.mmCaABbC( 1.0, 1.0, mi1, Si );					// mi1*Si
		//aux.printm( aux, "mi1*Si = \n", yellow );
		Tlim.set2Id( r, r );								// I
		Tlim.set2zero( 0, r, r, m );						// 0
		for( i = r; i < m; i++ )
		{
			for( j = 0; j < r; j++ )						// mi1*Si
				Tlim( i, j ) = aux( i - r, j );
			for( j = r; j < m; j++ )						// Vitilde
				Tlim( i, j ) = Vitilde( i - r, j - r );
		}
		//Tlim.printm( "Tlim = \n", yellow );
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the matrix mTui as follows:
 * 
 * 		mTui = ( I            ) * ( I  Si*mi2  ) = ( I   Si*mi2   )
 * 			   (    Uitilde^T )   (       I    )   ( 0  Uitilde^T )
 * 
 * with:
 * 		mi2 : r-by-(m-r) matrix
 * 		Si : r-by-r matrix
 * 		Uitilde : (m-r)-by-(m-r) matrix
 * 		mTui : m-by-m matrix
 * 
 * \param[in] mi2 The pointer to \a mi2, a matrix of \type T. 
 * \param[in] Si The pointer to \a Si, a matrix of \type T. 
 * \param[in] Uitilde The pointer to \a Uitilde, a matrix of \type T. 
 * \param[out] mTui The pointer to \a mTui, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computemTui( Matrix<T> &mi2, Matrix<T> &Si, 
		Matrix<T> &Uitilde, Matrix<T> &mTui )
{
	int i, j, r = Si.nrows(), m = mTui.nrows();
	
	try
	{
		Matrix<T> aux1( r, m - r, Si.dimT() ),				// aux1[r][m-r]
			aux2( Uitilde.ncols(), Uitilde.nrows(), Uitilde.dimT() );// Uitilde^T
		aux1.mmCaABbC( 1.0, 1.0, Si, mi2 );					// Si*mi2
		aux2 = Uitilde.transpm();							// transpose Uitilde
		mTui.set2Id( r, r );								// I
		mTui.set2zero( r, m, 0, r );						// 0
		for( i = 0; i < r; i++ )
			for( j = r; j < m; j++ )						// Si*mi2
				mTui( i, j ) = aux1( i, j - r );
		for( i = r; i < m; i++ )
			for( j = r; j < m; j++ )						// Uitilde^T
				mTui( i, j ) = aux2( i - r, j - r );
		//mTui.printm( "mTui = \n", yellow );
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the matrix Vi (V0 in the initialization of the matrix sequence and
 * Vtildem in the while-cycle).
 * 
 * No row or column is permuted.
 * 
 * \param[in] M The pointer to \a M, a matrix of \type T that is used in the computation. 
 * \param[out] Vi The pointer to \a Vi, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computeV( Matrix<T> &M, Matrix<T> &Vi )
{
	int i, j, d1 = Vi.nrows(), d2 = Vi.ncols(), m1 = M.nrows(), m2 = M.ncols();
	
	try
	{	
		Vi.set2Id( m1, d2 - m2 );							// first d2-n columns: fill out with I
		for( i = 0; i < m1; i++ )
			for( j = d2 - m2 + 1; j < d2; j++ )				// resting columns: fill out with M
				Vi( i, j ) = M( i, j - m2 );
		for( i = m1; i < d1; i++ )							// resting rows from aux, for all columns
		{
			for( j = 0; j < d2; j++ )
			{
				if(( i + (d1 - m1)) == ( j + (d2 - m2)))
					Vi( i, j ).set2const( 1.0 );			// set to the polynomial p(x) = 1
				else
					Vi( i, j ).set2zero();					// set to the polynomial p(x) = 0
			}
		}
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the matrix Vi (V0 in the initialization of the matrix sequence and
 * Vtildem in the while-cycle):
 * 
 * 		Vi = ( I  M )
 * 			 (    I )
 * 
 * Either the rows or the columns of Vi might be permuted. Thus, such matrices are possible:
 * 
 *		Vi = P * ( I  M )
 * 				 (    I )
 * 
 * with P the vector of rows permutations, or
 * 
 *		Vi = ( I  M ) * P^T
 * 			 (    I )
 * 
 * with P the vector of columns permutations.
 * 
 * \param[in] P The pointer to \a P, the vector of permutations over \a Vi. 
 * \param[in] col The binary parameter to indicate whether the rows or the columns of the
 * 		resulting matrix should be permuted (= 0, the rows; = 1, the columns).
 * \param[in] M The pointer to \a M, a matrix of \type T that is used in the computation. 
 * \param[out] Vi The pointer to \a Vi, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computeV( int *P, bool col, Matrix<T> &M, Matrix<T> &Vi )
{
	int i, j, v1 = Vi.nrows(), v2 = Vi.ncols(), m1 = M.nrows(), m2 = M.ncols();
	
	try
	{	
		int *PT = new int[ v2 ];							// P^T
		Matrix<T> aux( v1, v2, Vi.dimT() );					// auxiliary matrix
		
		aux.set2Id( m1, v2 - m2 );							// I
		for( i = 0; i < m1; i++ )							// first m1 rows from aux, for all columns
			for( j = v2 - m2; j < v2; j++ )					// resting columns: fill out with M
				aux( i, j ) = M( i, j - v2 + m2 );
		aux.set2zero( m1, v1, 0, v2 - m2 );					// 0
		aux.set2Id( m1, v1, v2 - m2, v2 );					// I
		//printm( v1, v2, aux, "Vi (before permuting)= \n", yellow );
		Vi = aux;
		if( col )											// consider COLUMNS PERMUTATIONS
		{
			//printiv( v2, P, "\nP = \n", magenta );
			//transpv( v2, P, PT );
			//printiv( v2, P, "\nP^T = \n", magenta );
			//Vi.cpermutem( PT );
			Vi.cpermutem( P, 1 );							// =1 to transpose P
		}
		else												// or ROWS PERMUTATIONS
			Vi.rpermutem( P );
		//printm( v1, v2, Vi, "Vi permuted = \n", yellow );
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the matrix Vi (V0 in the initialization of the matrix sequence and
 * Vtildem in the while-cycle):
 * 
 *		Vi = ( I  M ) * P^T
 * 			 (    I )
 * 
 * with P the vector of columns permutations.
 * 
 * \param[in] P The pointer to \a P, the vector of permutations over \a Vi. 
 * \param[in] M The pointer to \a M, a matrix of \type T that is used in the computation. 
 * \param[out] Vi The pointer to \a Vi, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computeV( int *P, Matrix<T> &M, Matrix<T> &Vi )
{
	int i, j, v1 = Vi.nrows(), v2 = Vi.ncols(), m1 = M.nrows(), m2 = M.ncols();
	
	try
	{	
		Matrix<T> aux( v1, v2, Vi.dimT() );					// auxiliary matrix
		aux.set2Id( m1, v2 - m2 );							// I
		for( i = 0; i < m1; i++ )							// first m1 rows from aux, for all columns
			for( j = v2 - m2; j < v2; j++ )					// resting columns: fill out with M
				aux( i, j ) = M( i, j - v2 + m2 );
		aux.set2zero( m1, v1, 0, v2 - m2 );					// 0
		aux.set2Id( m1, v1, v2 - m2, v2 );					// I
		//printm( v1, v2, aux, "Vi (before permuting)= \n", yellow );
		Vi = aux;
		Vi.cpermutem( P, 1 );								// columns permutation
		//printm( v1, v2, Vi, "Vi permuted = \n", yellow );
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the matrix Gm (i.e. G_i^{-1}):
 * 
 * Gm = V * ( S1m     m2    ) * Um
 *          ( m1   m1*S1*m2 )
 * 
 * with:
 * 		Gm : m-by-m matrix
 * 		V : m-by-m matrix
 * 		Um : m-by-m matrix
 * 		S1 : r-by-r matrix
 * 		S1m : r-by-r matrix
 * 		m1 : (m-r)-by-r matrix
 * 		m2 : r-by-(m-r) matrix
 * 
 * \param[in] m1 The pointer to \a m1, a matrix of \type T, free parameter. 
 * \param[in] m2 The pointer to \a m2, a matrix of \type T, free parameter. 
 * \param[in] V The pointer to \a V, a matrix of \type T that is used in the computation. 
 * \param[in] S1 The pointer to \a S1, a matrix of \type T that is used in the computation. 
 * \param[in] S1m The pointer to \a S1m, a matrix of \type T that is used in the computation. 
 * \param[in] U The pointer to \a U, a matrix of \type T that is used in the computation. 
 * \param[out] Gm The pointer to \a Gm, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computeGm( Matrix<T> &m1, Matrix<T> &m2, Matrix<T> &V, 
		Matrix<T> &S1, Matrix<T> &S1m, Matrix<T> &Um, Matrix<T> &Gm )
{
	int i, j, m = Gm.nrows(), r = S1m.ncols();
	
	try
	{	
		Matrix<T> aux1( m1.nrows(), S1.ncols(), Gm.dimT() ),// dimension for product m1*S1
			aux2( m1.nrows(), m2.ncols(), Gm.dimT() ),		// dimension for product m1*S1*m2
			aux3( m, m, Gm.dimT() );						// aux. matrix
		
		aux1.mmCaABbC( 1.0, 1.0, m1, S1 );					// m1 * S1
		aux2.mmCaABbC( 1.0, 1.0, aux1, m2 );				// * m2
		//printf( "\nMatrices involved in the calculation of Gm:\n" );
		//S1m.printm( "S1m = \n", green );
		//m1.printm( "m1 = \n", green );
		//m2.printm( "m2 = \n", green );
		//aux2.printm( "m1*S1*m2 = \n", green );
		for( i = 0; i < r; i++ )
		{
			for( j = 0; j < r; j++ )
				Gm( i, j ) = S1m( i, j );					// = S1m
			for( j = r; j < m; j++ )
				Gm( i, j ) = m2( i, j - r );				// = m2
		}
		for( i = r; i < m; i++ )
		{
			for( j = 0; j < r; j++ )
				Gm( i, j ) = m1( i - r, j );				// = m1
			for( j = r; j < m; j++ )
				Gm( i, j ) = aux2( i - r, j - r );			// = m1*S1*m2
		}
		//Gm.printm( "big matrix = \n", green );
		//Um.printm( "Um = \n", green );
		//V.printm( "V = \n", green );
		aux3.mmCaABbC( 1.0, 1.0, Gm, Um );					// * Um
		//aux3.printm( "big matrix * Um = \n", green );
		Gm.mmCaABbC( 1.0, 1.0, V, aux3 );					// V *
		//Gm.printm( "Gm = V * big matrix * Um = \n", green );
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the matrix mi1:
 * 
 * 		mi1 = ( y   m ) U^T
 * 
 * with:
 * 		y : m-by-m matrix
 * 		U : m-by-m matrix
 * 		m : m-by-m matrix
 * 		mi1 : r-by-r matrix
 * 
 * \param[in] y The pointer to \a y, a matrix of \type T. 
 * \param[in] U The pointer to \a U, a matrix of \type T. 
 * \param[in] m The pointer to \a m, a matrix of \type T, free parameter. 
 * \param[out] m1 The pointer to \a m1, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computemi1( Matrix<T> &y, Matrix<T> &U, Matrix<T> &m, Matrix<T> &mi1 )
{
	int i, j;
	
	try
	{	
		Matrix<T> aux1( y.nrows(), U.nrows(), U.dimT() ),	// U.nrows()=nr.cols. from U^T
			aux2( m.nrows(), U.nrows(), U.dimT() );
		
		//printf( "\naux1.nrows()=%d, aux1.ncols()=%d", aux1.nrows(), aux1.ncols() );
		//printf( "\naux2.nrows()=%d, aux2.ncols()=%d", aux2.nrows(), aux2.ncols() );
		aux1.mmCaABTbC( y.ncols(), 1, 1.0, 1.0, y, U );		// y * U1^T
		//aux1.printm( "aux1 = \n", yellow );
		//aux2.mmCaABTbC( m.ncols(), 0, 1.0, 1.0, m, U );	// m * U2^T
		//aux2.printm( "aux2 = \n", yellow );
		for( j = 0; j < U.nrows(); j++ )					// form mi1
		{
			for( i = 0; i < y.nrows(); i++ )
				mi1( i, j ) = aux1( i, j );
			//for( i = 0; i < m.nrows(); i++ )
			//	mi1( i + y.nrows(), j ) = aux2( i, j );
		}
		//mi1.printm( "mi1 = \n", yellow );
	} //try
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}

/**
 * Computes the matrix mi1:
 * 
 * 		mi1 = y * U^T
 * 
 * with:
 * 		y : m-by-m matrix
 * 		U : m-by-m matrix
 * 		mi1 : r-by-r matrix
 * 
 * \param[in] y The pointer to \a y, a matrix of \type T. 
 * \param[in] U The pointer to \a U, a matrix of \type T. 
 * \param[out] m1 The pointer to \a m1, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computemi1( Matrix<T> &y, Matrix<T> &U, Matrix<T> &mi1 )
{
	mi1.mmCaABTbC( 1.0, 1.0, y, U );						// y * U^T
	//mi1.printm( "mi1 = \n", yellow );
	return 0;
}

/**
 * Computes the matrix Dm (i.e. D_^{-1}):
 * 
 * 		Dm = U_D * BM * PI_D^T
 * 
 * where BM is of the form:
 * 
 * 		BM (or "big matrix") = 	( BM11  BM12 )
 * 								( BM21  BM22 )
 * 
 * with:
 * 
 * 		BM11 = Y^T after solving the system SD1 * Y = SD2 * m2D^T
 * 		BM12 = m2D
 * 		BM21 = m1D * SD1^T * Y^T
 * 		BM22 = m1D * SD1^T * m2D
 * 
 * and:
 * 
 * 		M11 = SA1 * SD1^T  +  SA2 * SD2^T
 * 		M12 = SA1 * 0  +  SA2 * I
 * 
 * in order to solve the system of equations:
 * 
 * 		M11 * m2D = M12
 * 
 * The most important matrix dimensions are:
 * 
 * 		SA : m-by-n matrix
 * 		SA1 : r-by-r matrix; SA2 : r-by-(n-r) matrix
 * 		SD : m-by-n matrix, calculated for D^T since the original D is a n-by-m matrix
 * 		SD1 : r-by-r matrix; SD2 : r-by-(n-r) matrix
 * 		PA : n vector
 * 		PD : n vector
 * 		UD : m-by-m matrix, calculated for D^T
 * 		Dm : m-by-n matrix
 * 		BM : m-by-n matrix
 * 		BM11 : r-by-r submatrix; BM12 : r-by-(n-r) submatrix
 * 		BM21 : (m-r)-by-r submatrix; BM22 : (m-r)-by-(n-r) submatrix
 * 		m1D : (m-r)-by-r matrix
 * 		m2D : r-by-(n-r) matrix
 * 		I (used in calculation of M12) : (n-r)-by-(n-r) matrix
 * 		I (used in calculation of BM11) : r-by-r matrix
 * 
 * \param[in] r The number of rows and columns that delimit submatrices.
 * \param[in] SA The pointer to \a SA, a matrix of \type T. 
 * \param[in] SD The pointer to \a SD, a matrix of \type T. 
 * \param[in] PA The pointer to \a PA, the vector of column permutations in matrix \a SA. 
 * \param[in] PD The pointer to \a PD, the vector of column permutations in matrix \a SD. 
 * \param[in] UD The pointer to \a UD, a matrix of \type T. 
 * \param[out] Dm The pointer to \a Dm, the resulting matrix, of \type T. 
 * \return The error code.
 * 
 */
template <typename T> int computeDm( int r, Matrix<T> &SA, Matrix<T> &SD, 
		int *PA, int *PD, Matrix<T> &UD, Matrix<T> &Dm )
{
	int i, j, m = Dm.nrows(), n = Dm.ncols(), dimT = Dm.dimT();
	
	try
	{	
		//printf( "\nEntering computeDm...\n" );
		Matrix<T> M11( n, n, dimT ), 						// M11[r][r] (max.nxn)
			M12( n, n, dimT ),								// M12[r][n-r] (max. nxn)
			SA1( m, n, dimT ), 								// SA1[r][r], submatrix of SA (max.mxn)
			SA2( m, n, dimT ),								// SA2[r][n-r], submatrix of SA (max.mxn)
			SD1( m, n, dimT ),								// SD1[r][r], submatrix of SD (max.mxn)
			SD2( m, n, dimT ),								// SD2[r][n-r], submatrix of SD (max.mxn)
			SL( n, n, dimT ), SR( n, n, dimT ),				// intermediate matrices
			m1D( m, n, dimT ),								// m1D[m-r][r], free parameter (max.mxn)
			aux( n, n, dimT ),								// aux[r][r] (max.nxn)
			aux1( m, m, dimT ),								// aux1[r][r] (max.mxm)
			aux2( m, m, dimT ),								// aux2[m-r][r] (max.mxm)
			aux3( m, m, dimT ), 							// aux3[m-r][r] (max.mxm)
			aux4( m, n, dimT ), 							// aux4[m-r][n-r] (max.mxn)
			aux5( m, n, dimT ),								// aux5[m][n] (max.mxn)
			I( r, r, dimT );

		if( r != 0 )										// SA1 and SD1 exist
		{
			SA1.setdim( r, r, dimT );
			SD1.setdim( r, r, dimT );
			M11.setdim( r, r, dimT );
		}
		if( n - r != 0 )									// SA2 and SD2 exist
		{
			SA2.setdim( r, n - r, dimT );
			SD2.setdim( r, n - r, dimT );
			M12.setdim( r, n - r, dimT );
		}
		decompS( r, SA, PA, SA1, SA2 );						// form SA1 and SA2 from S_A
		decompS( r, SD, PD, SD1, SD2 );						// form SD1 and SD2 from S_D
		formSL( r, SA1, SA2, SL );							// ( SA1  SA2 )
															// (  0    I  )
		formSR( r, SD1, SD2, SR );							// ( SD1^T  0 )
															// ( SD2^T  I )
		//
		// Mhat = ( SA1  SA2 )*PA^T * PD*( SD1^T  0 )
		//        (  0    I  )           ( SD2^T  I )
		//
		SL.cpermutem( PA, 1 );								// columns permutation, transposed
		SR.rpermutem( PD );									// rows permutation
		aux.mmCaABbC( 1.0, 1.0, SL, SR );					// obtain Mhat
		decompM( r, aux, M11, M12 );						// M21 and M22 are not interesting
		if( n - r != 0 )									// M12 exists
		{
			M11.utsolve( M12 );								// m2D = M11m*M12  <==>  M11*m2D = M12
															// on output: M12 is overwritten with 
															// the results (i.e. M12 = m2D)
		}
		//
		// Dm = U_D * BM * PI_D^T
		//
		if( m - r != 0 )									// m1D exists
		{
			aux2.setdim( m - r, r, dimT );					// new dim. aux2[m-r][r]
			m1D.setdim( m - r, r, dimT );					// new dim. m1D[m-r][r]
			m1D.set2zero();									// m1D = 0, free parameter
			aux2.mmCaABTbC( 1.0, 1.0, m1D, SD1 );			// aux2 = m1D*SD1^T
			//aux2.printdm( "\nm1D*SD1^T = \n", yellow );
		}
		else
		{
			// m1D, aux2 do not exist!!
		}
		I.set2Id();
		if( n - r != 0 )									// SA2 and SD2 exist
		{
			aux1.mmCaABTbC( 1.0, 1.0, SD2, M12 );			// aux1 = SD2*m2D^T
			//aux1.printdm( "\nSD2*m2D^T = \n", yellow );
			aux1 = I - aux1;
			//aux1.printdm( "\nI - SD2*m2D^T = \n", yellow );
		}
		else
		{
			aux1 = I;										// aux1 does not exist! Then, set it to I
			//aux1.printdm( "\nI = \n", yellow );
		}
		SD1.utsolve( aux1 );								// SD1 * Y = aux1
															// on output: aux1 is overwritten 
															// with the results
		//aux1.printdm( "\nY = \n", yellow );
		M11 = aux1.transpm();								// calc. Y^T using M11 as aux. matrix!
		//M11.printdm( "\nBM11 = Y^T = \n", yellow );
		aux5.setdim( m, n, dimT );							// new dim. aux5[m][n]
		if( m - r != 0 )									// aux2 = m1D*SD1^T exists, (m-r)xr
		{
			aux3.setdim( m - r, r, dimT );					// new dim. aux3[m-r][r]
			aux4.setdim( m - r, r, dimT );					// new dim. aux4[m-r][r]
			aux3.mmCaABbC( 1.0, 1.0, aux2, M11 );			// aux3 = BM21 = m1D*SD1^T * Y^T
			//aux3.printdm( "\nBM21 = m1D*SD1^T * Y^T = \n", yellow );
			if( n - r != 0 )								// m2D (i.e., M12) exists
			{
				aux4.mmCaABbC( 1.0, 1.0, aux2, M12 );		// aux4 = m1D*SD1^T * m2D
				//aux4.printdm( "\nBM22 = m1D*SD1^T * m2D = \n", yellow );
				// BM = ( BM11  BM12 ) = ( M11   M12 )
				//      ( BM21  BM22 )   ( aux3 aux4 )
				formBM( r, M11, M12, aux3, aux4, aux5 );	// form "big matrix"
			}
			else											// there is no m2D, aux4, BM22
			{
				// BM12 and BM22 do not exist and not will be needed to construct BM
				// BM = ( BM11 ) = ( M11  )
				//      ( BM21 )   ( aux3 )
				formBM( r, 1, M11, aux3, aux5 );			// form "big matrix", "vertically"=1
			}
		}
		else
		{
			// BM21 and BM22 do not exist and not will be needed to construct BM
			// BM = ( BM11  BM12 ) = ( M11  M12 )
			formBM( r, 0, M11, M12, aux5 );					// form "big matrix", "horizontally"=0
		}
		//aux5.printdm( "\nBM = \n", yellow );
		Dm.mmCaABbC( 1.0, 1.0, UD, aux5 );					// Dm = UD * BM
		//Dm.printdm( "\nDm = \n", yellow );
		Dm.cpermutem( PD, 1 );								// permute Dm
		//Dm.printdm( "\nDm (permuted) = \n", yellow );
		
	} //try



/***OLD VERSION***
		//
		// Initializations
		//
		//printf( "\nComputing Dm...\n" );
		//printf( "\nm = %d    n = %d    r = %d\n", m, n, r );
		SA1.setdim( r, r, dimT );							// new dim. SA1[r][r]
		SA2.setdim( r, n - r, dimT );						// new dim. SA2[r][n-r]
		SD1.setdim( r, r, dimT );							// new dim. SD1[r][r]
		SD2.setdim( r, n - r, dimT );						// new dim. SD2[r][n-r]
		decompS( r, SA, PA, SA1, SA2 );						// form SA1 and SA2 from S_A
		decompS( r, SD, PD, SD1, SD2 );						// form SD1 and SD2 from S_D
		//SA1.printtpm( "\nSA1 = \n", yellow );
		//SD1.printtpm( "\nSD1 = \n", yellow );
		//SA2.printtpm( "\nSA2 = \n", yellow );
		//SD2.printtpm( "\nSD2 = \n", yellow );
		aux1.setdim( r, r, dimT );							// new dim. aux1[r][r]
//Now: use SA1 columns permuted^T and SD1^T rows permuted!!!!!
		aux1.mmCaABTbC( 1.0, 1.0, SA1, SD1 );				// SA1 * SD1^T
		//aux1.printtpm( "\nSA1*SD1^T = \n", yellow );	
		M11.setdim( r, r, dimT );							// new dim. M11[r][r]
		M11.set2zero();
		M12.setdim( r, n - r, dimT );						// new dim. M12[r][n-r]
		M12.set2zero();
		if( n - r != 0 )									// SA2 and SD2 exist
		{
//Now: use SA2 columns permuted^T and SD2^T rows permuted!!!!!
			M11.mmCaABTbC( 1.0, 1.0, SA2, SD2 );			// SA2 * SD2^T
//Now: use SA2 columns permuted^T and verify that I was ok rows permuted!!!!!
			M12.mmCaAIbC( 1.0, 1.0, SA2, PD, 0 );			// SA2*I (I rows permuted, =0)
			M11 += aux1;									// M11 = SA1*SD1^T + SA2*SD2^T
			//M11.printtpm( "\nM11 = \n", yellow );
			//M12.printtpm( "\nM12 = \n", yellow );
			M11.utsolve( M12 );								// m2D = M11m*M12  <==>  M11*m2D = M12
															// on output: M12 is overwritten with 
															// the results (i.e. m2D)
			//M12.printtpm( "\nBM12 = m2D = \n", yellow );
		}
		else												// there is no SA2, SD2, M12, m2D
		{
			M11 += aux1;									// M11 = SA1*SD1^T
		}
		if( m - r != 0 )									// m1D exists
		{
			aux2.setdim( m - r, r, dimT );					// new dim. aux2[m-r][r]
			m1D.setdim( m - r, r, dimT );					// new dim. m1D[m-r][r]
			m1D.set2zero();									// m1D = 0, free parameter
			aux2.mmCaABTbC( 1.0, 1.0, m1D, SD1 );			// aux2 = m1D*SD1^T
			//aux2.printtpm( "\nm1D*SD1^T = \n", yellow );
		}
		else
		{
			// m1D, aux2 do not exist!!
		}
		Matrix<T> I( r, r, dimT );
		I.set2Id();
		if( n - r != 0 )									// SA2 and SD2 exist
		{
			aux1.mmCaABTbC( 1.0, 1.0, SD2, M12 );			// aux1 = SD2*m2D^T
			//aux1.printtpm( "\nSD2*m2D^T = \n", yellow );
			aux1 = I - aux1;
			//aux1.printtpm( "\nI - SD2*m2D^T = \n", yellow );
		}
		else
		{
			aux1 = I;										// aux1 does not exist! Then, set it to I
			//aux1.printtpm( "\nI = \n", yellow );
		}
		SD1.utsolve( aux1 );								// SD1 * Y = aux1
															// on output: aux1 is overwritten 
															// with the results
		//aux1.printtpm( "\nY = \n", yellow );
		M11 = aux1.transpm();								// calculate Y^T using auxiliary M11 !
		//M11.printtpm( "\nBM11 = Y^T = \n", yellow );
		aux5.setdim( m, n, dimT );							// new dim. aux5[m][n]
		if( m - r != 0 )									// aux2 = m1D*SD1^T exists, (m-r)xr
		{
			aux3.setdim( m - r, r, dimT );					// new dim. aux3[m-r][r]
			aux4.setdim( m - r, r, dimT );					// new dim. aux4[m-r][r]
			aux3.mmCaABbC( 1.0, 1.0, aux2, M11 );			// aux3 = BM21 = m1D*SD1^T * Y^T
			//aux3.printtpm( "\nBM21 = m1D*SD1^T * Y^T = \n", yellow );
			if( n - r != 0 )								// m2D (i.e., M12) exists
			{
				aux4.mmCaABbC( 1.0, 1.0, aux2, M12 );		// aux4 = m1D*SD1^T * m2D
				//aux4.printtpm( "\nBM22 = m1D*SD1^T * m2D = \n", yellow );
				// BM = ( BM11  BM12 ) = ( M11   M12 )
				//      ( BM21  BM22 )   ( aux3 aux4 )
				formBM( r, M11, M12, aux3, aux4, aux5 );	// form "big matrix"
			}
			else											// there is no m2D, aux4, BM22
			{
				// BM12 and BM22 do not exist and not will be needed to construct BM
				// BM = ( BM11 ) = ( M11  )
				//      ( BM21 )   ( aux3 )
				formBM( r, 1, M11, aux3, aux5 );			// form "big matrix", "vertically"=1
			}
		}
		else
		{
			// BM21 and BM22 do not exist and not will be needed to construct BM
			// BM = ( BM11  BM12 ) = ( M11  M12 )
			formBM( r, 0, M11, M12, aux5 );					// form "big matrix", "vertically"=0
		}
		//aux5.printtpm( "\nBM = \n", yellow );
		Dm.mmCaABbC( 1.0, 1.0, UD, aux5 );					// Dm = UD * BM
		//Dm.printtpm( "\nDm = \n", yellow );
		Dm.cpermutem( PD, 1 );								// permute Dm
		//Dm.printtpm( "\nDm (permuted) = \n", yellow );
		
	} //try
*/	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( out_of_range e )
	{
		printf( red );
		printf( "\n***Exception out_of_range found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 5;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	return 0;
}


int max( double *v, bool *selected, int *P, int a, int b );
int findIndex( int *v, int t );
template TPolyn dpc( int r, int m, TPolyn *x, TPolyn *y );
template int constructR<TPolyn>( Matrix<TPolyn> &A, int *piv, Matrix<TPolyn> &R, bool permuted );
template int hhlh<TPolyn>( int pos, TPolyn beta, Matrix<TPolyn> &A, int *piv, TPolyn *v );
template int hhlr<TPolyn>( int pos, TPolyn beta, Matrix<TPolyn> &A, int *piv, TPolyn *v );
template int hhr<TPolyn>( int rowpos, int colpos, TPolyn beta, Matrix<TPolyn> &A, TPolyn *v );
template int hm( int r, double eps, Matrix<TPolyn> &A, int *P, TPolyn *v, TPolyn beta );
template int hhv<TPolyn>( int r, Matrix<TPolyn> &A, int *P, TPolyn *v, double eps, TPolyn &beta );
template int hqrfcp<TPolyn>( double eps, Matrix<TPolyn> &A, int *P, Matrix<TPolyn> &Q, int& rank, 
		bool permuted, int printWhat, int printWhere, FILE * fn );
template int hqrfcp<TPolyn>( double eps, Matrix<TPolyn> &A, int *P, Matrix<TPolyn> &Q, 
		Matrix<TPolyn> &R, int& rank, bool permuted, int printWhat, int printWhere, FILE * fn );
template int pwq<TPolyn>( Matrix<TPolyn> &G, Matrix<TPolyn> &Gm, Matrix<TPolyn> &P, 
		Matrix<TPolyn> &W, Matrix<TPolyn> &Q );
template int decompS<TPolyn>( int r, Matrix<TPolyn> &R, int *P, Matrix<TPolyn> &R1, 
		Matrix<TPolyn> &R2 );
template int decompV<TPolyn>( int r, Matrix<TPolyn> &V, Matrix<TPolyn> &V1, Matrix<TPolyn> &V2 );
template int decompW<TPolyn>( int r, Matrix<TPolyn> &W, Matrix<TPolyn> &W1, Matrix<TPolyn> &W2 );
template int formBM<TPolyn>( int d, Matrix<TPolyn> &BM11, Matrix<TPolyn> &BM12, 
		Matrix<TPolyn> &BM21, Matrix<TPolyn> &BM22, Matrix<TPolyn> &BM );
template int formBM<TPolyn>( int d, bool vert, Matrix<TPolyn> &BM1, Matrix<TPolyn> &BM2, 
		Matrix<TPolyn> &BM );
template int filloutS<TPolyn>( int r, Matrix<TPolyn> &S1, Matrix<TPolyn> &S );
template int filloutM<TPolyn>( int d1, int d2, Matrix<TPolyn> &V, Matrix<TPolyn> &M );
template int computeTlim<TPolyn>( Matrix<TPolyn> &mi1, Matrix<TPolyn> &Si, 
		Matrix<TPolyn> &Vitilde, Matrix<TPolyn> &Tlim );
template int computemTui<TPolyn>( Matrix<TPolyn> &mi2, Matrix<TPolyn> &Si, 
		Matrix<TPolyn> &Uitilde, Matrix<TPolyn> &mTui );
template int computeV<TPolyn>( Matrix<TPolyn> &M, Matrix<TPolyn> &Vi );
template int computeV<TPolyn>( int *P, bool col, Matrix<TPolyn> &M, Matrix<TPolyn> &Vi );
template int computeV<TPolyn>( int *P, Matrix<TPolyn> &M, Matrix<TPolyn> &Vi );
template int computeGm<TPolyn>( Matrix<TPolyn> &m1, Matrix<TPolyn> &m2, Matrix<TPolyn> &V, 
		Matrix<TPolyn> &S1, Matrix<TPolyn> &S1m, Matrix<TPolyn> &Um, Matrix<TPolyn> &Gm );
template int computemi1( Matrix<TPolyn> &y, Matrix<TPolyn> &U, Matrix<TPolyn> &m, 
		Matrix<TPolyn> &mi1 );
template int computemi1( Matrix<TPolyn> &y, Matrix<TPolyn> &U, Matrix<TPolyn> &mi1 );
template int computeDm<TPolyn>( int r, Matrix<TPolyn> &SA, Matrix<TPolyn> &SD, 
		int *PA, int *PD, Matrix<TPolyn> &UD, Matrix<TPolyn> &Dm );
