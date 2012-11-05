/*!
 * \file TPolyn.cpp
 * \brief Class for Taylor arithmetic.
 * \author D. Monett Díaz
 * \date July, 2006
 * 
 * This class implements some functions used when working with Taylor arithmetic. 
 * (Taylor arithmetic is different from polynomial arithmetic!!)
 * 
 * It is used in the index determination in systems of ordinary differential
 * equations (ODEs) using Automatic/Algorithmic Differentiation (AD).
 * 
 * Included files:
 * 	\file TPolyn.h header file.
 * 
 * \todo
 * - not optimized for sparse polynomials
 * 
 */

#include "TPolyn.h"

using namespace std;


/**
 * Default constructor for the class. Creates the object.
 * 
 * It is a Taylor polynomial of order zero, i.e., it has a constant value:
 * 
 * 		p(x) = 1
 * 
 */
TPolyn::TPolyn()
{
	try
	{
		itsOrder = 0;
		_constant = false;
		itsCoeff = new double[ nrcoeff() ];					// coefficients
		if( (itsCoeff == 0)  ||  (itsCoeff == NULL))
			throw IDException( "Memory allocation failure.", 3 );
		itsCoeff[ 0 ] = 1;									// some initialization
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// other exceptions
	{
		IDException e( "Error when allocating a polynomial.", 32 );
        e.report();
		throw e.getErrCode();
	}
}

/**
 * Constructor for the class with a derivative order as parameter. Creates the object.
 * 
 * Example, order=3 <==> 4 coefficients:
 * 
 * 		p(x) = p_0 + p_1*x + p_2*x^2 + p_3*x^3
 * 
 * \param[in] order The derivative order of the Taylor polynomial.
 * 
 */
TPolyn::TPolyn( int order )
{
	try
	{
		itsOrder = order;
		_constant = false;
		itsCoeff = new double[ nrcoeff() ];					// coefficients
		if( (itsCoeff == 0)  ||  (itsCoeff == NULL))
			throw IDException( "Memory allocation failure.", 3 );
		for( int i = 0; i < nrcoeff(); i++ )
			itsCoeff[ i ] = 0.0;							// some initialization
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// other exceptions
	{
		IDException e( "Error when allocating a polynomial.", 32 );
        e.report();
		throw e.getErrCode();
	}
}

/**
 * Copy constructor.
 * 
 * 		TPolyn *newtp = new TPolyn( (*tp) );
 * 
 * \param[in] tp The Taylor polynomial object to copy from.
 * 
 */
TPolyn::TPolyn( const TPolyn &tp )
{
	try
	{
		itsOrder = tp.itsOrder;								// derivative order
		_constant = tp._constant;
		itsCoeff = new double[ nrcoeff() ];					// coefficients
		if( (itsCoeff == 0)  ||  (itsCoeff == NULL))
			throw IDException( "Memory allocation failure.", 3 );
		for( int i = 0; i < nrcoeff(); i++ )
			itsCoeff[ i ] = tp.itsCoeff[ i ];
	}
	catch( bad_alloc e )
	{
		printf( red );
		printf( "\n***Exception bad_alloc found:");
        printf( "\n***%s" , e.what() );
		printf( normal );
        throw 4;
	}
	catch( IDException e )
	{
        e.report();
		throw e.getErrCode();
	}
	catch(...)												// other exceptions
	{
		IDException e( "Error when allocating a polynomial.", 32 );
        e.report();
		throw e.getErrCode();
	}
}

/**
 * Destructor. Cleans up the object.
 * 
 */
TPolyn::~TPolyn()
{
	delete [] itsCoeff;
}

/**
 * Prints out the coefficients of a Taylor polynomial, starting by the independent term.
 * 
 */
void TPolyn::print()
{
	for( int i = 0; i < order(); i++ )
		printf( "%.16lg%c", itsCoeff[ i ], '\t' );
		//printf( " %.16lg%s%d %c", itsCoeff[ i ], "*x^", i, '+' );
	printf( "%.16lg", itsCoeff[ order() ] );
	//printf( " %.16lg%s%d", itsCoeff[ order() ], "*x^", order() );
}

/**
 * Prints out to a file the coefficients of a Taylor polynomial, starting by the 
 * independent term.
 * 
 * \param[in] fn The output file to write the polynomial to.
 * 
 */
void TPolyn::print( FILE * fn )
{
	for( int i = 0; i < order(); i++ )
		fprintf( fn, "%.16lg%c", itsCoeff[ i ], '\t' );
		//fprintf( fn, " %.16lg%s%d %c", itsCoeff[ i ], "*x^", i, '+' );
	fprintf( fn, "%.16lg", itsCoeff[ order() ] );
	//fprintf( fn, " %.16lg%s%d", itsCoeff[ order() ], "*x^", order() );
}

/**
 * Evaluates a Taylor polynomial at a given value with a point of expansion.
 * 
 * \param[in] x The value to evaluate the polynomial at.
 * \param[in] alpha The point of expansion.
 * \return The result of the evaluation.
 * 
 */
double TPolyn::eval( double x, double alpha )
{
	double result = itsCoeff[ order() ],
		t = x - alpha;
	
	for( int i = order() - 1; i >= 0; i-- )
		result = t * result + itsCoeff[ i ];
	
	return result;
}

/**
 * Returns the first coefficient of the Taylor polynomial, i.e., the evaluation of the function.
 * 
 * \return The evaluation of the function at the initial point.
 * 
 */
double TPolyn::feval()
{
	return itsCoeff[ 0 ];
}

/**
 * Implements the SHIFT operator to calculate the derivative of a Taylor polynomial.
 * 
 * The new coefficients are shifted to the left and the last one is zeroed.
 * 
 * E.g.:
 * 		y(t) = sum_{j=0}^{d} y_j * t^j + O(t^d+1)
 * 			 = y_0 + y_1*t + y_2*t^2 + ... + y_d*t^d
 * 
 * 		y'(t) = y_1 + 2*y_2*t + 3*y_3*t^2 + ... + d*y_d*t^d-1 + 0
 * 
 */
void TPolyn::shift()
{
	for( int i = 1; i < nrcoeff(); i++ )
		itsCoeff[ i - 1 ] = itsCoeff[ i ] * i; 
	itsCoeff[ nrcoeff() - 1 ] = 0.0;						// new last coeff. is zeroed
}

/**
 * Implements the [] operator (subscript operator).
 * 
 * (double & TPolyn::... with '&' allows p[2] = 7, i.e., '[ ]' also in the left side!!)
 * 
 * \param[in] index The index to be analyzed.
 * \return The coefficient at that index.
 * 
 */
double & TPolyn::operator[]( int index )
{ 
	return itsCoeff[ index ];
}

/**
 * Assignment operator.
 * 
 * 		newtp = tp;
 * 
 * \param[in] tp The Taylor polynomial object to assign to.
 * \return The resulting polynomial.
 * 
 */
TPolyn TPolyn::operator=( const TPolyn &tp )
{
	int i;
	
	if( this == &tp )
  		return *this;
	delete [] itsCoeff;										// deallocate the old Taylor polynomial
	itsOrder = tp.itsOrder;									// derivative order
	_constant = tp._constant;
	itsCoeff = new double[ nrcoeff() ];						// coefficients
	for( i = 0; i < nrcoeff(); i++ )
		itsCoeff[ i ] = tp.itsCoeff[ i ]; 
	
	return *this;
}

/**
 * Implements the == operator.
 * Compares two Taylor polynomials.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the Taylor polynomials are equal; 0 if they are different.
 * 
 */
int TPolyn::operator==( const TPolyn &tp )
{
	bool idem = true;
	int i;
	
	if( itsOrder == tp.itsOrder )
	{
		for( i = 0; i < nrcoeff(); i++ )
		{
			if( itsCoeff[ i ] != tp.itsCoeff[ i ] )
				idem = false;
		}
	}
	else
	{
		idem = false;
	}
	
	return idem;
}

/**
 * Implements the != operator.
 * Compares two Taylor polynomials.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the Taylor polynomials are different; 0 if they are equal.
 * 
 */
int TPolyn::operator!=( const TPolyn &tp )
{
	bool diff = false;
	int i;
	
	if( itsOrder == tp.itsOrder )
	{
		for( i = 0; i < nrcoeff(); i++ )
		{
			if( itsCoeff[ i ] == tp.itsCoeff[ i ] )
				diff = false;
		}
	}
	else
	{
		diff = true;
	}
	
	return diff;
}

/**
 * Implements the < operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
int TPolyn::operator<( const TPolyn &tp )
{
	return( itsCoeff[ 0 ] < tp.itsCoeff[ 0 ] );
}

/**
 * Implements the <= operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
int TPolyn::operator<=( const TPolyn &tp )
{
	return( itsCoeff[ 0 ] <= tp.itsCoeff[ 0 ] );
}

/**
 * Implements the > operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
int TPolyn::operator>( const TPolyn &tp )
{
	return( itsCoeff[ 0 ] > tp.itsCoeff[ 0 ] );
}

/**
 * Implements the >= operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
int TPolyn::operator>=( const TPolyn &tp )
{
	return( itsCoeff[ 0 ] >= tp.itsCoeff[ 0 ] );
}

/**
 * Adds up two Taylor polynomials. Implements the + operator for Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = u_k + w_k
 * 
 * for k = 1...d and v(t) = u(t) + w(t), u, v, w being Taylor polynomials, and d being 
 * the derivative degree.
 * 
 * It is assumed that all three Taylor polynomials have the same derivative degree d.
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 * \param[in] w The pointer to the Taylor polynomial to be added to.
 * \return A pointer to the resulting polynomial.
 * 
 */
TPolyn TPolyn::operator+( const TPolyn &w )
{
	TPolyn v( itsOrder ), aux( itsOrder );

	aux = w;
	if( isConst() )											// 'this' is a constant polynomial
	{
		v.itsCoeff[ 0 ] = itsCoeff[ 0 ] + w.itsCoeff[ 0 ];
		for( int k = 1; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = w.itsCoeff[ k ];
	}
	else if( aux.isConst() )								// 'w' is a constant polynomial
	{
		v.itsCoeff[ 0 ] = itsCoeff[ 0 ] + w.itsCoeff[ 0 ];
		for( int k = 1; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ k ];
	}
	else													// general case
	{
		for( int k = 0; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ k ] + w.itsCoeff[ k ];
	}

	return v;
}

/**
 * Implements the += operator.
 * Substracts a Taylor polynomial from the current one using pointers to arrays that store 
 * the coefficients.
 * 
 * \param[in] w The pointer to the Taylor polynomial to be added to.
 * \return A pointer to the resulting polynomial.
 * 
 */
TPolyn TPolyn::operator+=( const TPolyn &w )
{
	TPolyn aux( itsOrder );

	aux = w;
	if( aux.isConst() )										// 'w' is a constant polynomial
		itsCoeff[ 0 ] += w.itsCoeff[ 0 ];
	else													// 'this' is a constant polynomial...
		for( int k = 0; k < nrcoeff(); k++ )				// ...or general case
			itsCoeff[ k ] += w.itsCoeff[ k ];

	return *this;
}

/**
 * Substracts two Taylor polynomials. Implements the - operator for Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = u_k - c*w_k
 * 
 * for k = 1...d and v(t) = u(t) - c*w(t), c being a real constant, u, v, w being
 * Taylor polynomials, and d being the derivative degree.
 * 
 * It is assumed that all three Taylor polynomials have the same derivative degree d.
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 * \param[in] w The pointer to the Taylor polynomial to be substracted from.
 * \return A pointer to the resulting polynomial.
 * 
 */
TPolyn TPolyn::operator-( const TPolyn &w )
{
	TPolyn v( itsOrder ), aux( itsOrder );

	aux = w;
	if( isConst() )											// 'this' is a constant polynomial
	{
		v.itsCoeff[ 0 ] = itsCoeff[ 0 ] - w.itsCoeff[ 0 ];
		for( int k = 1; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = - w.itsCoeff[ k ];
	}
	else if( aux.isConst() )								// 'w' is a constant polynomial
	{
		v.itsCoeff[ 0 ] = itsCoeff[ 0 ] - w.itsCoeff[ 0 ];
		for( int k = 1; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ k ];
	}
	else													// general case
	{
		for( int k = 0; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ k ] - w.itsCoeff[ k ];
	}
	
	return v;
}

/**
 * Implements the -= operator.
 * Substracts a Taylor polynomial from the current one using pointers to arrays that 
 * store the coefficients.
 * 
 * \param[in] w The pointer to the Taylor polynomial to be substracted from.
 * \return A pointer to the resulting polynomial.
 * 
 */
TPolyn TPolyn::operator-=( const TPolyn &w )
{
	TPolyn aux( itsOrder );
	
	aux = w;
	if( aux.isConst() )										// 'w' is a constant polynomial
		itsCoeff[ 0 ] -= w.itsCoeff[ 0 ];
	else													// 'this' is a constant polynomial...
		for( int k = 0; k < nrcoeff(); k++ )				// ...or general case
			itsCoeff[ k ] -= w.itsCoeff[ k ];

	return *this;
}

/**
 * Implements the unary - operator.
 * 
 * \return A pointer to the resulting polynomial.
 * 
 */
TPolyn TPolyn::operator-()
{
	TPolyn v( itsOrder );
	for( int k = 0; k < nrcoeff(); k++ )
		v.itsCoeff[ k ] = - itsCoeff[ k ];
	
	return v;
}

/**
 * Multiplies two Taylor polynomials. Implements the * operator for Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = sum_{j=0}{k}  u_j * w_{k-j}
 * 
 * for k = 1...d and v(t) = u(t) * w(t), u, v, w being Taylor polynomials, and d being 
 * the derivative degree.
 * 
 * It is assumed that all three Taylor polynomials have the same derivative degree d.
 * 
 * Three different cases are distinguished here: when at least one of the polynomials
 * is a constant polynomial and when both polynomials are not.
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 * \param[in] w The pointer to the Taylor polynomial to be multiplied by.
 * \return A pointer to the resulting polynomial.
 * 
 */
TPolyn TPolyn::operator*( const TPolyn &w )
{
	TPolyn v( itsOrder ), aux( itsOrder );
	
	aux = w;
	if( isConst() )											// 'this' is a constant polynomial
	{
		for( int k = 0; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ 0 ] * w.itsCoeff[ k ];
	}
	else if( aux.isConst() )								// 'w' is a constant polynomial
	{
		for( int k = 0; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ k ] * w.itsCoeff[ 0 ];
	}
	else													// general case
	{
		for( int k = 0; k < nrcoeff(); k++ )
		{
			v.itsCoeff[ k ] = 0.0;
			for( int j = 0; j < k + 1; j++ )
				v.itsCoeff[ k ] += itsCoeff[ j ] * w.itsCoeff[ k - j ];
		}
	}
	
	return v;
}

/**
 * Implements the *= operator.
 * Multiplies two Taylor polynomials.
 * 
 * \param[in] w The pointer to the polynomial to be multiplied by.
 * \return A pointer to the resulting polynomial.
 * 
 */
TPolyn TPolyn::operator*=( const TPolyn &w )
{
	TPolyn v( itsOrder ), aux( itsOrder );

	aux = w;
	if( isConst() )											// 'this' is a constant polynomial
	{
		for( int k = 0; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ 0 ] * w.itsCoeff[ k ];
		_constant = false;
	}
	else if( aux.isConst() )								// 'w' is a constant polynomial
	{
		for( int k = 0; k < nrcoeff(); k++ )
			v.itsCoeff[ k ] = itsCoeff[ k ] * w.itsCoeff[ 0 ];
	}
	else													// general case
	{
		for( int k = 0; k < nrcoeff(); k++ )
		{
			v.itsCoeff[ k ] = 0.0;
			for( int j = 0; j < k + 1; j++ )
				v.itsCoeff[ k ] += itsCoeff[ j ] * w.itsCoeff[ k - j ];
		}
	}
	*this = v;
	
	return *this;
}

/**
 * Implements the * operator.
 * Multiplies a Taylor polynomial by a scalar.
 * 
 * \param[in] alpha The scalar value to multiply by.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
TPolyn TPolyn::operator*( double alpha )
{
	TPolyn v( itsOrder );
	for( int j = 0; j < nrcoeff(); j++ )
		v.itsCoeff[ j ] = itsCoeff[ j ] * alpha;
	
	return v;
}

/**
 * Divides a Taylor polynomial (dividend) by another Taylor polynomial (divisor). 
 * Implements the / operator for Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = 1 / w_0 * [ u_k - sum_{j=0}{k-1}  v_j * w_{k-j} ]
 * 
 * for k = 1...d and v(t) = u(t) / w(t), u, v, w being Taylor polynomials, and d being 
 * the derivative degree.
 * 
 * It is assumed that all three Taylor polynomials have the same derivative degree d.
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 * \param[in] w The pointer to the divisor.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
TPolyn TPolyn::operator/( const TPolyn &w )
{
	double sum;
	TPolyn v( itsOrder );
	for( int k = 0; k < nrcoeff(); k++ )
	{
		sum = 0.0;
		for( int j = 0; j < k; j++ )
			sum += v.itsCoeff[ j ] * w.itsCoeff[ k - j ];
		v.itsCoeff[ k ] = ( itsCoeff[ k ] - sum ) / w.itsCoeff[ 0 ];
	}
	
	return v;
}

/**
 * Implements the /= operator.
 * Divides a Taylor polynomial (dividend) by another Taylor polynomial (divisor). 
 * 
 * \param[in] w The pointer to the divisor.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
TPolyn TPolyn::operator/=( const TPolyn &w )
{
	double sum;
	TPolyn v( itsOrder );
	
	for( int k = 0; k < nrcoeff(); k++ )
	{
		sum = 0.0;
		for( int j = 0; j < k; j++ )
			sum += v.itsCoeff[ j ] * w.itsCoeff[ k - j ];
		v.itsCoeff[ k ] = ( itsCoeff[ k ] - sum ) / w.itsCoeff[ 0 ];
	}
	*this = v;

	return *this;
}

/**
 * Calculates the square of a Taylor polynomial. Implements the square function for 
 * Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = sum_{j=0}{k}  u_j * u_{k-j}
 * 
 * for k = 1...d and v(t) = u(t)^2, u and v being a Taylor polynomials, and d being 
 * the derivative degree.
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
TPolyn TPolyn::sqr()
{
	TPolyn v( itsOrder );
	for( int k = 0; k < nrcoeff(); k++ )
	{
		v.itsCoeff[ k ] = 0.0;
		for( int j = 0; j < k + 1; j++ )
			v.itsCoeff[ k ] += itsCoeff[ j ] * itsCoeff[ k - j ];
	}
	*this = v;
	
	return *this;
}

/**
 * Calculates the square root of a Taylor polynomial. Implements the square root function 
 * for Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = 1 / 2*v_0 * [ u_k - sum_{j=1}{k-1}  v_j * v_{k-j} ]
 * 
 * for k = 1...d and v(t) = sqrt(u(t)), u and v being a Taylor polynomials, and d being 
 * the derivative degree. In particular, v_0 = sqrt( u_0 ).
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
TPolyn TPolyn::sqrt()
{
	double sum;
	TPolyn v( itsOrder );
	v.itsCoeff[ 0 ] = ::sqrt( itsCoeff[ 0 ] );				// '::' to use sqrt of math.h
	for( int k = 1; k < nrcoeff(); k++ )
	{
		sum = 0.0;
		for( int j = 1; j < k; j++ )
			sum += v.itsCoeff[ j ] * v.itsCoeff[ k - j ];
		v.itsCoeff[ k ] = ( itsCoeff[ k ] - sum ) / (2*v.itsCoeff[ 0 ]);
	}
	*this = v;

	return *this;
}

/**
 * Returns \a true in case it is a constant Taylor polynomial; 
 * \a false otherwise.
 * 
 * \return \a true if it is a constant Taylor polynomial; \a false otherwise.
 * 
 */
bool TPolyn::isConst()
{
	bool c = true;

	for( int i = 1; i < nrcoeff(); i++ )
	{
		if( itsCoeff[ i ] != 0.0 )							// at least for one coefficient
		{													// holds p_i <> 0
			c = false;
			break;
		}
		if( !c )
			break;
	}
	if( c )													// it is constant!
		if( !_constant )
			_constant = true;
			
	return ( c );
}

/**
 * Returns \a true in case it is near a constant Taylor polynomial; 
 * \a false otherwise.
 * 
 * \param[in] eps The threshold value to compare with.
 * \return \a true if it is a constant Taylor polynomial; \a false otherwise.
 * 
 */
bool TPolyn::isConst( double eps )
{
	bool c = true;

	for( int i = 1; i < nrcoeff(); i++ )
	{
		if( fabs( itsCoeff[ i ] ) > eps )					// at least for one coefficient
		{													// holds p_i <> 0
			c = false;
			break;
		}
		if( !c )
			break;
	}
	if( c )													// it is constant!
		if( !_constant )
			_constant = true;
			
	return ( c );
}

/**
 * Returns \a true in case it is a constant Taylor polynomial with value 1; 
 * \a false otherwise.
 * 
 * \return \a true if it is a constant Taylor polynomial with value 1; \a false otherwise.
 * 
 */
bool TPolyn::isId()
{
	return ( isConst()  &&  itsCoeff[ 0 ] == 1.0 );
}

/**
 * Returns \a true in case it is near a constant Taylor polynomial with value 1; 
 * \a false otherwise.
 * 
 * \param[in] eps The threshold value to compare with.
 * \return \a true if it is a constant Taylor polynomial with value 1; \a false otherwise.
 * 
 */
bool TPolyn::isId( double eps )
{
	return ( isConst( eps )  &&  fabs( itsCoeff[ 0 ] - 1.0 ) < eps );
}

/**
 * Returns \a true in case all coefficients of the Taylor polynomial are zeroed; 
 * \a false otherwise.
 * 
 * \return \a true if all coefficients are zeroed; \a false otherwise.
 * 
 */
bool TPolyn::isZero()
{
	bool z = true;

	for( int i = 0; i < nrcoeff(); i++ )
	{
		if( itsCoeff[ i ] != 0.0 )							// at least for one coefficient
		{													// holds p_i <> 0
			z = false;
			break;
		}
		if( !z )
			break;
	}
	
	return z;
}

/**
 * Returns \a true in case all coefficients of the Taylor polynomial are lower or equal than
 * a threshold given as parameter; \a false otherwise.
 * 
 * \param[in] eps The threshold value to compare with.
 * \return \a true if all coefficients are almost null; \a false otherwise.
 * 
 */
bool TPolyn::isZero( double eps )
{
	bool z = true;

	for( int i = 0; i < nrcoeff(); i++ )
	{
		if( fabs( itsCoeff[ i ] ) > eps )					// at least for one coefficient
		{													// holds |p_i| > eps
			z = false;
			break;
		}
	if( !z )
		break;
	}
	
	return z;
}

/**
 * Sets all coefficients of a Taylor polynomial to zero.
 * 
 * \return The error code.
 * 
 */
int TPolyn::set2zero()
{
	for( int i = 0; i < nrcoeff(); i++ )					// p(x) = 0
		itsCoeff[ i ] = 0.0;
	_constant = true;
	
	return 0;
}

/**
 * Sets the coefficients of a Taylor polynomial to zero, from the order given as parameter on.
 * 
 * \param[in] ord Derivative order from which to start on (increasingly).
 * \return The error code.
 * 
 */
int TPolyn::set2zero( int ord )
{
	for( int i = ord; i < nrcoeff(); i++ )
		itsCoeff[ i ] = 0.0;								// p_i = 0, for i >= ord
	
	return 0;
}

/**
 * Sets a Taylor polynomial to the constant given as parameter.
 * 
 * \param[in] c The constant value of type \a double to set the Taylor polynomial to.
 * \return The error code.
 * 
 */
int TPolyn::set2const( double c )
{
	itsCoeff[ 0 ] = c;										// p(x) = c
	set2zero( 1 );											// set to zero the resting coefficients
	_constant = true;
	
	return 0;
}

/**
 * Sets the coefficients of a Taylor polynomial to the ones given as parameter.
 * 
 * \param[in] c A vector of coefficients of type \type double.
 * \return The error code.
 * 
 */
int TPolyn::setCoeffs( double *c )
{
	for( int i = 0; i < nrcoeff(); i++ )
		itsCoeff[ i ] = c [ i ];

	return 0;
}
