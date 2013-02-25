#include "TaylorPLib.h"

using namespace std;
using namespace LibMatrix;

int Polynomial::unsetConstCount = 0;
int Polynomial::isConstCount = 0;


/**
 * Default constructor for the class. Creates the object.
 * 
 * It is a Taylor polynomial of order zero, i.e., it has a constant value:
 * 
 * 		p(x) = 1
 * 
 */
Polynomial::Polynomial():
	_constant(1),
	_order(0),
	_coeffs(0)
{
	allocateMemory(false);
	_coeffs[0] = 1.0;
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
Polynomial::Polynomial(int order, bool initialize):
	_constant(-1),
	_order(order),
	_coeffs(0)
{
	allocateMemory(initialize);
}

/**
 * Copy constructor.
 * 
 * 		Polynomial *newtp = new Polynomial((*tp));
 * 
 * \param[in] tp The Taylor polynomial object to copy from.
 * 
 */
Polynomial::Polynomial(const Polynomial &p):
	_constant(0),
	_order(-1),
	_coeffs(0)
{
	copyFrom(p);
}

/**
 * Short construtor for small hardcoded polynomials.
 *
 * \param[in] order Order of the new polynomial.
 * \param[in] constVal Constant part of the polynomial.
 * \param[in] values All non constant coefficients of the new polynomial (order ascending)
 *
 */
Polynomial::Polynomial(int order, double constVal, ...):
	_constant(order == 0 ? 1 : -1),
	_order(order),
	_coeffs(0)
{
	allocateMemory(false);

	va_list args;
	va_start(args, constVal);

	_coeffs[0] = constVal;
	for( int i = 1; i <= _order; i++ )
	{
		_coeffs[i] = va_arg(args, double);
	}	

	va_end(args);
}

/**
 * Short constructor for python port.
 * \param[in] coeffs All coefficients of the new polynomial, from low to high order and the constan coefficient first.
 *
 */
Polynomial::Polynomial(std::vector<double> coeffs):
	_constant(coeffs.size() == 0 ? 1 : -1),
	_order((int) coeffs.size() - 1),
	_coeffs(0)
{
	allocateMemory(false);

	for( int i = 0; i <= _order; i++)
	{
		_coeffs[i] = coeffs[i];
	}
}

/**
 * Destructor. Cleans up the object.
 * 
 */
Polynomial::~Polynomial()
{
	deallocateMemory();
}

/**
 * Gets a single coeffizient of the polynomial. Index 0 is the constant part of the polynomial.
 * 
 * \param[in] index The index of the coeffizient to get.
 * \return The value of the desired element.
 *
 */
double Polynomial::get(int index) const
{ 
	if (index < 0 || index > _order)
	{
		throw MathException("%d is not a valid coefficient index of this %d order polynomial.", index, _order);
	}
	return _coeffs[index];
}

/**
 * Sets a single coeffizient of the polynomial. Index 0 is the constant part of the polynomial.
 * 
 * \param[in] index The index of the coeffizient to get.
 * \return The value of the desired element.
 *
 */
void Polynomial::set(int index, double value)
{ 
	if (index < 0 || index > _order)
	{
		throw MathException("%d is not a valid coefficient index of this %d order polynomial.", index, _order);
	}
	
	_coeffs[index] = value;

	if (index > 0)
	{
		_constant = -1;
	}
}

/**
 * Implements the [] operator (subscript operator).
 * 
 * (double & Polynomial::... with '&' allows p[2] = 7, i.e., '[]' also in the left side!)
 * 
 * \param[in] index The index to be analyzed.
 * \return The coefficient at that index.
 * 
 */
double& Polynomial::operator[](int index)
{ 
	if (index < 0 || index > _order)
	{
		throw MathException("%d is not a valid coefficient index of this %d order polynomial.", index, _order);
	}
	
	if (index > 0)
	{
		// the constant state of the polynomial is no longer kown
		_constant = -1;
	}

	return _coeffs[index];
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
Polynomial Polynomial::operator=(const Polynomial &p)
{
	copyFrom(p);

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
bool Polynomial::operator==(const Polynomial &p) const
{
	if (_order != p._order)
	{
		return false;
	}

	for (int i = 0; i <= _order; i++)
	{
		if (_coeffs[i] != p._coeffs[i])
		{
			return false;
		}
	}
	
	return true;
}

/**
 * Implements the != operator.
 * Compares two Taylor polynomials.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the Taylor polynomials are different; 0 if they are equal.
 * 
 */
bool Polynomial::operator!=(const Polynomial &p) const
{
	return !(*this == p);
}

/**
 * Implements the < operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
bool Polynomial::operator<(const Polynomial &p)
{
	return _coeffs[0] < p._coeffs[0];
}

/**
 * Implements the <= operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
bool Polynomial::operator<=(const Polynomial &p)
{
	return _coeffs[0] <= p._coeffs[0];
}

/**
 * Implements the > operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
bool Polynomial::operator>(const Polynomial &p)
{
	return _coeffs[0] > p._coeffs[0];
}

/**
 * Implements the >= operator.
 * Compares two Taylor polynomials according to the value of the first coefficient.
 * 
 * \param[in] tp The pointer to the Taylor polynomial to be compared with.
 * \return 1 if the inequation holds; 0 otherwise.
 * 
 */
bool Polynomial::operator>=(const Polynomial &p)
{
	return _coeffs[0] >= p._coeffs[0];
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
Polynomial Polynomial::operator+(const Polynomial &p) const
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	// if one of the polynomial is constant we can optimize a bit
	// Polynomial v(_order, false);
	Polynomial v(_order);
	if (isConst())
	{
		v = p;
		v._coeffs[0] += _coeffs[0];
	}
	else if (p.isConst())
	{
		v = *this;
		v._coeffs[0] += p._coeffs[0];
	}
	else
	{
		// general case
		for (int i = 0; i <= _order; i++)
		{
			v._coeffs[i] = _coeffs[i] + p._coeffs[i];
		}
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
Polynomial Polynomial::operator+=(const Polynomial &p)
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	if (p.isConst())
	{
		_coeffs[0] += p._coeffs[0];
	}
	else
	{
		// constant state is undefined now
		_constant = -1;

		for (int i = 0; i <= _order; i++)
		{
			_coeffs[i] += p._coeffs[i];
		}
	}
	return *this;
}

/**
 * Implements the unary - operator.
 * 
 * \return A pointer to the resulting polynomial.
 * 
 */
Polynomial Polynomial::operator-() const
{
	Polynomial p(_order, false);

	for (int i = 0; i <= _order; i++)
	{
		p._coeffs[i] = - _coeffs[i];
	}
	
	return p;
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
Polynomial Polynomial::operator-(const Polynomial &p) const
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	// if one of the polynomial is constant we can optimize a bit
	Polynomial v(_order, false);
	if (isConst())
	{
		v = -p;
		v._coeffs[0] = _coeffs[0] - p._coeffs[0];
	}
	else if (p.isConst())
	{
		v = *this;
		v._coeffs[0] -= p._coeffs[0];
	}
	else
	{
		// general case
		for (int i = 0; i <= _order; i++)
		{
			v._coeffs[i] = _coeffs[i] - p._coeffs[i];
		}
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
Polynomial Polynomial::operator-=(const Polynomial &p)
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	if (p.isConst())
	{
		_coeffs[0] -= p._coeffs[0];
	}
	else
	{
		// constant state is undefined now
		_constant = -1;

		for (int k = 0; k <= _order; k++)
		{
			_coeffs[k] -= p._coeffs[k];
		}
	}
	return *this;
}

/**
 * Implements the * operator.
 * Multiplies a Taylor polynomial by a scalar.
 * 
 * \param[in] d The scalar value to multiply by.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::operator*(double d) const
{
	Polynomial v(_order, false);
	
	for (int i = 0; i <= _order; i++)
	{
		v._coeffs[i] = d * _coeffs[i];
	}

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
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::operator*(const Polynomial &p) const
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	// if one of the polynomial is constant we can optimize a bit
	// _coeffs[0] is the term with x^0

	// to calculate the real values incomment the lines with numbers in front
	// 1. Polynomial v(_order * 2);
	Polynomial v(_order);

	if (isConst())
	{
		v = p;
		for (int i = 0; i <= _order; i++)
		{
			v._coeffs[i] *= _coeffs[0];
		}
	}
	else if (p.isConst())
	{
		v = *this;
		for (int i = 0; i <= _order; i++)
		{
			v._coeffs[i] *= p._coeffs[0];
		}
	}
	else
	{
		// 2. for (int i = 0; i < _order * 2; i++)
		for (int i = 0; i <= _order; i++)
		{
			v._coeffs[i] = 0.0;

			for (int j = 0; j < i + 1; j++)
			{
  				v._coeffs[i] += _coeffs[j] * p._coeffs[i - j];
			}

		}
		// 3. v._coeffs[_order*2] = _coeffs[_order] * p._coeffs[_order];
	}
	return v;
}

/**
 * Implements the *= operator.
 * Multiplies a Taylor polynomial by a scalar.
 * 
 * \param[in] d The scalar value to multiply by.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::operator*=(double d)
{
	for (int i = 0; i <= _order; i++)
	{
		_coeffs[i] *= d;
	}

	return *this;
}

/**
 * Implements the *= operator.
 * Multiplies two Taylor polynomials.
 * 
 * \param[in] w The pointer to the polynomial to be multiplied by.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::operator*=(const Polynomial &p)
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	*this = *this * p;

	return *this;

}

/**
 * Implements the / operator.
 * Divides a Taylor polynomial by a scalar.
 * 
 * \param[in] d The scalar value to divide by.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::operator/(double d) const
{
	Polynomial v(_order, false);
	
	for (int i = 0; i <= _order; i++)
	{
		v._coeffs[i] = _coeffs[i] / d;
	}

	return v;
}

/**
 * Divides a Taylor polynomial (dividend) by another Taylor polynomial (divisor). 
 * Implements the / operator for Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = 1 / w_0 * [u_k - sum_{j=0}{k-1}  v_j * w_{k-j}]
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
Polynomial Polynomial::operator/(const Polynomial &p) const
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	Polynomial v(_order, false);

	for (int i = 0; i <= _order; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < i; j++)
		{
			sum += v._coeffs[j] * p._coeffs[i - j];
		}
		v._coeffs[i] = (_coeffs[i] - sum) / p._coeffs[0];
	}
	
	return v;
}

/**
 * Implements the /= operator.
 * Divides a Taylor polynomial by a scalar.
 * 
 * \param[in] d The scalar value to divide by.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::operator/=(double d)
{
	for (int i = 0; i <= _order; i++)
	{
		_coeffs[i] /= d;
	}

	return *this;
}

/**
 * Implements the /= operator.
 * Divides a Taylor polynomial (dividend) by another Taylor polynomial (divisor). 
 * 
 * \param[in] w The pointer to the divisor.
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::operator/=(const Polynomial &p)
{
	if (_order != p._order)
	{
		throw MathException("The order of both Taylor Polynoms should match.");
	}

	Polynomial v(_order, false);
	for (int i = 0; i <= _order; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < i; j++)
		{
			sum += v._coeffs[j] * p._coeffs[i - j];
		}
		v._coeffs[i] = (_coeffs[i] - sum) / p._coeffs[0];
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
 * \return The resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::sqr() const
{
	Polynomial v(_order, false);
	
	for (int i = 0; i <= _order; i++)
	{
		v._coeffs[i] = 0.0;
		for (int j = 0; j < i + 1; j++)
		{
			v._coeffs[i] += _coeffs[j] * _coeffs[i - j];
		}
	}

	return v;
}

/**
 * Sets this Taylor polynomial to its square.
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
 */
void Polynomial::setSqr()
{
	*this = sqr();
	/*
	Polynomial v(_order, false);
	
	for (int i = 0; i <= _order; i++)
	{
		v._coeffs[i] = 0.0;
		for (int j = 0; j < i + 1; j++)
		{
			v._coeffs[i] += _coeffs[j] * _coeffs[i - j];
		}
	}

	*this = v;
	*/
}

/**
 * Calculates the square root of a Taylor polynomial. Implements the square root function 
 * for Taylor arithmetic.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = 1 / 2*v_0 * [u_k - sum_{j=1}{k-1}  v_j * v_{k-j}]
 * 
 * for k = 1...d and v(t) = sqrt(u(t)), u and v being a Taylor polynomials, and d being 
 * the derivative degree. In particular, v_0 = sqrt(u_0).
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 * \return A pointer to the resulting Taylor polynomial.
 * 
 */
Polynomial Polynomial::sqrt() const
{
	Polynomial v(_order, false);
	
	v._coeffs[0] = ::sqrt(_coeffs[0]);
	
	for (int i = 1; i <= _order; i++)
	{
		double sum = 0.0;
		for (int j = 1; j < i; j++)
		{
			sum += v._coeffs[j] * v._coeffs[i - j];
		}
		v._coeffs[i] = (_coeffs[i] - sum) / (2 * v._coeffs[0]);
	}

	return v;
}

/**
 * Sets this Taylor polynomial to its square root.
 * 
 * The following coefficient propagation rule is applied:
 * 
 * 		v_k = 1 / 2*v_0 * [u_k - sum_{j=1}{k-1}  v_j * v_{k-j}]
 * 
 * for k = 1...d and v(t) = sqrt(u(t)), u and v being a Taylor polynomials, and d being 
 * the derivative degree. In particular, v_0 = sqrt(u_0).
 * 
 * (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques
 * of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, 
 * Philadelphia, PA, 2000)
 * 
 */
void Polynomial::setSqrt()
{
	*this = sqrt();
	/*
	Polynomial v(_order, false);
	
	v._coeffs[0] = ::sqrt(_coeffs[0]);
	
	for (int i = 1; i <= _order; i++)
	{
		double sum = 0.0;
		for (int j = 1; j < i; j++)
		{
			sum += v._coeffs[j] * v._coeffs[i - j];
		}
		v._coeffs[i] = (_coeffs[i] - sum) / (2 * v._coeffs[0]);
	}

	*this = v;
	*/
}

/**
 * Prints out the coefficients of a Taylor polynomial, starting by the independent term.
 * 
 */
void Polynomial::print()
{
	for (int i = 0; i < _order; i++)
	{
		printf("%.16lg\t", _coeffs[i]);
	}
	printf("%.16lg", _coeffs[_order]);
}

/**
 * Prints out to a file the coefficients of a Taylor polynomial, starting by the 
 * independent term.
 * 
 * \param[in] fn The output file to write the polynomial to.
 * 
 */
void Polynomial::print(FILE * fn)
{
	for (int i = 0; i < _order; i++)
	{
		fprintf(fn, "%.16lg\t", _coeffs[i]);
	}	
	fprintf(fn, "%.16lg", _coeffs[_order]);
}

/**
 * Evaluates a Taylor polynomial at a given value with a point of expansion.
 * 
 * \param[in] x The value to evaluate the polynomial at.
 * \param[in] alpha The point of expansion.
 * \return The result of the evaluation.
 * 
 */
double Polynomial::eval(double x, double alpha)
{
	double result = _coeffs[_order];
	double t = x - alpha;
	
	for (int i = _order - 1; i >= 0; i--)
	{
		result = t * result + _coeffs[i];
	}

	return result;
}

/**
 * Returns the first coefficient of the Taylor polynomial, i.e., the evaluation of the function.
 * 
 * \return The evaluation of the function at the initial point.
 * 
 */
double Polynomial::feval()
{
	return _coeffs[0];
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
void Polynomial::shift()
{
	for (int i = 1; i <= _order; i++)
	{
		_coeffs[i - 1] = i * _coeffs[i]; 
	}

	// set last coefficient to zero
	_coeffs[_order] = 0.0;
}

/**
 * Returns \a true in case it is a constant Taylor polynomial; 
 * \a false otherwise.
 * 
 * \return \a true if it is a constant Taylor polynomial; \a false otherwise.
 * 
 */
bool Polynomial::isConst() const
{
	isConstCount++;

	if (_constant == -1)
	{
		_constant = 1;
		for (int i = 1; i <= _order; i++)
		{
			if (_coeffs[i] != 0.0)
			{
				_constant = 0;
				return false;				
			}
		}
	}

	return _constant == 1;
}

/**
 * Returns \a true in case it is near a constant Taylor polynomial; 
 * \a false otherwise.
 * 
 * \param[in] eps The threshold value to compare with.
 * \return \a true if it is a constant Taylor polynomial; \a false otherwise.
 * 
 */
bool Polynomial::isConst(double eps)
{
	for (int i = 1; i <= _order; i++)
	{
		if (abs(_coeffs[i]) > eps)
		{
			return false;
		}
	}

	return true;
}

/**
 * Returns \a true in case it is a constant Taylor polynomial with value 1; 
 * \a false otherwise.
 * 
 * \return \a true if it is a constant Taylor polynomial with value 1; \a false otherwise.
 * 
 */
bool Polynomial::isId() const
{
	return (isConst()  &&  _coeffs[0] == 1.0);
}

/**
 * Returns \a true in case it is near a constant Taylor polynomial with value 1; 
 * \a false otherwise.
 * 
 * \param[in] eps The threshold value to compare with.
 * \return \a true if it is a constant Taylor polynomial with value 1; \a false otherwise.
 * 
 */
bool Polynomial::isId(double eps)
{
	return (isConst(eps)  &&  abs(_coeffs[0] - 1.0) < eps);
}

/**
 * Returns \a true in case all coefficients of the Taylor polynomial are zeroed; 
 * \a false otherwise.
 * 
 * \return \a true if all coefficients are zeroed; \a false otherwise.
 * 
 */
bool Polynomial::isZero() const
{
	return _coeffs[0] == 0 && isConst();
}

/**
 * Returns \a true in case all coefficients of the Taylor polynomial are lower or equal than
 * a threshold given as parameter; \a false otherwise.
 * 
 * \param[in] eps The threshold value to compare with.
 * \return \a true if all coefficients are almost null; \a false otherwise.
 * 
 */
bool Polynomial::isZero(double eps)
{
	for (int i = 0; i <= _order; i++)
	{
		if (abs(_coeffs[i]) > eps)
		{
			return false;
		}
	}
	
	return true;
}

/**
 * Sets a Taylor polynomial to the constant given as parameter.
 * 
 * \param[in] c The constant value of type \a double to set the Taylor polynomial to.
 * 
 */
void Polynomial::set2Const(double c)
{
	_coeffs[0] = c;

	for (int i = 1; i <= _order; i++)
	{
		_coeffs[i] = 0.0;
	}

	_constant = 1;
}

/**
 * Sets all coefficients of a Taylor polynomial to zero.
 * 
 */
void Polynomial::set2Zero()
{
	for (int i = 0; i <= _order; i++)	
	{
		_coeffs[i] = 0.0;
	}

	_constant = 1;
}

/**
 * Sets the coefficients of a Taylor polynomial to zero, from the order given as parameter on.
 * 
 * \param[in] ord Derivative order from which to start on (increasingly).
 * 
 */
void Polynomial::set2Zero(int ord)
{
	for (int i = ord; i <= _order; i++)
	{
		_coeffs[i] = 0.0;
	}
	
	_constant = ord < 2 ? 1 : -1;
}

/**
 * Sets the coefficients of a Taylor polynomial to the ones given as parameter.
 * 
 * \param[in] c A vector of coefficients of type \type double.
 * 
 */
void Polynomial::setCoeffs(double *coeffs)
{
	for (int i = 0; i <= _order; i++)
	{
		_coeffs[i] = coeffs[i];
	}
}


void Polynomial::setCoeffs(vector<double> coeffs)
{
	for (int i = 0; i <= _order; i++)
	{
		_coeffs[i] = coeffs[i];
	}
}


/***************
  P R I V A T E
  **************/

void Polynomial::allocateMemory(double *&coeffs, int order, bool initialize)
{
	if (order < 0)
	{
		throw MathException("The order of a polynomial cannot be negative.");
	}

	try
	{
		coeffs = new double[order + 1];
		if (coeffs == 0  ||  coeffs == NULL)
		{
			throw MathException("Memory allocation failure.");
		}
		if (initialize)
		{
			for (int i = 0; i <= order; i++)
			{
				coeffs[i] = 0.0;
			}
		}
	}
	catch (bad_alloc e)
	{
		throw MathException(e.what(), 4);
	}
	catch (...)
	{
		throw MathException("Error when allocating memory for a polynomial.");
	}
}

void Polynomial::deallocateMemory(double *&coeffs)
{
	delete [] coeffs;
}

void Polynomial::copyFrom(const Polynomial &p)
{
	if (p._order == _order)
	{
		// if p constant is defined (0 or 1) we can copy, otherwise it has to be undefined anyway
		_constant = p._constant;
		
		for (int i = 0; i <= _order; i++)
		{
			_coeffs[i] = p._coeffs[i];
		}
	}
	else
	{
		deallocateMemory();

		_order = p._order;
		_constant = p._constant;

		try
		{
			_coeffs = new double[_order + 1];
			if (_coeffs == 0  ||  _coeffs == NULL)
			{
				throw MathException("Memory allocation failure.");
			}
			for (int i = 0; i <= _order; i++)
			{
				_coeffs[i] = p._coeffs[i];
			}
		}
		catch (bad_alloc e)
		{
			throw MathException(e.what(), 4);
		}
		catch (...)
		{
			throw MathException("Error when allocating memory for a polynomial.");
		}
	}
}

void Polynomial::unsetConst()
{
	unsetConstCount++;
	_constant = false;
}

std::ostream& operator<<(std::ostream &out, const Polynomial &p)
{
	out << setiosflags(ios::fixed) << setprecision(1);

	out << '(' << p.get(0);
	for (int i = 1; i <= p.order(); i++)
 	{
 		out << " + " << p.get(i) << "x^" << i;
 	}
 	out<< ')';
	
	return out; 
}