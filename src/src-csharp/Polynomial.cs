using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace LibMatrix
{
    /// <summary>
    ///
    /// (Taylor) Polynomial with derivate degree n (n+1 coefficients):<para/>
    /// <para/>
    /// P_n(x) = f(a) + (x-a)*f'(a)/1! + (x-a)^2*f''(a)/2! + (x-a)^3*f'''(a)/3! <para/>
    ///			+ ... + (x-a)^n*f^{(n)}(a)/n!<para/>
    /// <para/>
    ///		  = sum{k=0}{n} (x-a)^k*f^{(k)}(a)/k!<para/>
    /// <para/>
    /// being 'f' the function to be approximated by P_n(x) at point 'a', with its<para/>
    /// first n derivatives existing on a closed interval I, so that<para/>
    /// <para/>
    /// f(x) = P_n(x) + R_n(x)<para/>
    /// <para/>
    /// the remainder term being R_n(x) = (x-a)^{n+1}*f^{(n+1)}(c)/(n+1)! for some<para/>
    /// 'c' between 'x' and 'a'.<para/>
    /// <para/>
    /// Another often used form:<para/>
    /// <para/>
    /// f(x0+h) = f(x0) + h*f'(x0)/1! + h^2*f''(x0)/2! + h^3*f'''(x0)/3! <para/>
    ///			+ ... + h^n*f^{(n)}(x0)/n!<para/>
    /// <para/>
    /// </summary>
    public class Polynomial
    {
        #region private vars

        private int _constant;
        private int _order;
        private double[] _coeffs;

        private int unsetConstCount = 0,
                    isConstCount    = 0;

        #endregion

        #region public getter and setter

        /// <summary>
        /// returns order<para/>
        /// </summary>
        public int order {
            get { return _order; }
        }

        /// <summary>
        /// returns the number of coeffs<para/>
        /// </summary>
        public int ncoeff {
            get { return _order + 1; }
        }

        #endregion

        #region ctors

        /// <summary>
        /// Default constructor for the class. Creates the object.<para/>
        /// It is a Taylor polynomial of order zero, i.e., it has a constant value: <para/>
        /// <para/>
        /// p(x) = 1
        /// </summary>
        public Polynomial() : this( 0, new double[] { 1 }, 1)
        {
        }

        /// <summary>
        /// Constructor for the class with a derivative order as parameter. Creates the object.<para/>
        /// Example, order=3 &lt;==/&gt; 4 coefficients:<para/>
        /// <para/>
        /// p(x) = p_0 + p_1*x + p_2*x^2 + p_3*x^3<para/>
        /// </summary>
        /// <param name="order">The derivative order of the Taylor polynomial.</param>
        public Polynomial(int order) : this(order, new double[order + 1])
        {
        }

        /// <summary>
        /// Constructor for the class with given params<para/>
        /// </summary>
        /// <param name="order">number of items (zerobased)</param>
        /// <param name="coeffs">must match order<para/>
        /// e.g.: order = 3 -> coeffs = new double[4] with values</param>
        /// <param name="constant">[optional] 1 = constant,<para /> zero = not constant<para />, -1 = unknown<para/> </param>
        public Polynomial(int order, double[] coeffs, int constant = -1)
        {
            initializePolynomial(order, coeffs, constant);
        }

        /// <summary>
        /// Copy constructor.<para/>
        /// Polynomial newtp = new Polynomial(oldtp);<para/>
        /// Equivalent to Polynomial newtp = oldtp;
        /// </summary>
        /// <param name="p">The Taylor polynomial object to copy from.</param>
        public Polynomial(Polynomial p)
        {
            initializePolynomial(p.order, p._coeffs);
        }

        /// <summary>
        /// Destructor. Cleans up the object.
        /// </summary>
        ~Polynomial()
        {
            _order = 0;
            _constant = 0;
            _coeffs = new double[] { 1 };
        }

        #endregion

        #region private functions

        /// <summary>
        /// fior unseting the constCount
        /// </summary>
        private void unsetConst()
        {
            unsetConstCount++;
            _constant = 0;
        }

        /// <summary>
        /// for initializing the Polynomial<para/>
        /// </summary>
        /// <param name="order">The order of Polynomial</param>
        /// <param name="coeffs">The Values </param>
        /// <param name="constant">default = -1</param>
        private void initializePolynomial(int order, double[] coeffs, int constant = -1)
        {
            this._order = order;
            this._constant = constant;

            // this must be done, if a copy of a Polynomial is created.
            // e.g. Polynomial b = new Polynomial(a) where a is a Polynomial
            // if it would not copied, for example b[0] = 7 would mean, that also a[0] = 7
            double[] tempCoeffs = new double[order + 1];
            for (int i = 0; i <= order; i++)
                tempCoeffs[i] = coeffs[i];
            this._coeffs = tempCoeffs;
        }

        #endregion

        #region public functions

        #region Operators

        /// <summary>
        /// poly[2] = 7, i.e., '[]' also in the left side!)
        /// </summary>
        /// <param name="i">The index to be analyzed.</param>
        /// <returns>The coefficient at that index.</returns>
        public double this[int i]
        {
            get
            {
                if (i > _order || i < 0)
                    throw new MathException("Invalid coefficient index of polynomial.");
                else
                    return _coeffs[i];
            }
            set
            {
                if (i > _order || i < 0)
                    throw new MathException("Invalid coefficient index of polynomial.");
                else
                {
                    _coeffs[i] = value;
                    _constant = -1;
                }

            }
        }

        /// <summary>
        /// Implements the == operator.<para/>
        /// Compares two Taylor polynomials.<para/>
        /// </summary>
        /// <param name="a">Polynomial on the left side</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>true if equal</returns>
        public static bool operator ==(Polynomial a, Polynomial b)
        {
            if ((Object)a == null || (Object)b == null)
            {
                if ((Object)a == null && (Object)b == null)
                    return true;
                else
                    return false;
            }
            return ( ( Enumerable.SequenceEqual<double>(a._coeffs, b._coeffs) ) && (a.isConst() == b.isConst()) && (a._order == b._order));
        }

        /// <summary>
        /// Implements the != operator.<para/>
        /// Compares two Taylor polynomials.<para/>
        /// </summary>
        /// <param name="a">Polynomial on the left side</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>true if not equal</returns>
        public static bool operator !=(Polynomial a, Polynomial b) { return !(a == b); }

        /// <summary>
        /// Overrides the base function Equals
        /// </summary>
        /// <param name="obj">object of Polynomial</param>
        /// <returns>true if coeffs[] are equal</returns>
        public override bool Equals(object obj)
        {
            // If parameter cannot be cast to ThreeDPoint return false:
            Polynomial p = obj as Polynomial;
            if ((object)p == null)
            {
                return false;
            }

            return Enumerable.SequenceEqual<double>(_coeffs, p._coeffs);
        }

        /// <summary>
        /// Overrides the base function GetHashCode()
        /// </summary>
        /// <returns>hash of coeffs[]</returns>
        public override int GetHashCode()
        {
            return _coeffs.GetHashCode();
        }

        /// <summary>
        /// Implements the &lt; operator.<para/>
        /// Compares two Taylor polynomials according to the value of the first coefficient.
        /// </summary>
        /// <param name="a">Polynomial on the left side</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>true if a &lt; b</returns>
        public static bool operator <(Polynomial a, Polynomial b)
        {
            if (a._order < b._order) 
                return true;
            else if (a._order > b._order)
                return false;
            else {
                for (int i = 0; i < a._order; i++)
                {
                    if (a[i] < b[i])
                        return true;
                }
                return false;
            }
        }

        /// <summary>
        /// Implements the &lt;= operator.<para/>
        /// Compares two Taylor polynomials according to the value of the first coefficient.
        /// </summary>
        /// <param name="a">Polynomial on the left side</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>true if a &lt;= b</returns>
        public static bool operator <=(Polynomial a, Polynomial b)
        {
            if (a._order < b._order)
                return true;
            else if (a._order > b._order)
                return false;
            else
            {
                for (int i = 0; i < a._order; i++)
                {
                    if (a[i] > b[i])
                        return false;
                }
                return true;
            }
        }

        /// <summary>
        /// Implements the > operator.<para/>
        /// Compares two Taylor polynomials according to the value of the first coefficient.
        /// </summary>
        /// <param name="a">Polynomial on the left side</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>true if a > b</returns>
        public static bool operator >(Polynomial a, Polynomial b)
        {
            if (a._order < b._order)
                return false;
            else if (a._order > b._order)
                return true;
            else
            {
                for (int i = 0; i < a._order; i++)
                {
                    if (a[i] > b[i])
                        return true;
                }
                return false;
            }
        }

        /// <summary>
        /// Implements the >= operator.<para/>
        /// Compares two Taylor polynomials according to the value of the first coefficient.
        /// </summary>
        /// <param name="a">Polynomial on the left side</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>true if a >= b</returns>
        public static bool operator >=(Polynomial a, Polynomial b)
        {
            if (a._order < b._order)
                return false;
            else if (a._order > b._order)
                return true;
            else
            {
                for (int i = 0; i < a._order; i++)
                {
                    if (a[i] < b[i])
                        return false;
                }
                return true;
            }
        }

        /// <summary>
        /// Implements the + operator
        /// </summary>
        /// <param name="a">Polynomial on the left side</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>a + b (if order mathches)</returns>
        public static Polynomial operator +(Polynomial a, Polynomial b)
        {
            if (a._order != b._order)
            {
                throw new MathException("The order of both Taylor Polynoms should match. ");
            }

            Polynomial retval = new Polynomial(a._order);
            if (a.isConst())
            {
                retval = new Polynomial(b);
                retval[0] += a[0];
            }
            else if (b.isConst())
            {
                retval = new Polynomial(a);
                retval[0] += b[0];
            }
            // general case
            else
            {
                for (int i = 0; i <= a._order; i++)
                    retval[i] = a[i] + b[i];
            }
            return retval;
        }

        /// <summary>
        /// Implements the unary - Operator
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static Polynomial operator -(Polynomial a)
        {
            Polynomial p = new Polynomial(a._order);
            for (int i = 0; i <= a._order; i++)
            {
	            p[i] = - a[i];
            }
            return p;
        }

        /// <summary>
        /// Implements the -= operator.<para/>
        /// Substracts a Taylor polynomial from the current one using pointers to arrays that <para/>
        /// store the coefficients.<para/>
        /// </summary>
        /// <param name="a">the Taylor polynomial to be substracted from.</param>
        /// <param name="b">Polynomial on the right side</param>
        /// <returns>a - b (if order mathches)</returns>
        public static Polynomial operator -(Polynomial a, Polynomial b)
        {
            if (a._order != b._order)
            {
                throw new MathException("The order of both Taylor Polynoms should match. ");
            }

            Polynomial retval = new Polynomial(a._order);
            if (a.isConst())
            {
                retval = new Polynomial(-b);
                retval[0] = a[0] - b[0];
            }
            else if (b.isConst())
            {
                retval = new Polynomial(a);
                retval[0] -= b[0];
            }
            // general case
            else
            {
                for (int i = 0; i <= a._order; i++)
                    retval[i] = a[i] - b[i];
            }
            return retval;
        }

        /// <summary>
        /// Multiplies two Taylor polynomials. Implements the * operator for Taylor arithmetic.<para/>
        /// The following coefficient propagation rule is applied:<para/>
        /// <para/>
        /// 		v_k = sum_{j=0}{k}  u_j * w_{k-j}<para/>
        /// <para/>
        /// for k = 1...d and v(t) = u(t) * w(t), u, v, w being Taylor polynomials, and d being <para/>
        /// the derivative degree.<para/>
        /// <para/>
        /// It is assumed that all three Taylor polynomials have the same derivative degree d.<para/>
        /// Three different cases are distinguished here: when at least one of the polynomials<para/>
        /// is a constant polynomial and when both polynomials are not.<para/>
        /// <para/>
        /// (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques<para/>
        /// of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM,<para/>
        /// Philadelphia, PA, 2000)<para/>
        /// </summary>
        /// <param name="a">the Taylor polynomial to be multiplied with.</param>
        /// <param name="b">the Taylor polynomial to be multiplied by.</param>
        /// <returns>The resulting Polynomail (without changing the order)</returns>
        public static Polynomial operator *(Polynomial a, Polynomial b)
        {
            if (a._order != b._order)
            {
                throw new MathException("The order of both Taylor Polynoms should match.");
            }

            Polynomial retVal = new Polynomial(a._order);

            if (a.isConst())
            {
                retVal = new Polynomial(b);
                for (int i = 0; i <= a._order; i++)
                {
                    retVal[i] *= a[0];
                }
            }
            else if (b.isConst())
            {
                retVal = new Polynomial(a);
                for (int i = 0; i <= b._order; i++)
                {
                    retVal[i] *= b[0];
                }
            }
            else
            {
                // 2. for (int i = 0; i < _order * 2; i++)
                for (int i = 0; i <= a._order; i++)
                {
                    retVal[i] = 0.0;

                    for (int j = 0; j < i + 1; j++)
                    {
                        retVal[i] += a[j] * b[i - j];
                    }

                }
                // 3. v._coeffs[_order*2] = _coeffs[_order] * p._coeffs[_order];
            }
            return retVal;
        }

        /// <summary>
        /// Implements the * operator.<para/>
        /// Multiplies a Taylor polynomial by a scalar.<para/>
        /// </summary>
        /// <param name="a">The Polynomial value to multiply with</param>
        /// <param name="d">The scalar value to multiply by</param>
        /// <returns></returns>
        public static Polynomial operator *(Polynomial a, double d)
        {
            Polynomial retVal = new Polynomial(a);
            for (int i = 0; i <= a._order; i++)
            {
                retVal[i] *= d;
            }

            return retVal;
        }

        /// <summary>
        /// Divides a Taylor polynomial (dividend) by another Taylor polynomial (divisor).<para/>
        /// Implements the / operator for Taylor arithmetic.<para/>
        /// <para/>
        /// 		v_k = 1 / w_0 * [u_k - sum_{j=0}{k-1}  v_j * w_{k-j}]<para/>
        /// <para/>
        /// for k = 1...d and v(t) = u(t) / w(t), u, v, w being Taylor polynomials, and d being <para/>
        /// the derivative degree.<para/>
        /// <para/>
        /// It is assumed that all three Taylor polynomials have the same derivative degree d.<para/>
        /// <para/>
        /// (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques<para/>
        /// of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, <para/>
        /// Philadelphia, PA, 2000)<para/>
        /// </summary>
        /// <param name="a">Polynomial as divisor</param>
        /// <param name="b">Polynomial as divident</param>
        /// <returns>the resulting Taylor polynomial.</returns>
        public static Polynomial operator /(Polynomial a, Polynomial b)
        {

            if (a._order != b._order)
	        {
		        throw new MathException("The order of both Taylor Polynoms should match.");
	        }

	        Polynomial retval = new Polynomial(a._order);

	        for (int i = 0; i <= a._order; i++)
	        {
		        double sum = 0.0;
		        for (int j = 0; j < i; j++)
		        {
			        sum += retval[j] * b[i - j];
		        }
		        retval[i] = (a[i] - sum) / b[0];
	        }
	
	        return retval;
        }

        #endregion

        #region Functions

        /// <summary>
        /// Polynomial PExpect = new Polynomial(3, new double[] { 1, 4, 10, 20 });<para/>
        /// Taylor arithmetic.<para/>
        /// <para/>
        /// The following coefficient propagation rule is applied:<para/>
        /// <para/>
        ///     v_k = sum_{j=0}{k}  u_j * u_{k-j}<para/>
        /// <para/>
        /// for k = 1...d and v(t) = u(t)^2, u and v being a Taylor polynomials, and d being<para/>
        /// the derivative degree.<para/>
        /// <para/>
        /// (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques<para/>
        /// of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, <para/>
        /// Philadelphia, PA, 2000)<para/>
        /// </summary>
        /// <returns>The resulting Taylor polynomial.</returns>
        public Polynomial sqr()
        {
            Polynomial retval = new Polynomial(_order);
            for (int i = 0; i <= _order; i++)
            {
                retval[i] = 0.0;
                for (int j = 0; j < i + 1; j++)
                    retval[i] += this[j] * this[i - j];
            }
            return retval;
        }

        /// <summary>
        /// Sets this Taylor polynomial to its square.<para/>
        /// <para/>
        /// The following coefficient propagation rule is applied:<para/>
        /// <para/>
        /// 		v_k = sum_{j=0}{k}  u_j * u_{k-j}<para/>
        /// <para/>
        /// for k = 1...d and v(t) = u(t)^2, u and v being a Taylor polynomials, and d being<para/> 
        /// the derivative degree.<para/>
        /// <para/>
        /// (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques<para/>
        /// of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, <para/>
        /// Philadelphia, PA, 2000)<para/>
        /// </summary>
        public void setSqr()
        {
            Polynomial temp = this.sqr();
            this.initializePolynomial(order, temp._coeffs, temp._constant);
        }

        /// <summary>
        /// Calculates the square root of a Taylor polynomial. Implements the square root function <para/>
        /// for Taylor arithmetic.<para/>
        /// <para/>
        /// The following coefficient propagation rule is applied:<para/>
        /// <para/>
        /// 		v_k = 1 / 2*v_0 * [u_k - sum_{j=1}{k-1}  v_j * v_{k-j}]<para/>
        /// <para/>
        /// for k = 1...d and v(t) = sqrt(u(t)), u and v being a Taylor polynomials, and d being <para/>
        /// the derivative degree. In particular, v_0 = sqrt(u_0).<para/>
        /// <para/>
        /// (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques<para/>
        /// of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, <para/>
        /// Philadelphia, PA, 2000)<para/>
        /// </summary>
        /// <returns></returns>
        public Polynomial sqrt()
        {
            Polynomial retval = new Polynomial(_order);

            retval[0] = Math.Sqrt(this[0]);

            for (int i = 1; i <= _order; i++)
            {
                double sum = 0.0;
                for (int j = 1; j < i; j++)
                    sum += retval[j] * retval[i - j];

                retval[i] = ( this[i] - sum ) / (2 * retval[0]);
            }

            return retval;
        }

        /// <summary>
        /// Sets this Taylor polynomial to its square root.<para/>
        /// <para/>
        /// The following coefficient propagation rule is applied:<para/>
        /// <para/>
        /// 		v_k = 1 / 2*v_0 * [u_k - sum_{j=1}{k-1}  v_j * v_{k-j}]<para/>
        /// <para/>
        /// for k = 1...d and v(t) = sqrt(u(t)), u and v being a Taylor polynomials, and d being <para/>
        /// the derivative degree. In particular, v_0 = sqrt(u_0).<para/>
        /// <para/>
        /// (See Griewank's book, p.222 from "Evaluating Derivatives: Principles and Techniques<para/>
        /// of Algorithmic Differentiation". In Frontiers in Applied Mathematics Nr. 19, SIAM, <para/>
        /// Philadelphia, PA, 2000)<para/>
        /// </summary>
        public void setSqrt()
        {
            Polynomial temp = this.sqrt();
            this.initializePolynomial(order, temp._coeffs, temp._constant);
        }

        /// <summary>
        /// Prints out the coefficients of a Taylor polynomial, starting by the independent term.<para/>
        /// Prints out to Standard Output (Console)<para/>
        /// </summary>
        [Obsolete("Use ToString() and print where you want.")] 
        public void print()
        {
            String output = String.Empty;
            for (int i = 0; i <= _order; i++)
                output += this[i] + "\t";

            Console.WriteLine(output);
        }

        /// <summary>
        /// Prints out to a file the coefficients of a Taylor polynomial, starting by the <para/>
        /// independent term.<para/>
        /// </summary>
        /// <param name="filenameWithPath">The output file to write the polynomial to.</param>
        public void print(String filenameWithPath)
        {
            String output = String.Empty;
            for (int i = 0; i <= _order; i++)
                output += this[i] + "\t";

            File.WriteAllText(filenameWithPath, output + System.Environment.NewLine);
        }

        /// <summary>
        /// Evaluates a Taylor polynomial at a given value with a point of expansion.<para/>
        /// </summary>
        /// <param name="x">The value to evaluate the polynomial at.</param>
        /// <param name="alpha">alpha The point of expansion.</param>
        /// <returns>The result of the evaluation.</returns>
        public double eval(double x, double alpha)
        {
            double result = this[_order];
            double t = x - alpha;

            for (int i = _order - 1; i >= 0; i--)
            {
                result = t * result + this[i];
            }

            return result;
        }

        /// <summary>
        /// Returns the first coefficient of the Taylor polynomial, i.e., the evaluation of the function.
        /// </summary>
        /// <returns>The evaluation of the function at the initial point.</returns>
        public double feval()
        {
	        return this[0];
        }

        /// <summary>
        /// Implements the SHIFT operator to calculate the derivative of a Taylor polynomial.<para/>
        /// <para/>
        /// The new coefficients are shifted to the left and the last one is zeroed.<para/>
        /// E.g.:<para/>
        /// 		y(t) = sum_{j=0}^{d} y_j * t^j + O(t^d+1)<para/>
        /// 			 = y_0 + y_1*t + y_2*t^2 + ... + y_d*t^d<para/>
        /// <para/>
        /// 		y'(t) = y_1 + 2*y_2*t + 3*y_3*t^2 + ... + d*y_d*t^d-1 + 0<para/>
        /// <para/>
        /// </summary>
        public void shift()
        {
	        for (int i = 1; i <= _order; i++)
	        {
		        this[i - 1] = i * this[i]; 
	        }

	        // set last coefficient to zero
	        this[_order] = 0.0;
        }

        /// <summary>
        /// Checks if Polynomial is constant<para/>
        /// (if _constant unknown then it will be set)<para/>
        /// </summary>
        public bool isConst()
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

        /// <summary>
        /// Returns true in case it is near a constant Taylor polynomial; <para/>
        /// false otherwise.<para/>
        /// </summary>
        /// <param name="eps">The threshold value to compare with.</param>
        /// <returns>\a true if it is a constant Taylor polynomial; \a false otherwise.</returns>
        public bool isConst(double eps)
        {
            for (int i = 1; i <= _order; i++)
            {
                if (Math.Abs(this[i]) > eps)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Returns \a true in case it is a constant Taylor polynomial with value 1; <para/>
        /// \a false otherwise.<para/>
        /// </summary>
        /// <returns>\a true if it is a constant Taylor polynomial with value 1; \a false otherwise.</returns>
        public bool isId()
        {
            return (isConst() && this[0] == (double)1);
        }

        /// <summary>
        /// Returns \a true in case it is near a constant Taylor polynomial with value 1; <para/>
        /// \a false otherwise.<para/>
        /// </summary>
        /// <param name="eps">The threshold value to compare with.</param>
        /// <returns>\a true if it is a constant Taylor polynomial with value 1; \a false otherwise.</returns>
        public bool isId(double eps)
        {
            return (isConst(eps) && Math.Abs(this[0] - (double)1) < eps);
        }

        /// <summary>
        /// Returns \a true in case all coefficients of the Taylor polynomial are zeroed; <para/>
        /// \a false otherwise.<para/>
        /// </summary>
        /// <returns>\a true if all coefficients are zeroed; \a false otherwise.</returns>
        public bool isZero()
        {
            return this[0] == 0 && isConst();
        }

        /// <summary>
        ///  Returns \a true in case all coefficients of the Taylor polynomial are lower or equal than<para/>
        ///  a threshold given as parameter; \a false otherwise.<para/>
        /// </summary>
        /// <param name="eps">The threshold value to compare with.</param>
        /// <returns>\a true if all coefficients are almost null; \a false otherwise.</returns>
        public bool isZero(double eps)
        {
            for (int i = 0; i <= _order; i++)
            {
                if (Math.Abs(this[i]) > eps)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Sets all coefficients of a Taylor polynomial to zero.<para/>
        /// </summary>
        public void set2Zero()
        {
            for (int i = 0; i <= _order; i++)
            {
                this[i] = 0.0;
            }
            _constant = 1;
        }

        /// <summary>
        /// Sets the coefficients of a Taylor polynomial to zero, from the order given as parameter on.
        /// </summary>
        /// <param name="order">Derivative order from which to start on (increasingly).</param>
        public void set2Zero(int order)
        {
            for (int i = order; i <= _order; i++)
            {
                this[i] = 0.0;
            }
            _constant = -1;
        }

        /// <summary>
        /// Sets a Taylor polynomial to the constant given as parameter.
        /// </summary>
        /// <param name="c">The constant value of type \a double to set the Taylor polynomial to.</param>
        public void set2const(double c)
        {
            this[0] = c;
            for (int i = 1; i <= _order; i++)
                this[i] = 0;

            _constant = 1;
        }

        /// <summary>
        /// Sets the coefficients of a Taylor polynomial to the ones given as parameter.
        /// </summary>
        /// <param name="c">A vector of coefficients of type \type double.</param>
        public void setCoeffs(double[] c)
        {
            if (c.Length - 1 > _order)
                throw new MathException("Wrong Dimension of values...");

            for (int i = 0; i <= _order; i++)
            {
                if (c.Length - 1 < i)
                    this[i] = 0;
                else
                    this[i] = c[i];
            }
        }

        /// <summary>
        /// sets the Polynomial to id
        /// </summary>
        public void set2Id() { 
            set2const(1.0); 
        }

        /// <summary>
        /// Returns value at the given index from the array
        /// </summary>
        /// <param name="index">The index to be analyzed.</param>
        /// <returns>The coefficient at that index.</returns>
        public double getValueAt(int index)
        {
            if (index > _order || index < 0)
                throw new MathException("Invalid coefficient index of polynomial.");
            else
                return _coeffs[index];
        }

        /// <summary>
        /// Returns a String of the Polynomial
        /// </summary>
        /// <returns>String in the Format 2x^2 + -1x^1 + 7</returns>
        override public String ToString()
        {
            String retval = String.Empty;
            for (int i = _order; i > 0; i--)
            {
                if (i != _order)
                    retval += " + ";

                retval += this[i] + "x^" + i + "\t";
            }
            if (_order > 0)
                retval += " + " + this[0];
            else
                retval += this[0];
            return retval;
        }
        
        #endregion

        #endregion
    }
}
