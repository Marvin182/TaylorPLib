#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <stdexcept>
#include <math.h>
#include "CustomException.h"

#define abs(x) ((x) > 0 ? (x) : -(x))

namespace LibMatrix {

	class __declspec(dllexport) Polynomial
	{
	//
	// (Taylor) Polynomial with derivate degree n (n+1 coefficients):
	//
	// P_n(x) = f(a) + (x-a)*f'(a)/1! + (x-a)^2*f''(a)/2! + (x-a)^3*f'''(a)/3! 
	//			+ ... + (x-a)^n*f^{(n)}(a)/n!
	//
	//		  = sum{k=0}{n} (x-a)^k*f^{(k)}(a)/k!
	//
	// being 'f' the function to be approximated by P_n(x) at point 'a', with its
	// first n derivatives existing on a closed interval I, so that
	//
	// f(x) = P_n(x) + R_n(x)
	//
	// the remainder term being R_n(x) = (x-a)^{n+1}*f^{(n+1)}(c)/(n+1)! for some
	// 'c' between 'x' and 'a'.
	//
	// Another often used form:
	//
	// f(x0+h) = f(x0) + h*f'(x0)/1! + h^2*f''(x0)/2! + h^3*f'''(x0)/3! 
	//			+ ... + h^n*f^{(n)}(x0)/n!
	//
	private:
    	mutable char _constant;		// whether it is a constant (Taylor) polynomial or not
    								//  0 = not constant
    								//  1 = constant
    								// -1 = unkown
		int _order;					// the derivative order of the (Taylor) polynomial
    	double *_coeffs;			// coefficients, number of coeffientk = order + 1

    	void allocateMemory(bool initialize) { allocateMemory(_coeffs, _order, initialize); }
    	void allocateMemory(double *&coeffs, int order, bool initialize);
    	void deallocateMemory() { deallocateMemory(_coeffs); }
    	void deallocateMemory(double *&coeffs);
    	void copyFrom(const Polynomial &p);

    	static int unsetConstCount;
    	static int isConstCount;
    	void unsetConst();

	public:
		//
		// Constructors, destructor
		//
		Polynomial();
		Polynomial(int order, bool initialize = true);
		Polynomial(const Polynomial &p);
		~Polynomial();

		//
		// Accessing properties
		//
		int order() const { return _order; }
		int ncoeff() const { return _order + 1; }

		//
		// Overloaded operators for Taylor arithmetic
		//
		double & operator[](int index) const;
		Polynomial operator=(const Polynomial &p);
		bool operator==(const Polynomial &p) const;
		bool operator!=(const Polynomial &p) const;
		bool operator<(const Polynomial &p);
		bool operator<=(const Polynomial &p);
		bool operator>(const Polynomial &p);
		bool operator>=(const Polynomial &p);
		
		Polynomial operator+(const Polynomial &p) const;
		Polynomial operator+=(const Polynomial &p);
		
		Polynomial operator-() const;
		Polynomial operator-(const Polynomial &p) const;
		Polynomial operator-=(const Polynomial &p);
		
		Polynomial operator*(const double d) const;
		Polynomial operator*(const Polynomial &p) const;
		Polynomial operator*=(const double d);
		Polynomial operator*=(const Polynomial &p);

		Polynomial operator/(const Polynomial &p) const;
		Polynomial operator/=(const Polynomial &p);

		//
		// Special operators
		//
		Polynomial sqr() const;
		void setSqr();
		Polynomial sqrt() const;
		void setSqrt();

		//
		// Other functions
		//
		void print();
		void print(FILE * fn);
		double eval(double x, double alpha);
		double feval();
		void shift();
		bool isConst() const;
		bool isConst(double eps);
		bool isId() const;
		bool isId(double eps);
		bool isZero() const;
		bool isZero(double eps);
		void set2zero();
		void set2zero(int ord);
		void set2const(double c);
		void setCoeffs(double *c);
	};
};

#endif /*POLYNOMIAL_H_*/
