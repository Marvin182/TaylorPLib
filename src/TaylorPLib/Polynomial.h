#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <stdexcept>
#include <math.h>
#include "IDException.h"

namespace LibMatrix {

	class __declspec(dllexport) class Polynomial
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
		int itsOrder;				// The derivative order of the Taylor polynomial
									// Nr. of coeff = itsGrade + 1
    	double *itsCoeff;			// The coefficients
    	bool _constant;				// whether it is a constant Taylor polynomial or not

	public:
		//
		// Constructors, destructor
		//
		Polynomial();
		Polynomial(int order);
		Polynomial(const Polynomial &p);
		~Polynomial();

		//
		// Accessing properties
		//
		double* coeffs() { return itsCoeff; }
		int order() { return itsOrder; }
		int nrcoeff() { return itsOrder + 1; }
		char* typeName() { return "Polynomial"; }

		//
		// Overloaded operators for Taylor arithmetic
		//
		double & operator[](int index);
		Polynomial operator=(const Polynomial &p);
		bool operator==(const Polynomial &p) const;
		bool operator!=(const Polynomial &p) const;
		int operator<(const Polynomial &p);
		int operator<=(const Polynomial &p);
		int operator>(const Polynomial &p);
		int operator>=(const Polynomial &p);
		Polynomial operator+(const Polynomial &w);
		Polynomial operator+=(const Polynomial &w);
		Polynomial operator-(const Polynomial &w);
		Polynomial operator-=(const Polynomial &w);
		Polynomial operator-();
		Polynomial operator*(const Polynomial &w);
		Polynomial operator*=(const Polynomial &w);
		Polynomial operator*(double alpha);
		Polynomial operator/(const Polynomial &w);
		Polynomial operator/=(const Polynomial &w);

		//
		// Special operators
		//
		Polynomial sqr();
		Polynomial sqrt();

		//
		// Other functions
		//
		void print();
		void print(FILE * fn);
		double eval(double x, double alpha);
		double feval();
		void shift();
		bool isConst();
		bool isConst(double eps);
		bool isId();
		bool isId(double eps);
		bool isZero();
		bool isZero(double eps);
		int set2zero();
		int set2zero(int ord);
		int set2const(double c);
		int setCoeffs(double *c);
	};
};

#endif /*POLYNOMIAL_H_*/
