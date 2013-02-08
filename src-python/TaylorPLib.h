#ifndef _TAYLORPLIB_H_
#define _TAYLORPLIB_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdexcept>
#include <math.h>
#include <typeinfo>

#define abs(x) ((x) > 0 ? (x) : -(x))

#define MAX_MESSAGE_SIZE 256
// #define MATH_ERROR(format, ...) MathException("%s(%d): " format, __FILE__, __LINE__, __VA_ARGS__)

namespace LibMatrix {

#ifdef SWIG
	class Polynomial
#else
	class __declspec(dllexport) Polynomial
#endif
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
		void set2Const(double c);
		void set2Id() { set2Const(1.0); }
		void set2Zero();
		void set2Zero(int ord);
		void setCoeffs(double *c);
	};

#ifdef SWIG
	class Matrix
#else
	class __declspec(dllexport) Matrix
#endif
	{
	private:
		int _rows,											// The number of rows
			_cols,											// The number of columns
			_dimT;											// The dimension of the Taylor polynomials 
    	Polynomial **_data;										// The pointer to the allocated memory

    	void allocateMemory(bool initialize)
    		{ Matrix::allocateMemory(_data, _rows, _cols, _dimT, initialize); }
    	static void allocateMemory(Polynomial **&data, int rows, int cols, int dimT, bool initialize);
	    void deallocateMemory()
	    	{ Matrix::deallocateMemory(_data, _rows, _cols); }
	    static void deallocateMemory(Polynomial **&data, int rows, int cols);
		void copyFrom(const Matrix &m);

	public:
		//
		// Constructors, destructor
		//
		Matrix();														// Default constructor
		Matrix(int rows, int cols, bool initialize = true);				// Regular constructor
		Matrix(int rows, int cols, int dimT, bool initialize = true);	// Regular constructor
		Matrix(const Matrix &m);										// Copy constructor

		Matrix(int rows, int cols, Polynomial *values);						// Test construtor		

		~Matrix();														// Destructor

		//
		// Accessing properties
		//
		// double** data() { return _data; }			// Returns the ptr to the alloc. memory
		int nrows() const { return _rows; }			// Returns the number of rows
		int ncols() const { return _cols; }			// Returns the number of columns
		int dimT() const { return _dimT; }			// Returns the dimension of the type T
		const Polynomial* get(int row, int col) const;			// Returns a single element from the matrix

		//
		// Overloaded operators
		//
		Polynomial & operator()(int row, int col);			// Element
		Matrix operator=(const Matrix &m);				// Assignment operators
		bool operator==(const Matrix &m) const;			// Comparison operators
		bool operator!=(const Matrix &m) const;
		Matrix operator+(const Matrix &m) const;		// Arithmetic operators
		Matrix operator+=(const Matrix &m);
		Matrix operator-(const Matrix &m) const;
		Matrix operator-=(const Matrix &m);
		Matrix operator-();
		Matrix operator*(double alpha) const;
		Matrix operator*=(double alpha);
		Matrix operator*(const Matrix &m) const;		// simple matrix multiplication

		//
		// Especial matrix multiplications
		//
		void mmCaABbC(double alpha, double beta, const Matrix &A, const Matrix &B);
		void bmmCaABbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
		void mmCasABbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
		void mmCaAsBbC(int r, double alpha, double beta, const Matrix&A, const Matrix&B);
		void mmCaAUTBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
		void mmCaAATbC(double alpha, double beta, const Matrix&A);
		void mmCaATAbC(double alpha, double beta, const Matrix&A);
		void mmCaATBbC(double alpha, double beta, const Matrix&A, const Matrix&B);
		void mmCaATBPbC(double alpha, double beta, const Matrix&A, const Matrix&B, int *piv);
		void mmCaABTbC(double alpha, double beta, const Matrix&A, const Matrix&B);
		void mmCaABTbC(int r, bool up, double alpha, double beta, const Matrix&A, const Matrix&B);
		void bmmCaABTbC(int r, int c, double alpha, double beta, const Matrix&A, const Matrix&B);
		void mmCaIBbC(double alpha, double beta, const Matrix&B);
		void mmCaIBbC(double alpha, double beta, int *piv, bool rows, const Matrix&B);
		void mmCaAIbC(double alpha, double beta, const Matrix&A);
		void mmCaAIbC(double alpha, double beta, const Matrix &A, int *piv, bool rows);

		//
		// Solving systems of equations
		//
		// int utsolve(Matrix &B);					// U X = B, back-substitution
		// int utsolve(Matrix &B, Matrix &X, int *piv);// U X = B, back-substitution, pivoting
		// int utsolve(T *b);						// U x = b, back-substitution
		// int utxsolve(Matrix &B);					// X U = B, back-substitution
		// int gsolve(Matrix &B);					// A X = B, Gaussian elimination
		// int gsolve(T *b);							// A x = b, Gaussian elimination

		//
		// Other functions
		//
		// double fnorm();								// Frobenius norm
		// int tcfnorm(int nrTC, double *fn);		// vector of Frobenius norms
		// int colnorm(double *c);
		// int colnormdown(int pos, int *piv, double *c);
		void cpermutem(int *piv, bool trans = false);
		void rpermutem(int *piv);
		void transpose();
		Matrix asTranspose() const;
		void shift();
		bool isId() const;
		// bool isId(double eps);
		// bool isId(int m1, int m2, int n1, int n2, double eps);
		bool isZero() const;
		// bool isZero(double eps);
		void set2Id();
		void set2Id(int top, int bottom, int left, int right);
		void set2IdFromIndices(int firstRow, int lastRow, int firstCol, int lastCol);
		void set2Zero();
		void set2Zero(int top, int bottom, int left, int right);
		void set2ZeroFromIndices(int firstRow, int lastRow, int firstCol, int lastCol);
		void set2Val(double v);
		void set2Val(int top, int bottom, int left, int right, double v);
		void set2ValFromIndices(int firstRow, int lastRow, int firstCol, int lastCol, double v);
		// int trinvm(Matrix &Am);
		// int trinvm(int r, Matrix &Am);
		// bool mcompare(Matrix &B);
		// bool mcompare(Matrix &B, double eps);
		// bool mcompare(Matrix &B, int r, double eps);

		//
		// Print out functions
		//
		void print();
		void print(const char *name);
		// void fprintm(FILE * fn, char *str);
		// void printm(char *str, const char *const color);
		// void printdm(char *str);
		// void fprintdm(FILE * fn, char *str);
		// void printdm(char *str, const char *const color);
		// void printdm(char *str, const char *const color, double eps);
		// void fprintdm(FILE * fn, char *str, double eps);
		// void printm(int *piv, char *str, const char *const color);
		// void printdm(int *piv, char *str, const char *const color, double eps);
		// void fprintm(FILE * fn, int *piv, char *str);
		// void fprintdm(FILE * fn, int *piv, char *str, double eps);
		// void printtpm(char *str, const char *const color, double eps);
		// void fprinttpm(FILE * fn, char *str, double eps);
		// void printtpm(int m1, int m2, int n1, int n2,
		// 		char *str, const char *const color, double eps);
		// void fprinttpm(FILE * fn, int m1, int m2, int n1, int n2,
		// 		char *str, double eps);
		// void printtpm(int *piv, char *str, const char *const color, double eps);
		// void fprinttpm(FILE * fn, int *piv, char *str, double eps);
	};

};

#ifdef SWIG
	class MathException
#else
	class __declspec(dllexport) MathException
#endif
{
private:
	char* _message;
	
public:
	MathException(const char* format, ...): 
		_message(0)
	{
		va_list args;
		va_start(args, format);

		// format message
		_message = (char*) malloc(MAX_MESSAGE_SIZE);
		vsnprintf_s(_message, MAX_MESSAGE_SIZE, _TRUNCATE, format, args);
			
		va_end(args);
	}
	MathException(const MathException &me)
	{
		// copy message
		_message = (char*) malloc(MAX_MESSAGE_SIZE);
		memcpy(_message, me._message, MAX_MESSAGE_SIZE);
	}
	~MathException() { free(_message); }

	const char* what() const { return _message; }
};

#endif // _TAYLORPLIB_H_