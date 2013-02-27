/* File : TaylorPLib.i */
%module TaylorPLib

%{
#include "..\src\TaylorPLib.h"
%}

/* Support for vector<T> */

%include "std_vector.i"

namespace std {
	%template(DoubleArray) vector<double>;
	%template(IntArray) vector<int>;
    %template(PolynomialArray) vector<LibMatrix::Polynomial>;
}

%ignore LibMatrix::Polynomial::operator=(const Polynomial &p);

%rename (_print) print();
%rename (printToFile) print(FILE*);
%rename (printWithName) print(const char*);

%ignore LibMatrix::Matrix::operator=(const Matrix &p);
%ignore LibMatrix::Polynomial::operator[](int);

%ignore operator<<(std::ostream &out, const LibMatrix::Polynomial &p);
%ignore operator<<(std::ostream &out, const LibMatrix::Matrix &m);


/* include orignial header file */

%include "..\src\TaylorPLib.h"