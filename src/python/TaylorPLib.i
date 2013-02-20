/* File : TaylorPLib.i */
%module TaylorPLib

%{
#include "..\src\TaylorPLib.h"
%}

%include "std_vector.i"

namespace std {
	%template(DoubleArray) vector<double>;
	%template(IntArray) vector<int>;
    %template(2dArray) vector< vector<double> >;
}

/* Let's just grab the original header file here */
%include "..\src\TaylorPLib.h"
