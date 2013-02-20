/* File : TaylorPLib.i */
%module TaylorPLib

%{
#include "..\src\TaylorPLib.h"
%}

%include "std_vector.i"

namespace std {
	%template(Line)  vector <double>;
    %template(Array) vector< vector<double> >;
}

/* Let's just grab the original header file here */
%include "..\src\TaylorPLib.h"
