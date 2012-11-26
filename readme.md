TaylorPLib is a library for calculations on matrices of Taylor Polynoms
=======================================================================

The project just started, but we plan to release the C++ version by the end of February.


Dependencies
------------
- Google C++ Testing Framework Versien 1.6.0 (newer version should work too, but weren't tested)

Compiling the Library
---------------------

### Windows MSVC 2010
1. Clone the repository to some folder F  
   *F shoud now contain the folders src/, src-old/, test/*
2. Download and extract GTest to F\GTest\
3. Open the F\GTest\msvc\gtest.sln in Visual Studio and compile it  
   *This should produce a F\GTest\msvc\gtest\Debug\ with the binaries (we need the lib files)*
4. Close the GTest project and open F\TaylorPLib\src\TaylorPLib.sln
5. Compile and test the library with Strg+F5

Using the Library
-----------------

Coming soon