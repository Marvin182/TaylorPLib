TaylorPLib is a library for calculations on matrices from Taylor Polynoms
=========================================================================

This projects aims to create a library to help researches working with Taylor Polynoms and big matrices from Taylor Polynoms. It contains many methods for faster multiplications of matrices with a special form (e.g. upper left matrix, column pivoting). The main library it written in C++ but it is also avaible as a Python module (big thanks to the SWIG project) and as a C# port.

The plan to release a stable version by March 1st.


C++
---
### Dependencies
- Google C++ Testing Framework version 1.6.0 or newer

### Compiling the libary using Windows MSVC 2010
1. Clone the repository to some folder F  
   *F shoud now contain the folder src/
2. Download and extract GTest to F\GTest\
3. Open the F\GTest\msvc\gtest.sln in Visual Studio and compile it  
   *This should produce a F\GTest\msvc\gtest\Debug\ with the binaries (we need the lib files)*
4. Close the GTest project and open F\TaylorPLib\src\TaylorPLib.sln
5. Compile and test the library with Strg + F5


Python
------

### Dependencies
- SWIG (we used version 2.0.9 but others versions might work, too)


C#
--
