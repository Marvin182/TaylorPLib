#ifndef _MATH_EXCEPTION_H_
#define _MATH_EXCEPTION_H_

#include <stdio.h>
#include <stdarg.h>

// this should also work, but Visual Studio doen't like the ... in tha macro parameter list
// #define MATH_ERROR(format,x...) MathException("%s(%d): " format, __FILE__, __LINE__, x)

class MathException
{
	private:
		char* msg;
	public:
		MathException(const char* format, ...): 
			msg(0)
		{
			int msgSize = 256;
			msg = (char*) malloc(msgSize);
			
			va_list args;
			va_start(args, format);
			vsnprintf_s(msg, msgSize, _TRUNCATE, format, args);
			va_end(args);
		}
		// ~MathException() { free(msg); }

		const char* what() const { return msg; }
};

#endif