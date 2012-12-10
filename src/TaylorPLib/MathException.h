#ifndef _MATH_EXCEPTION_H_
#define _MATH_EXCEPTION_H_

#include <stdio.h>
#include <stdarg.h>

#define MATH_ERROR(format, ...) MathException("%s(%d): " format, __FILE__, __LINE__, __VA_ARGS__)

#define MAX_MESSAGE_SIZE 256

class MathException
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

#endif