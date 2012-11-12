#ifndef _CUSTOMEXCEPTION_H
#define _CUSTOMEXCEPTION_H

#include <string>

class __declspec(dllexport) CustomException
{
public:
	CustomException();
	CustomException( char* msg );
	CustomException( char* msg, int err );
	// ~CustomException();
	void addReason( char* msg );
	void addReason( char* msg, int err );
	std::string what();
	int getErrCode();
private:
	char* reason;
	int errCode;
};

#endif // _CUSTOMEXCEPTION_H
