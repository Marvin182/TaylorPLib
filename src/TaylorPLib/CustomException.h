#ifndef _CUSTOMEXCEPTION_H
#define _CUSTOMEXCEPTION_H

#include <string>

class __declspec(dllexport) CustomException
{
public:
	CustomException();
	CustomException( const char *msg );
	CustomException( const char *msg, int err );
	~CustomException();
	void addReason( const char *msg );
	void addReason( const char *msg, int err );
	std::string what();
	int getErrCode();
private:
	const char* reason;
	int errCode;
};

#endif // _CUSTOMEXCEPTION_H
