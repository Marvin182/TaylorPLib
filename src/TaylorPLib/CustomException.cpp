#include "CustomException.h"

#include <sstream>
#include <iostream>

using namespace std;

CustomException::CustomException() : reason(""), errCode(0)
{
}

CustomException::CustomException( char* msg ) : reason(msg), errCode(0)
{
}

CustomException::CustomException( char* msg, int err ) : reason(msg), errCode(err)
{
}

// CustomException::~CustomException()
// {
// }

void CustomException::addReason( char* msg )
{
	this->reason = msg;
}

void CustomException::addReason( char* msg, int err )
{
	this->reason = msg;
	this->errCode = err;
}

string CustomException::what()
{
	stringstream s;
	s << "***The exception number " << errCode << " occurred.\n";
	s << "***Reason: " << reason;
	return s.str();
}

int CustomException::getErrCode()
{
	return errCode;
}
