#ifndef ExceptionUtil_h
#define ExceptionUtil_h

#include <string>
#include <iostream>

namespace pmswitch{
	void die(std::string message, int retCode = 1);
}

void pmswitch::die(std::string message, int retCode){
		std::cerr << message << std::endl;
		exit(retCode);
}

#endif