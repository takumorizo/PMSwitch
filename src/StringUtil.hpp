#ifndef StringUtil_h
#define StringUtil_h

// #include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

namespace pmswitch{
	namespace string {
		std::vector<std::string> split(const std::string &str, char delimiter);
	}
}

std::vector<std::string> pmswitch::string::split(const std::string &str, char delimiter){
	std::vector<std::string> ans;
	std::istringstream sstr(str);
	std::string token;
	while (getline(sstr, token, delimiter)) ans.push_back(token);
	return ans;
}


#endif