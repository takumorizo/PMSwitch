// #include "sim_test.h"
// #include "../tests/src/sim_test.h"
#include "cmdline.h"
#include "Parameters.hpp"
#include "D2Vector.hpp"
#include "FixedSizeMultiVector.hpp"
#include "ExceptionUtil.hpp"
#include "InputFileParser.hpp"

#include <vector>
#include <iostream>
#include <assert.h>

int main(int argc, const char *argv[]){
	using namespace pmswitch;
	Parameters<int> param(argc, argv);
	param.getFromCommandLineArguments(argc, argv);
	
}
	// // InputFileParser parser;
	// D2Vector<int,int> vec(3,6,3);
	// FixedSizeMultiVector<int, int> vec2(3,3,2);
	// FixedSizeMultiVector<int, int> vecvec;
	// vecvec = vec2;
	// std::cout << "vecvec" << std::endl;
	// vecvec.print();
	// std::vector<int> v(6,3);
	// // for(int i = 0; i < 6; i++){ vec[i] = i;}
	// // for(int i = 0; i < 6; i++){ std::cout << vec[i] << std::endl; }

	// // std::cout << " === " << std::endl;
	// // for(int i = 0; i < 2; i++)for(int j = 0; j < 3 ; j++){ std::cout << vec(i,j) << std::endl;}

	// // long long sum = 0;
	// // for(long long i = 0; i < 1e9; i++){
	// // 	for(int i = 0; i < 2; i++)for(int j = 0; j < 3 ; j++){
	// // 		sum += vec(i,j);
	// // 	}
	// // }
	// // std::cout << "D2Vector" << std::endl;
	// // std::cout << sum << std::endl;

	// long long sum = 0;
	// for(long long i = 0; i < 1e9; i++){
	// 	for(int i = 0; i < 3; i++)for(int j = 0; j < 2 ; j++){
	// 		sum += vec2(i,j);
	// 	}
	// }

	// std::cout << "FixedSizeMultiVector" << std::endl;
	// std::cout << sum << std::endl;

	// // long long sum = 0;
	// // for(long long i = 0; i < 1e9; i++){
	// // 	for(int i = 0; i < 2; i++)for(int j = 0; j < 3 ; j++){
	// // 		sum += v[2*i+j];
	// // 	}
	// // }
	// // std::cout << "std::vector" << std::endl;
	// // std::cout << sum << std::endl;

