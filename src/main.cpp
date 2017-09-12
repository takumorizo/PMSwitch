#include "cmdline.h"
#include "Parameters.hpp"
#include "D2Vector.hpp"
#include "FixedSizeMultiVector.hpp"
#include "ExceptionUtil.hpp"
#include "InputFileParser.hpp"
#include "Inference.hpp"
#include "FeatureData.hpp"
#include "DBData.hpp"

#include <vector>
#include <iostream>
#include <assert.h>

int main(int argc, const char *argv[]){
	using namespace pmswitch;

	Parameters<long long, long double> param(argc, argv);
	param.getFromCommandLineArguments(argc, argv);

	if(param.method == Parameters<long long, long double>::VB){
	    Inference<long long, long double> inference = InferenceCreator<long long, long double>::createInference(param);

	    InferenceData<long long, long double> ans = inference.vb(param.isDBUpdate, 1e-4);
	    std::cerr << "end VB" << std::endl;
	 	ans.printHyperParameters(param.pathOutputDir);
	 	ans.printDB(param.pathOutputDir);
	}
}