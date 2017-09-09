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
	Parameters<long long, double> param(argc, argv);
	param.getFromCommandLineArguments(argc, argv);

	if(param.method == Parameters<int>::VB){
	    Inference<long long, double> inference = InferenceCreator<long long, double>::createInference(param);
	    InferenceData<long long, double> ans0 = inference.vb(true, 1e-4);
	    InferenceData<long long, double> ans1 = inference.vb(true, 1e-5);
	    ans0 = ans1;
	}

}