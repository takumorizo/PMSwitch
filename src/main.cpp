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
	    if(not param.isDBUpdate){
			InferenceData<long long, long double> maxAns = inference.vbFull(false, 1e-4, param.pathOutputDir+"/err");
			for(int i = 0; i < 10-1; i++){
			    InferenceData<long long, long double> ans = inference.vbFull(false, 1e-4, param.pathOutputDir+"/err");
			    if(maxAns.getLq() < ans.getLq()) maxAns = ans;
			}
		    std::cerr << "end VB" << std::endl;
		 	maxAns.printHyperParameters(param.pathOutputDir);
		 	maxAns.printLatents(param.pathOutputDir);
		 	maxAns.printDB(param.pathOutputDir);

	    }else{
			InferenceData<long long, long double> maxAns = inference.vb(true, 1e-4, param.pathOutputDir+"/err");
			for(int i = 0; i < 10-1; i++){
			    InferenceData<long long, long double> ans = inference.vb(param.isDBUpdate, 1e-4, param.pathOutputDir+"/err");
			    if(maxAns.getLq() < ans.getLq()) maxAns = ans;
			}
		    std::cerr << "end VB" << std::endl;
		 	maxAns.printHyperParameters(param.pathOutputDir);
		 	maxAns.printLatents(param.pathOutputDir);
		 	maxAns.printDB(param.pathOutputDir);
	    }
	}
}