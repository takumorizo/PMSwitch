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

	if(param.method == Parameters<long long, double>::VB){
	    Inference<long long, double> inference = InferenceCreator<long long, double>::createInference(param);

		std::cerr << "start vb " << std::endl;
	    InferenceData<long long, double> ans = inference.vb(true, 1e-4);
	 	std::cerr << "end vb " << std::endl;

	 	FixedSizeMultiVector<double, long long> alpha = ans.getAlpha();
	 	FixedSizeMultiVector<double, long long> beta = ans.getBeta();

	 	long long sum = 0.0;
	 	for(long long i = 0; i < alpha.size(); i++){
	 		sum += alpha[i];
	 	}
	 	std::cerr << "alpha sum : " << sum << std::endl;

	 	alpha.print();

	 	std::cerr << "pathOutputDir: " << param.pathOutputDir << std::endl;
	 	ans.printHyperParameters(param.pathOutputDir);
	 	ans.printLatents(param.pathOutputDir);
	 	std::cerr << "end output " << std::endl;
	}
}