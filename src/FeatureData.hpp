#ifndef FeatureData_h
#define FeatureData_h

#include "FixedSizeMultiVector.hpp"

#include <unordered_map>
#include <vector>
#include <string>

// I     :: Int // sample size
// map<Int, Int>    idxToPosSize
// map<Int, string> idxToSample
// J_{i} :: Int // position size
// maxJ  :: Int
// maxMl :: Int
// x_{i,j,l, }    :: 1 out of m_{l} vector // x_{i,j,l,m}
// L     :: Int // dimension of feature vector
// M_{l} :: Int // # of possible state at l-th col of feature vector
namespace pmswitch{
	template<typename Int = long long, typename Real = double>
	class FeatureData{
	public:
		FeatureData(Int I, Int maxJ, FixedSizeMultiVector<Int> Js,
					Int L, Int maxMl, FixedSizeMultiVector<Int> Ml,
					FixedSizeMultiVector<Int> X,
					std::unordered_map<Int, std::string> sampleIdxToSampleName);
		const Int I;
		const Int maxJ;
		const FixedSizeMultiVector<Int> Js;
		const Int L;
		const Int maxMl;
		const FixedSizeMultiVector<Int> Ml;
		const FixedSizeMultiVector<Int> X;

		const std::unordered_map<Int, std::string> sampleIdxToSampleName;
		// std::unordered_map<Int, Int> idxToPosSize;
	};
}


template<typename Int, typename Real>
pmswitch::FeatureData<Int, Real>::FeatureData(
				Int I, Int maxJ, FixedSizeMultiVector<Int> Js,
				Int L, Int maxMl, FixedSizeMultiVector<Int> Ml,
				FixedSizeMultiVector<Int> X,
				std::unordered_map<Int, std::string> sampleIdxToSampleName) : I(I)
														, maxJ(maxJ), Js(Js)
														, L(L), maxMl(maxMl), Ml(Ml)
														, X(X)
														, sampleIdxToSampleName(sampleIdxToSampleName){
}

#endif