#ifndef InferenceData_h
#define InferenceData_h

#include "DBData.hpp"
#include "FeatureData.hpp"
#include "InferenceData.hpp"
#include "FixedSizeMultiVector.hpp"


#include <vector>
#include <iostream>
#include <assert.h>

namespace pmswitch{
	template<typename Int = long long, typename Real = double>
	class InferenceData{
	public:
		InferenceData(FixedSizeMultiVector<double, Int> _EZ, FixedSizeMultiVector<double, Int> _ES, FixedSizeMultiVector<double, Int> _EY,
					  FixedSizeMultiVector<double, Int> _alpha, FixedSizeMultiVector<double, Int> _eta,
					  FixedSizeMultiVector<double, Int> _beta, FixedSizeMultiVector<double, Int> _gamma,
					  FixedSizeMultiVector<double, Int> _g);
		~InferenceData();

		const FixedSizeMultiVector<double, Int> EZ;
		const FixedSizeMultiVector<double, Int> ES;
		const FixedSizeMultiVector<double, Int> EY;

		const FixedSizeMultiVector<double, Int> alpha;
		const FixedSizeMultiVector<double, Int> eta;
		const FixedSizeMultiVector<double, Int> beta;
		const FixedSizeMultiVector<double, Int> gamma;

		const FixedSizeMultiVector<double, Int> g;
	};
}

template<typename Int , typename Real >
pmswitch::InferenceData<Int, Real>::InferenceData(FixedSizeMultiVector<double, Int> _EZ, FixedSizeMultiVector<double, Int> _ES, FixedSizeMultiVector<double, Int> _EY,
					  FixedSizeMultiVector<double, Int> _alpha, FixedSizeMultiVector<double, Int> _eta,
					  FixedSizeMultiVector<double, Int> _beta, FixedSizeMultiVector<double, Int> _gamma,
					  FixedSizeMultiVector<double, Int> _g) : EZ(_EZ), ES(_ES), EY(_EY), alpha(_alpha), eta(_eta), beta(_beta), gamma(_gamma), g(g) {
}


/*
	public functions
*/





/*
	private functions
*/


#endif
