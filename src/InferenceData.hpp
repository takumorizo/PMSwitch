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
		InferenceData(FixedSizeMultiVector<Real, Int> _EZ, FixedSizeMultiVector<Real, Int> _ES, FixedSizeMultiVector<Real, Int> _EY,
					  FixedSizeMultiVector<Real, Int> _alpha, FixedSizeMultiVector<Real, Int> _eta,
					  FixedSizeMultiVector<Real, Int> _beta, FixedSizeMultiVector<Real, Int> _gamma,
					  FixedSizeMultiVector<Real, Int> _g);

		const FixedSizeMultiVector<Real, Int> EZ;
		const FixedSizeMultiVector<Real, Int> ES;
		const FixedSizeMultiVector<Real, Int> EY;

		const FixedSizeMultiVector<Real, Int> alpha;
		const FixedSizeMultiVector<Real, Int> eta;
		const FixedSizeMultiVector<Real, Int> beta;
		const FixedSizeMultiVector<Real, Int> gamma;

		const FixedSizeMultiVector<Real, Int> g;
	};
}

template<typename Int , typename Real >
pmswitch::InferenceData<Int, Real>::InferenceData(FixedSizeMultiVector<Real, Int> _EZ, FixedSizeMultiVector<Real, Int> _ES, FixedSizeMultiVector<Real, Int> _EY,
					  FixedSizeMultiVector<Real, Int> _alpha, FixedSizeMultiVector<Real, Int> _eta,
					  FixedSizeMultiVector<Real, Int> _beta, FixedSizeMultiVector<Real, Int> _gamma,
					  FixedSizeMultiVector<Real, Int> _g) : EZ(_EZ), ES(_ES), EY(_EY), alpha(_alpha), eta(_eta), beta(_beta), gamma(_gamma), g(_g) {
}


/*
	public functions
*/





/*
	private functions
*/


#endif
