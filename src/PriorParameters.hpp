#ifndef PriorParameters_h
#define PriorParameters_h

#include "FixedSizeMultiVector.hpp"
#include "FeatureData.hpp"
#include "DBData.hpp"
#include "InputFileParser.hpp"
#include "Inference.hpp"

#include <vector>
#include <iostream>
#include <assert.h>

/*
###EZ_{i,j, }###
	v_{i, }        :: T-dimensional real in [0,1] //v_{i, k}
	dig_v_{i,}     :: T-dimensional real
	\alpha{i, }    :: T-dimensional real

	f_{k, }        :: m_{l}-dimensional non-negative simplex // f_{k,l}
	dig_f_{k, }    :: m_{l}-dimensional real vector          // dig_f_{k, l}
	\eta_{k,}	   :: m_{l}-dimensional real vector 		 // eta_{k, l}
################

###ES_{i, }###
	\theta_{i}     :: real in [0,1]
	dig_\theta_{i,}:: 2-dimensional real vector
	\beta_{i, }    :: 2-dimensional real vector
################

###EY_{i, }###
	\pi_{i}        :: real in [0,1]
	dig_\pi_{i, }  :: N-dimensional real
	\gamma_{i, }   :: N-dimensional real

	g_{n, }        :: m_{l}-dimensional non-negative simplex // g_{n,l}
################
*/
namespace pmswitch{
	template<typename Int = long long, typename Real = double>
	class PriorParameters{
	public:
		PriorParameters(pmswitch::FixedSizeMultiVector<Real> alpha,
						pmswitch::FixedSizeMultiVector<Real> beta,
						pmswitch::FixedSizeMultiVector<Real> gamma,
						pmswitch::FixedSizeMultiVector<Real> eta);
		const pmswitch::FixedSizeMultiVector<Real> alpha;
		const pmswitch::FixedSizeMultiVector<Real> beta;
		const pmswitch::FixedSizeMultiVector<Real> gamma;
		const pmswitch::FixedSizeMultiVector<Real> eta;
	private:
	};

	template<typename Int = long long, typename Real = double>
	class PriorParametersCreator{
	public:
		PriorParametersCreator();
/*
	\alpha{i, }    :: T-dimensional real
	\eta_{k,}	   :: m_{l}-dimensional real vector 		 // eta_{k, l}
	\beta_{i, }    :: 2-dimensional real vector
	\gamma_{i, }   :: N-dimensional real
*/
		pmswitch::PriorParameters<Int, Real> testFunc(pmswitch::FeatureData<Int,Real> fvData){
			pmswitch::FixedSizeMultiVector<Real> alpha(0,1,1);
			pmswitch::FixedSizeMultiVector<Real> beta(0,1,1);
			pmswitch::FixedSizeMultiVector<Real> gamma(0,1,1);
			pmswitch::FixedSizeMultiVector<Real> eta(0,1,1);
			PriorParameters<Int, Real> prior(alpha,beta,gamma,eta);
			return prior;
		}

		pmswitch::PriorParameters<Int, Real> createPriorParameters(
														   Int I, Int N, Int T, pmswitch::FixedSizeMultiVector<Int> Ml,
														   Real _alpha,
														   Real _beta0, Real _beta1,
														   Real _gamma, Real _eta);
		pmswitch::PriorParameters<Int, Real> createPriorParameters(
														   pmswitch::FeatureData<Int,Real> fvData,
														   pmswitch::DBData<Int,Real> dbData,
														   Int T,
														   Real _alpha,
														   Real _beta0, Real _beta1,
														   Real _gamma, Real _eta);
	};
}
template<typename Int, typename Real>
pmswitch::PriorParameters<Int, Real>::PriorParameters(pmswitch::FixedSizeMultiVector<Real> _alpha,
													  pmswitch::FixedSizeMultiVector<Real> _beta,
													  pmswitch::FixedSizeMultiVector<Real> _gamma,
													  pmswitch::FixedSizeMultiVector<Real> _eta) : alpha(_alpha), beta(_beta), gamma(_gamma), eta(_eta){

}

template<typename Int, typename Real>
pmswitch::PriorParametersCreator<Int, Real>::PriorParametersCreator(){

}

template pmswitch::PriorParameters<long long, double> pmswitch::PriorParametersCreator<long long, double>::createPriorParameters(
														   long long I, long long N, long long T, pmswitch::FixedSizeMultiVector<long long> Ml,
														   double _alpha,
														   double _beta0, double _beta1,
														   double _gamma, double _eta);

template<typename Int, typename Real>
pmswitch::PriorParameters<Int, Real> pmswitch::PriorParametersCreator<Int, Real>::createPriorParameters(
														   Int I, Int N, Int T, pmswitch::FixedSizeMultiVector<Int> Ml,
														   Real _alpha,
														   Real _beta0, Real _beta1,
														   Real _gamma, Real _eta){
	assert( T-1 > 0 );
	Int maxMl = 0; for(Int l = 0; l < Ml.size(); l++) maxMl = std::max(maxMl, Ml[l]);
	pmswitch::FixedSizeMultiVector<Real, Int> alpha(_alpha, I, T-1, 2);
	pmswitch::FixedSizeMultiVector<Real, Int> beta(0.0, I, 2);
	pmswitch::FixedSizeMultiVector<Real, Int> gamma(_gamma, I, N);
	pmswitch::FixedSizeMultiVector<Real, Int> eta(_eta, T, Ml.size(), maxMl);
	{
		for(Int i = 0; i < I; i++)for(Int t = 0; t < T-1; t++){
			alpha(i, t, 0) = 1.0; alpha(i, t, 1) = _alpha;
		}
		for(Int i = 0; i < I; i++){
			beta(i,0) = _beta0; beta(i,1) = _beta1;
		}
		for(Int t = 0; t < T; t++)for(Int l = 0; l < Ml.size(); l++)for(Int m = 0; m < Ml(l); m++){
			eta(t,l,m) = _eta;
		}
	}
	pmswitch::PriorParameters<Int, Real> prior(alpha, beta, gamma, eta);
	return prior;
}


template<typename Int, typename Real>
pmswitch::PriorParameters<Int, Real> pmswitch::PriorParametersCreator<Int, Real>::createPriorParameters(
												   pmswitch::FeatureData<Int,Real> fvData,
												   pmswitch::DBData<Int,Real> dbData,
												   Int T,
												   Real _alpha, Real _beta0, Real _beta1, Real _gamma, Real _eta){
	return createPriorParameters(fvData.I, dbData.N, T, fvData.Ml, _alpha, _beta0, _beta1, _gamma, _eta);
}


#endif
