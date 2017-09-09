#ifndef Inference_h
#define Inference_h

#include "DBData.hpp"
#include "FeatureData.hpp"
#include "InferenceData.hpp"
#include "PriorParameters.hpp"
#include "InputFileParser.hpp"
#include "MathUtil.hpp"
#include "Parameters.hpp"


#include <vector>
#include <iostream>
#include <assert.h>
#include <string>
#include <boost/math/special_functions/digamma.hpp>
#include <math.h>
#include <cmath>

namespace pmswitch{
	template<typename Int = long long, typename Real = double>
	class Inference{
	public:
		Inference(FeatureData<Int, Real> _fvData, DBData<Int, Real> _dbData, PriorParameters<Int, Real> prior, Int _T);
		InferenceData<Int, Real> vb(bool updateDB = true, const Real epsilon = 1e-9);
		InferenceData<Int, Real> onlineVb();
	private:
		FeatureData<Int, Real> fvData;
		DBData<Int, Real> dbData;
		PriorParameters<Int, Real> prior;
		Int T;

		FixedSizeMultiVector<Real> initEZ(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int T);
		FixedSizeMultiVector<Real> initES(Int I, FixedSizeMultiVector<Int> Js, Int maxJ);
		FixedSizeMultiVector<Real> initEY(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int N);
	};

	template<typename Int = long long, typename Real = double>
	class InferenceCreator{
	public:
		InferenceCreator();
	  	static pmswitch::Inference<Int, Real> createInference(std::string pathFV, std::string pathDB,
													   Int T, Real _alpha, Real _beta0, Real _beta1,
													   Real _gamma, Real _eta,
													   std::string alphaPath = "", std::string betaPath = "",
													   std::string gammaPath = "", std::string etaPath = "",
													   Real (*strToE)(std::string) = NULL);
	  	static pmswitch::Inference<Int, Real> createInference(pmswitch::Parameters<Int, Real> param);
	 private:
	 	static Real strToReal(std::string str);
	};
}
/*
	public functions
*/
template<typename Int, typename Real>
pmswitch::Inference<Int, Real>::Inference(FeatureData<Int, Real> _fvData, DBData<Int, Real> _dbData, PriorParameters<Int, Real> _prior, Int _T)
		:fvData(_fvData), dbData(_dbData), prior(_prior), T(_T){
}


template<typename Int, typename Real>
pmswitch::InferenceData<Int, Real> pmswitch::Inference<Int, Real>::vb(bool updateDB, const Real epsilon ){
	using namespace pmswitch;
	Int N = dbData.N; // existing cluster Number
	Int I = fvData.I; // sample size
	FixedSizeMultiVector<Int> Js = fvData.Js; // position size list for each sample
	Int maxJ = fvData.maxJ;					  // position size max
	Int L = dbData.L;						  // signature size @ each position
	assert(dbData.L == fvData.L);
	Int maxMl = dbData.maxMl;				  // feature dimension max
	FixedSizeMultiVector<Int> ml = fvData.Ml; // feature dimension list
	FixedSizeMultiVector<Int> X  = fvData.X;  // i,maxJ,l,maxMl, each data

	FixedSizeMultiVector<Real, Int> EZ = initEZ(I, Js, maxJ, T); // I, maxJ, T
	FixedSizeMultiVector<Real, Int> ES = initES(I, Js, maxJ);    // I, maxJ, 2
	FixedSizeMultiVector<Real, Int> EY = initEY(I, Js, maxJ, N); // I, maxJ, N

	FixedSizeMultiVector<Real, Int> alpha = prior.alpha; // I, T-1, 2
	FixedSizeMultiVector<Real, Int> beta  = prior.beta;  // I, 2
	FixedSizeMultiVector<Real, Int> gamma = prior.gamma; // I, N
	FixedSizeMultiVector<Real, Int> eta   = prior.eta;   // T, Ml.size(), maxMl

	FixedSizeMultiVector<Real, Int> dig_v     = math::calDirExp<Int, Real>(alpha, 1); // I, T-1, 2
	FixedSizeMultiVector<Real, Int> dig_theta = math::calDirExp<Int, Real>(beta,  0); // I, 2
	FixedSizeMultiVector<Real, Int> dig_pi    = math::calDirExp<Int, Real>(gamma, 0); // I, N
	FixedSizeMultiVector<Real, Int> dig_f     = math::calDirExp<Int, Real>(eta,   1); // T, L, maxMl
	FixedSizeMultiVector<Real, Int> g         = dbData.g; // N, l, maxMl
	FixedSizeMultiVector<Real, Int> dig_v_sig(0.0, I, T); // I, T

	FixedSizeMultiVector<bool, Int> EZFilter(true, I, maxJ, T);   math::makeFilter(EZFilter,  I, Js, FixedSizeMultiVector<Int, Int>(T, maxJ) );
	FixedSizeMultiVector<bool, Int> ESFilter(true, I, maxJ, 2);   math::makeFilter(ESFilter,  I, Js, FixedSizeMultiVector<Int, Int>(2, maxJ) );
	FixedSizeMultiVector<bool, Int> EYFilter(true, I, maxJ, N);   math::makeFilter(EYFilter,  I, Js, FixedSizeMultiVector<Int, Int>(N, maxJ) );
	FixedSizeMultiVector<bool, Int> etaFilter(true, T, L, maxMl); math::makeFilter(etaFilter, T, FixedSizeMultiVector<Int, Int>(L, T), ml);
	FixedSizeMultiVector<bool, Int> gFilter(true, N, L, maxMl);   math::makeFilter(gFilter,   N, FixedSizeMultiVector<Int, Int>(L, N), ml);

	FixedSizeMultiVector<Real, Int> lnEZ = math::applied<Int, Real>(EZ, std::log, EZFilter) ; // I, maxJ, T
	FixedSizeMultiVector<Real, Int> lnES = math::applied<Int, Real>(ES, std::log, ESFilter) ; // I, maxJ, 2
	FixedSizeMultiVector<Real, Int> lnEY = math::applied<Int, Real>(EY, std::log, EYFilter) ; // I, maxJ, N
	FixedSizeMultiVector<Real, Int> lng  = math::applied<Int, Real>(g,  std::log,  gFilter) ; // I, maxJ, N

    Real beforeLq = -1.0 * 1e50;
    Real nextLq = -1.0 * 1e50;
    bool first = true;
    Int update = 0;

    Real INF = 1e300;

    while( first || std::abs(nextLq - beforeLq)/std::abs(beforeLq) >=  epsilon){
	    beforeLq = nextLq;
	    nextLq = 0.0;
	    // E-step
	    {// EZ
		    for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++)for(Int k = 0; k < T; k++) EZ(i,j,k) = 0.0;

		    for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T; k++){
	    		for(Int kk = 0; kk <= k-1; kk++){ EZ(i,j,k) += dig_v(i,kk,1);}
	    		if(k < T-1) {EZ(i,j,k) += dig_v(i,k,0);}
	    		for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
	    			EZ(i,j,k) += ES(i,j,0) * X(i,j,l,m) * dig_f(k, l, m);
	    		}
		    }
    		math::subMax<Int, Real>(EZ, 1);
    		math::apply<Int, Real>(EZ, std::exp);
    		math::filter<Int, Real>(EZ, EZFilter, 0.0);
    		math::norm<Int, Real>(EZ, 1);
    		lnEZ = math::applied<Int, Real>(EZ, std::log, EZFilter);
	    }
	    {// ES
	    	std::vector<Real> tmpXLogY = {0, 0};
	    	for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++){ES(i,j,0) = 0.0; ES(i,j,1) = 0.0;}
	    	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++ ){
	    		ES(i,j,0) += dig_theta(i, 0);
	    		for(Int k = 0; k < T; k++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
	    			ES(i,j,0) += EZ(i,j,k) * X(i,j,l,m) * dig_f(k, l, m);
	    		}
	    		ES(i,j,1) += dig_theta(i,1);
	    		for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
	    			tmpXLogY[1] = EY(i,j,n) * X(i,j,l,m) * lng(n,l,m);
	    			ES(i,j,1) += tmpXLogY[!std::isnan(tmpXLogY[1])];
	    		}
	    	}
    		math::subMax<Int, Real>(ES, 1);
    		math::apply<Int, Real>(ES, std::exp);
    		math::filter<Int, Real>(ES, ESFilter, 0.0);
    		math::norm<Int, Real>(ES, 1);
    		lnES = math::applied<Int, Real>(ES, std::log, ESFilter);
	    }
		{// EY // I, maxJ, N
			std::vector<Real> tmpXLogY = {0, 0};
	    	for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++)for(Int n = 0; n < N; n++){EY(i,j,n) = 0.0;}

	    	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int n = 0; n < N; n++){
	    		EY(i,j,n) += dig_pi(i,n);
	    		for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
	    			tmpXLogY[1] = X(i,j,l,m) * ES(i,j,1) * lng(n,l,m);
	    			EY(i,j,n) += tmpXLogY[!std::isnan(tmpXLogY[1])];
	    		}
	    	}
    		math::subMax<Int, Real>(EY, 1);
    		math::apply<Int, Real>(EY, std::exp);
    		math::filter<Int, Real>(EY, EYFilter, 0.0);
    		math::norm<Int, Real>(EY, 1);
    		lnEY = math::applied<Int, Real>(EY, std::log, EYFilter);
	    }

	    // M-step
	    {// FixedSizeMultiVector<Real, Int> eta   = prior.eta;   // T, Ml.size(), maxMl
	   		eta = prior.eta;
	   		for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T; k++){
	   			for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
	   				eta(k,l,m) += X(i,j,l,m) * EZ(i,j,k) * ES(i,j,0);
	   			}
	   		}
	   		dig_f = math::calDirExp<Int, Real>(eta, 1, etaFilter);
	    }
	    {// FixedSizeMultiVector<Real, Int> alpha = prior.alpha; // I, T-1, 2
			alpha = prior.alpha; // I, T-1, 2
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T-1; k++){
				alpha(i,k,0) += EZ(i,j,k);
				for(Int t = k+1; t < T; t++){ alpha(i,k,1) += EZ(i,j,t);}
			}
			dig_v = math::calDirExp<Int, Real>(alpha, 1); // I, T-1, 2
	    }
	    {// FixedSizeMultiVector<Real, Int> beta  = prior.beta;  // I, 2
	    	beta  = prior.beta;  // I, 2
	    	for(Int i = 0; i < I; i++){
	    		for(Int j = 0; j < Js(i); j++){
		    		beta(i, 0) += ES(i,j,0); beta(i, 1) += ES(i,j,1);
	    		}
	    	}
	    	dig_theta = math::calDirExp<Int, Real>(beta,  0);
	    }
		{// FixedSizeMultiVector<Real, Int> gamma = prior.gamma; // I, N
			gamma = prior.gamma; // I, N
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int n = 0; n < N; n++){
					gamma(i,n) += EY(i,j,n);
	    	}
    	 	dig_pi    = math::calDirExp<Int, Real>(gamma, 0); // I, N
		}

		if(updateDB){ // update < 275 -> pass ,update < 276 -> not pass
			for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < maxMl; m++){ g(n,l,m) = 0.0;}

			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++)for(Int n = 0; n < N; n++){
				g(n,l,m) += X(i,j,l,m) * EY(i,j,n) * ES(i,j,1);
			}
    		math::filter(g, gFilter, (Real)0);
    		math::norm<Int, Real>(g, 1);
    		lng  = math::applied<Int, Real>(g,  std::log,  gFilter) ; // I, maxJ, N
		}

		{// nextLq = 0.0
			nextLq = 0.0;
			nextLq += math::calELogDir<Int, Real>(dig_f,     prior.eta,    1, etaFilter); // T, L, maxMl
			nextLq += math::calELogDir<Int, Real>(dig_v,     prior.alpha,  1); // I, T-1, 2
			nextLq += math::calELogDir<Int, Real>(dig_theta, prior.beta,   0); // I, 2
			nextLq += math::calELogDir<Int, Real>(dig_pi,    prior.gamma,  0); // I, N

	    	std::vector<Real> tmpXLogY = {0, 0};
			{//dig_v_sig
				for(Int i = 0; i < I; i++)for(Int k = 0; k < T; k++){ dig_v_sig(i,k) = 0.0; }
				for(Int i = 0; i < I; i++)for(Int k = 0; k < T; k++){
					if(k < T-1){ dig_v_sig(i, k) += dig_v(i, k, 0);}
					for(Int t = 0; t <= k-1; t++ ){ dig_v_sig(i, k) += dig_v(i, t, 1);}
				}
			}

			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				for(Int k = 0; k < T; k++){
					nextLq += EZ(i,j,k) * dig_v_sig(i,k);
				}
			}
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				nextLq += ES(i,j,0) * dig_theta(i,0) + ES(i,j,1) * dig_theta(i,1);
			}
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				for(Int n = 0; n < N; n++){
					nextLq += EY(i,j,n) * dig_pi(i,n);
				}
			}
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T; k++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
					nextLq += EZ(i,j,k) * ES(i,j,0) * X(i,j,l,m) * dig_f(k,l,m);
			}

			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
					tmpXLogY[1] = EY(i,j,n) * ES(i,j,1) * X(i,j,l,m) * lng(n,l,m);
					nextLq     += tmpXLogY[ ! (std::isnan(tmpXLogY[1]) or std::isinf(tmpXLogY[1]) ) ];
			}

			nextLq -= math::calELogDir<Int, Real>(dig_f,     eta,    1, etaFilter);  // T, L, maxMl
			nextLq -= math::calELogDir<Int, Real>(dig_v,     alpha,  1); 			 // I, T-1, 2
			nextLq -= math::calELogDir<Int, Real>(dig_theta, beta,   0); 			 // I, 2
			nextLq -= math::calELogDir<Int, Real>(dig_pi,    gamma,  0); 			 // I, N

			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				for(Int k = 0; k < T; k++){
					tmpXLogY[1] = EZ(i,j,k) * lnEZ(i,j,k);;
					nextLq -= tmpXLogY[!std::isnan(tmpXLogY[1])];
				}
			}

			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				tmpXLogY[1] = ES(i,j,0) * lnES(i,j,0);
				nextLq -= tmpXLogY[!std::isnan(tmpXLogY[1])];
				tmpXLogY[1] = ES(i,j,1) * lnES(i,j,1) ;
				nextLq -= tmpXLogY[!std::isnan(tmpXLogY[1])];
				if(std::isnan(nextLq)){
					std::cerr << "(i,j) : " << i << ", " << j << std::endl;
			    	std::cerr << "	ES(i,j,0): " << ES(i,j,0) << std::endl;
			    	std::cerr << "	ES(i,j,1): " << ES(i,j,1) << std::endl;
			    	std::cerr << "	lnES(i,j,0): " << lnES(i,j,0) << std::endl;
			    	std::cerr << "	lnES(i,j,1): " << lnES(i,j,1) << std::endl;
			    	assert(false);
				}
			}
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				for(Int n = 0; n < N; n++){
					tmpXLogY[1] = EY(i,j,n) * lnEY(i,j,n);
					nextLq -= tmpXLogY[!std::isnan(tmpXLogY[1])];
				}
			}
		}

		// if(updateDB){std::cerr << "ok vb, with g update, update : ";}
		// else{std::cerr << "ok vb, update : ";}
    	// std::cerr << update << ", beforeLq :" << beforeLq << " nextLq : "  << nextLq << ", nextLq-beforeLq/beforeLq :" << (nextLq - beforeLq)/std::abs(beforeLq) << std::endl;
    	assert( !(std::isnan(beforeLq) || std::isnan(nextLq) ) ) ;
    	assert( !(std::isinf(beforeLq) || std::isinf(nextLq) ) ) ;
		if( not first && (beforeLq - nextLq)/std::abs(beforeLq) >=  epsilon ){
	    	std::cerr << "decrease lower bound @ vb, update : " << update << std::endl;
	    	std::cerr << "vb, update : " << update << ", beforeLq :" << beforeLq << " nextLq : "  << nextLq << ", nextLq-beforeLq/beforeLq :" << (nextLq - beforeLq)/std::abs(beforeLq) << std::endl;
		    assert( (beforeLq - nextLq)/std::abs(beforeLq) < epsilon  );
	    }
		if(first) first = false;
		update++;
	}
	pmswitch::InferenceData<Int, Real> ans(EZ, ES, EY, alpha, eta, beta, gamma, g, nextLq);
	return ans;
}

// EZ_{i,j, }     :: T-dimensional non-negative simplex // EZ_{i,j,k   }
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::initEZ(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int T){
	pmswitch::FixedSizeMultiVector<Real> EZ(0.0, I, maxJ, T);
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++) EZ(i,j,0) = 1.0;
	return EZ;
}

// ES_{i, }	   :: 2-dimensional non-negative simplex // ES_{i,{0,1} }
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::initES(Int I, FixedSizeMultiVector<Int> Js, Int maxJ){
	pmswitch::FixedSizeMultiVector<Real> ES(0.0, I, maxJ, 2);
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++) ES(i,j,1) = 1.0;
	return ES;
}

// EY_{i, }       :: N-dimensional non-negative simplex // EY_{i,n}
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::initEY(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int N){
	pmswitch::FixedSizeMultiVector<Real> EY(0.0, I, maxJ, N);
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++) EY(i,j,0) = 1.0;
	return EY;
}

/*
	public functions
*/
template<typename Int, typename Real>
pmswitch::InferenceCreator<Int, Real>::InferenceCreator(){

}

template<typename Int, typename Real>
pmswitch::Inference<Int, Real> pmswitch::InferenceCreator<Int, Real>::createInference(std::string pathFV, std::string pathDB,
													   Int _T, Real _alpha, Real _beta0, Real _beta1,
													   Real _gamma, Real _eta,
   													   std::string alphaPath, std::string betaPath,
													   std::string gammaPath, std::string etaPath,
													   Real (*strToE)(std::string)){
	using namespace pmswitch;
	InputFileParser<Int, Real> parser;

	FeatureData<Int, Real> fvData = parser.parseFeatureFile(pathFV);
	DBData<Int, Real> dbData = parser.parseDBFile(pathDB);
	// PriorParametersCreator<Int, Real> priorCreator;
	PriorParameters<Int, Real> prior = PriorParametersCreator<Int, Real>::createPriorParameters(fvData,  dbData, _T, _alpha, _beta0, _beta1, _gamma, _eta,
 																								alphaPath,  betaPath, gammaPath,  etaPath, strToE);
	Inference<Int, Real> inference(fvData, dbData, prior, _T);
	return inference;
}


template<typename Int, typename Real>
pmswitch::Inference<Int, Real> pmswitch::InferenceCreator<Int, Real>::createInference(pmswitch::Parameters<Int, Real> param){
	return pmswitch::InferenceCreator<Int, Real>::createInference(param.pathFV, param.pathDB, param.T, param.alpha, param.beta0, param.beta1,
																  param.gamma,  param.eta, param.alphaPath, param.betaPath, param.gammaPath, param.etaPath,
																  pmswitch::InferenceCreator<Int, Real>::strToReal);
}

template<typename Int, typename Real>
Real pmswitch::InferenceCreator<Int, Real>::strToReal(std::string str){
	return (Real) std::stold(str);
}


/*
	private functions
*/
#endif
