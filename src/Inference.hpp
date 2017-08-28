#ifndef Inference_h
#define Inference_h

#include "DBData.hpp"
#include "FeatureData.hpp"
#include "InferenceData.hpp"
#include "PriorParameters.hpp"
#include "InputFileParser.hpp"

#include <vector>
#include <iostream>
#include <assert.h>
#include <string>
#include <boost/math/special_functions/digamma.hpp>
#include <math.h>

namespace pmswitch{
	namespace math{
		template<typename Int = long long, typename Real = double>
		FixedSizeMultiVector<Real> initEZ(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int T);

		template<typename Int = long long, typename Real = double>
		FixedSizeMultiVector<Real> initES(Int I, FixedSizeMultiVector<Int> Js, Int maxJ);

		template<typename Int = long long, typename Real = double>
		FixedSizeMultiVector<Real> initEY(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int N);

		template<typename Int = long long, typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int baseIdx, Int dim);

		template<typename Int = long long, typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec);

		template<typename Int = long long, typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int baseIdx, Int dim,
				 const pmswitch::FixedSizeMultiVector<bool,Int> &filter);

		template<typename Int = long long, typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec,
				 const pmswitch::FixedSizeMultiVector<bool,Int> &filter);

		template<typename Int = long long, typename Real = double>
		Real dot(const pmswitch::FixedSizeMultiVector<Real,Int> &a, const pmswitch::FixedSizeMultiVector<Real,Int> &b);

		template<typename Int = long long, typename Real = double>
		Real dot(const pmswitch::FixedSizeMultiVector<Real,Int> &a, const pmswitch::FixedSizeMultiVector<Real,Int> &b,
				 const pmswitch::FixedSizeMultiVector<Real,Int> &filter);

		template<typename Int = long long, typename Real = double>
		void subMax(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim, Real INF = 1e300);

		template<typename Int = long long, typename Real = double>
		void subMin(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim);

		template<typename Int = long long, typename Real = double>
		void norm(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim);

		template<typename Int = long long, typename Real = double>
		void norm(pmswitch::FixedSizeMultiVector<Real,Int> &vec);

		template<typename Int = long long, typename Real = double>
		void apply(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real) );

		template<typename Int = long long, typename Real = double>
		pmswitch::FixedSizeMultiVector<Real,Int> applied(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real) );

		template<typename Int = long long, typename Real = double>
		void apply(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real),
				   const pmswitch::FixedSizeMultiVector<bool,Int> &filter);

		template<typename Int = long long, typename Real = double>
		pmswitch::FixedSizeMultiVector<Real,Int> applied(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real),
				   									   const pmswitch::FixedSizeMultiVector<bool,Int> &filter);

		template<typename Int = long long, typename Real = double>
		void filter(pmswitch::FixedSizeMultiVector<Real,Int> &vec, const FixedSizeMultiVector<bool, Int> &filterVec, Real defaultVal);

		template<typename Int = long long, class... Args>
		void makeFilter(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int range, const Args&... args);

		template<typename Int = long long, class... Args>
		void makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int, Int> &range, const Args&... args);

		template<typename Int = long long, class... Args>
		void makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int, Int> &range);


		template<typename Int = long long, typename Real = double>
		pmswitch::FixedSizeMultiVector<Real,Int> calDirExp(const pmswitch::FixedSizeMultiVector<Real,Int> &alpha, Int dim);

		template<typename Int = long long, typename Real = double>
		pmswitch::FixedSizeMultiVector<Real,Int> calDirExp(const pmswitch::FixedSizeMultiVector<Real,Int> &alpha);

		template<typename Int = long long, typename Real = double>
		pmswitch::FixedSizeMultiVector<Real,Int> calDirExp(const pmswitch::FixedSizeMultiVector<Real,Int> &alpha, Int dim,
														   const pmswitch::FixedSizeMultiVector<bool,Int> &filter);

		template<typename Int = long long, typename Real = double>
		pmswitch::FixedSizeMultiVector<Real,Int> calDirExp(const pmswitch::FixedSizeMultiVector<Real,Int> &alpha,
														   const pmswitch::FixedSizeMultiVector<bool,Int> &filter);

		template<typename Int = long long, typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param, Int dim);

		template<typename Int = long long, typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param);

		template<typename Int = long long, typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param, Int dim,
						const pmswitch::FixedSizeMultiVector<bool,Int> &filter);

		template<typename Int = long long, typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param,
						const pmswitch::FixedSizeMultiVector<bool,Int> &filter);
	}

	template<typename Int = long long, typename Real = double>
	class Inference{
	public:
		Inference(FeatureData<Int, Real> _fvData, DBData<Int, Real> _dbData, PriorParameters<Int, Real> prior);
		InferenceData<Int, Real> vb(Int T, bool updateDB = true);
		InferenceData<Int, Real> onlineVb();
	private:
		FeatureData<Int, Real> fvData;
		DBData<Int, Real> dbData;
		PriorParameters<Int, Real> prior;
	};

	template<typename Int = long long, typename Real = double>
	class InferenceCreator{
	public:
		InferenceCreator();
		pmswitch::Inference<Int, Real> createInference(std::string pathFV, std::string pathDB,
													   Int T, Real _alpha, Real beta0, Real beta1,
													   Real _gamma, Real _eta );
	};
}



/*
################
namespace pmswitch{
	namespace math{
		template<typename Int = long long, typename Real = double>
		FixedSizeMultiVector<Real> initEZ(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int T);

		template<typename Int = long long, typename Real = double>
		FixedSizeMultiVector<Real> initES(Int I, Int N);

		template<typename Int = long long, typename Real = double>
		FixedSizeMultiVector<Real> initEY(Int I, Int N);
	}
}
################
*/

// EZ_{i,j, }     :: T-dimensional non-negative simplex // EZ_{i,j,k   }
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::initEZ(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int T){
	pmswitch::FixedSizeMultiVector<Real> EZ(0.0, I, maxJ, T);
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++) EZ(i,j,0) = 1.0;
	return EZ;
}

// ES_{i, }	   :: 2-dimensional non-negative simplex // ES_{i,{0,1} }
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::initES(Int I, FixedSizeMultiVector<Int> Js, Int maxJ){
	pmswitch::FixedSizeMultiVector<Real> ES(0.0, I, maxJ, 2);
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++) ES(i,j,1) = 1.0;
	return ES;
}

// EY_{i, }       :: N-dimensional non-negative simplex // EY_{i,n}
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::initEY(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int N){
	pmswitch::FixedSizeMultiVector<Real> EY(0.0, I, maxJ, N);
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++) EY(i,j,0) = 1.0;
	return EY;
}

template<typename Int = long long, typename Real = double>
void subMax(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim, Real INF){
	Int outN = vec.size() / vec.volume(dim);
	Int inN  = vec.volume(dim);
	Real localMax = -INF;
	Int baseIdx = 0;
	for(Int i = 0; i < outN; i++){
		baseIdx = i * inN;
		localMax = -INF;
		for(Int j = 0; j < inN; j++){ localMax = std::max(vec[baseIdx + j], localMax) ;}
		for(Int j = 0; j < inN; j++){ vec[baseIdx + j] -= localMax;}
	}
}

template<typename Int = long long, typename Real = double>
void subMin(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim, Real INF){
	Int outN = vec.size() / vec.volume(dim);
	Int inN  = vec.volume(dim);
	Real localMin = INF;
	Int baseIdx = 0;
	for(Int i = 0; i < outN; i++){
		baseIdx = i * inN;
		localMin = INF;
		for(Int j = 0; j < inN; j++){ localMin = std::min(vec[baseIdx + j], localMin) ;}
		for(Int j = 0; j < inN; j++){ vec[baseIdx + j] -= localMin;}
	}
}

template<typename Int, typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int baseIdx, Int dim){
	Int inN = vec.volume(dim);
	Real sum = 0.0;
	for(Int j = 0; j < inN; j++){ sum += vec[baseIdx + j];}
	return sum;
}

template<typename Int, typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec){
	Real sum = 0.0;
	for(Int i = 0; i < vec.size(); i++){ sum += vec[i];}
	return sum;
}

template<typename Int, typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int baseIdx, Int dim,
							const pmswitch::FixedSizeMultiVector<bool,Int> &filter){
	assert(vec.size() == filter.size() && vec.dim() == filter.dim());
	#ifdef NDEBUG
	#else
		for(Int d = 0; d < vec.dim(); d++){ assert(vec.size(d) == filter.size(d)); }
	#endif

	Int inN = vec.volume(dim);
	Real sum = 0.0;
	for(Int j = 0; j < inN; j++)if(filter[baseIdx + j]){ sum += vec[baseIdx + j];}
	return sum;
}

template<typename Int, typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real,Int> &vec,
						 const pmswitch::FixedSizeMultiVector<bool,Int> &filter){
	assert(vec.size() == filter.size() && vec.dim() == filter.dim());
	#ifdef NDEBUG
	#else
		for(Int d = 0; d < vec.dim(); d++){ assert(vec.size(d) == filter.size(d)); }
	#endif

	Real sum = 0.0;
	for(Int i = 0; i < vec.size(); i++){ sum += (Int)(filter[i]) * vec[i];}
	return sum;
}


template<typename Int, typename Real>
Real pmswitch::math::dot(const pmswitch::FixedSizeMultiVector<Real,Int> &a, const pmswitch::FixedSizeMultiVector<Real,Int> &b){
	assert(a.size() == b.size() && a.dim() == b.dim());
	#ifdef NDEBUG
	#else
		for(Int d = 0; d < a.dim(); d++){ assert(a.size(d) == b.size(d)); }
	#endif
	Real dot = 0.0;
	for(Int i = 0; i < a.size(); i++){ dot += a(i) * b(i);}
	return dot;
}

template<typename Int, typename Real>
Real pmswitch::math::dot(const pmswitch::FixedSizeMultiVector<Real,Int> &a, const pmswitch::FixedSizeMultiVector<Real,Int> &b,
							const pmswitch::FixedSizeMultiVector<Real,Int> &filter){
	assert(a.size() == b.size()      && a.dim() == b.dim());
	assert(a.size() == filter.size() && a.dim() == filter.dim());
	#ifdef NDEBUG
	#else
		for(Int d = 0; d < a.dim(); d++){ assert(a.size(d) == b.size(d)); }
		for(Int d = 0; d < a.dim(); d++){ assert(a.size(d) == filter.size(d)); }
	#endif
	Real dot = 0.0;
	for(Int i = 0; i < a.size(); i++){ dot += (Int)(filter[i]) * a[i] * b[i];}
	return dot;
}

template<typename Int, typename Real>
void pmswitch::math::norm(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim){
	Int outN = vec.size() / vec.volume(dim);
	Int inN  = vec.volume(dim);
	Real sum = 0.0;
	Int baseIdx = 0;
	for(Int i = 0; i < outN; i++){
		baseIdx = i * inN;
		sum = 0.0;
		for(Int j = 0; j < inN; j++){ sum += vec[baseIdx + j];}
		for(Int j = 0; j < inN; j++){ vec[baseIdx + j] /= sum;}
	}
}

template<typename Int, typename Real>
void pmswitch::math::norm(pmswitch::FixedSizeMultiVector<Real,Int> &vec){
	Real sum = 0.0;
	for(Int i = 0; i < vec.size(); i++){ sum += vec[i];}
	for(Int i = 0; i < vec.size(); i++){ vec[i] /= sum;}
}

template<typename Int, typename Real>
void pmswitch::math::apply(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real) ){
	for(Int i = 0; i < vec.size(); i++){ vec[i] = func(vec[i]);}
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real,Int> pmswitch::math::applied(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real) ){
	pmswitch::FixedSizeMultiVector<Real,Int> ans = vec;
	apply(ans, func);
	return ans;
}

template<typename Int, typename Real>
void pmswitch::math::apply(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real),
							  const FixedSizeMultiVector<bool,Int> &filter){
	{// size, dim check in debug
		assert(vec.size() == filter.size() && vec.dim() == filter.dim());
		#ifdef NDEBUG
		#else
			for(Int d = 0; d < vec.dim(); d++){ assert(vec.size(d) == filter.size(d)); }
		#endif
	}
	for(Int i = 0; i < vec.size(); i++)if(filter[i]){ vec[i] = func(vec[i]);}
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real,Int> pmswitch::math::applied(const pmswitch::FixedSizeMultiVector<Real,Int> &vec, Real (*func)(Real),
							  									  const FixedSizeMultiVector<bool,Int> &filter){
	pmswitch::FixedSizeMultiVector<Real,Int> ans = vec;
	apply(ans, func, filter);
	return ans;
}

template<typename Int, typename Real>
void pmswitch::math::filter(pmswitch::FixedSizeMultiVector<Real,Int> &vec, const FixedSizeMultiVector<bool,Int> &filterVec, Real defaultVal){
	assert(vec.size() == filterVec.size());
	for(Int i = 0; i < vec.size(); i++)if( not filterVec[i] ){ vec[i] = defaultVal; }
}


template<typename Int, class... Args>
void pmswitch::math::makeFilter(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int range, const Args&... args){
	for(Int i = 0; i < vec.size(0); i++){
		pmswitch::math::makeFilterImpl(vec, (long long)1, i*vec.volume(0), i, (i < range), args...);
	}
}

template<typename Int, class... Args>
void pmswitch::math::makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx,  bool val,
  									   const pmswitch::FixedSizeMultiVector<Int, Int> &range, const Args&... args){
	for(Int i = 0; i < vec.size(n); i++){
		pmswitch::math::makeFilterImpl(vec, n+1, idx+i*vec.volume(n), i, val && (i < range(prevIdx)), args...);
	}
}

template<typename Int, class... Args>
void pmswitch::math::makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int, Int> &range){
	for(Int i = 0; i < vec.size(n); i++){
		vec[idx+i*vec.volume(n)] = ( val && ( i < range(prevIdx) ) );
	}
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real,Int> pmswitch::math::calDirExp(const pmswitch::FixedSizeMultiVector<Real,Int> &alpha, Int dim){

	pmswitch::FixedSizeMultiVector<Real,Int> ansDig = alpha;
	Int baseIdx = 0;
	Real sum = 0.0, digSum = 0.0;
	Int outN = alpha.size() / alpha.volume(dim);
	Int inN	= alpha.volume(dim);
	for(Int i = 0; i < outN; i++){
		sum = 0.0;
		baseIdx = i * inN;
		for(Int j = 0; j < inN; j++){ sum += alpha[baseIdx + j];}
		digSum =  boost::math::digamma(sum);
		for(Int j = 0; j < inN; j++){
			ansDig[baseIdx + j] = boost::math::digamma(ansDig[baseIdx + j]) - digSum;
		}
	}
	return ansDig;
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real,Int> pmswitch::math::calDirExp( const pmswitch::FixedSizeMultiVector<Real,Int> &alpha){
	pmswitch::FixedSizeMultiVector<Real,Int> ansDig = alpha;
	Real sum = 0.0, digSum = 0.0;
	sum = 0.0;
	for(Int j = 0; j < alpha.size(); j++){ sum += alpha[j];}
	digSum =  boost::math::digamma(sum);
	for(Int j = 0; j < alpha.size(); j++){ ansDig[j] = boost::math::digamma(ansDig[j]) - digSum;}
	return ansDig;
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real,Int> pmswitch::math::calDirExp(const pmswitch::FixedSizeMultiVector<Real,Int> &alpha, Int dim,
																   const pmswitch::FixedSizeMultiVector<bool,Int> &filter){

	pmswitch::FixedSizeMultiVector<Real,Int> ansDig = alpha;
	Int baseIdx = 0;
	Real sum = 0.0, digSum = 0.0;
	Int outN = alpha.size() / alpha.volume(dim);
	Int inN	= alpha.volume(dim);
	for(Int i = 0; i < outN; i++){
		sum = 0.0;
		baseIdx = i * inN;
		for(Int j = 0; j < inN; j++)if(filter[baseIdx + j]){ sum += alpha[baseIdx + j]; }
		digSum =  boost::math::digamma(sum);
		for(Int j = 0; j < inN; j++)if(filter[baseIdx + j]){
			ansDig[baseIdx + j] =  boost::math::digamma(ansDig[baseIdx + j]) - digSum;
		}
	}
	return ansDig;
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real,Int> pmswitch::math::calDirExp( const pmswitch::FixedSizeMultiVector<Real,Int> &alpha, const pmswitch::FixedSizeMultiVector<bool,Int> &filter){
	pmswitch::FixedSizeMultiVector<Real,Int> ansDig = alpha;
	Real sum = 0.0, digSum = 0.0;
	sum = 0.0;
	for(Int j = 0; j < alpha.size(); j++)if(filter[j]){ sum += alpha[j];}
	digSum =  boost::math::digamma(sum);
	for(Int j = 0; j < alpha.size(); j++)if(filter[j]){ ansDig[j] = boost::math::digamma(ansDig[j]) - digSum;}
	return ansDig;
}

template<typename Int, typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param, Int dim){
	Int outN = param.size() / param.volume(dim);
	Int inN  = param.volume(dim);
	{// size, dim check in debug
		assert(logTheta.size() == param.size() && logTheta.dim() == param.dim());
		#ifdef NDEBUG
		#else
			for(Int d = 0; d < logTheta.dim(); d++){ assert(logTheta.size(d) == param.size(d)); }
		#endif
	}
	Real ans = 0.0;
	Int baseIdx = 0;
	for(Int i = 0; i < outN; i++){
		baseIdx = i * param.volume(dim);
		Real local_lgam_sum_param = lgamma(pmswitch::math::sum(param, baseIdx, dim));
		Real local_sum_lgam_param = 0.0;
		Real local_sum_par_lt     = 0.0;
		{
			for(Int j = 0; j < inN; j++){ local_sum_lgam_param -= lgamma(param[baseIdx + j]);}
			for(Int j = 0; j < inN; j++){ local_sum_par_lt += (param[baseIdx + j]-1.0) * logTheta[baseIdx + j]; }
		}
		ans += local_lgam_sum_param + local_sum_lgam_param + local_sum_par_lt;
	}
	return ans;
}

template<typename Int, typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param){
	assert(logTheta.size() == param.size());
	Real ans = 0.0;
	Real lgam_sum_param = lgamma(pmswitch::math::sum(param));
	Real sum_lgam_param = 0.0;
	Real sum_par_lt     = 0.0;
	{
		for(Int i = 0; i < param.size(); i++){ sum_lgam_param -= lgamma(param[i]);}
		for(Int i = 0; i < param.size(); i++){ sum_par_lt += (param[i]-1.0) * logTheta[i]; }
	}
	ans += lgam_sum_param + sum_lgam_param + sum_par_lt;
	return ans;
}

template<typename Int, typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param, Int dim,
								   const pmswitch::FixedSizeMultiVector<bool,Int> &filter){
	Int outN = param.size() / param.volume(dim);
	Int inN  = param.volume(dim);
	{// size, dim check in debug
		assert(logTheta.size() == param.size() && logTheta.dim() == param.dim());
		assert(filter.size() == param.size() && filter.dim() == param.dim());
		#ifdef NDEBUG
		#else
			for(Int d = 0; d < logTheta.dim(); d++){ assert(logTheta.size(d) == param.size(d)); }
			for(Int d = 0; d < filter.dim();   d++){ assert(filter.size(d)   == param.size(d)); }
		#endif
	}
	Real ans = 0.0;
	Int baseIdx = 0;
	for(Int i = 0; i < outN; i++){
		baseIdx = i * param.volume(dim);
		{/*
			std::cerr << "baseIdx: " << baseIdx << std::endl;
			for(Int j = 0; j < inN; j++){
				std::cerr << "j: " << j << ", param[baseIdx+j]: " << param[baseIdx+j] << ", filter[baseIdx+j]: " << filter[baseIdx+j] << std::endl;
			}
			std::cerr << pmswitch::math::sum(param, baseIdx, dim, filter) << std::endl;
		*/}
		Real local_lgam_sum_param = lgamma(pmswitch::math::sum(param, baseIdx, dim, filter));
		Real local_sum_lgam_param = 0.0;
		Real local_sum_par_lt     = 0.0;
		{
			for(Int j = 0; j < inN; j++)if(filter[baseIdx + j]){ local_sum_lgam_param -= lgamma(param[baseIdx + j]);}
			for(Int j = 0; j < inN; j++)if(filter[baseIdx + j]){ local_sum_par_lt     +=  (param[baseIdx + j]-1.0) * logTheta[baseIdx + j]; }
		}
		ans += local_lgam_sum_param + local_sum_lgam_param + local_sum_par_lt;
	}
	return ans;
}

template<typename Int, typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real,Int> &logTheta, const pmswitch::FixedSizeMultiVector<Real,Int> &param,
								   const pmswitch::FixedSizeMultiVector<bool,Int> &filter){
	assert(logTheta.size() == param.size());
	Real ans = 0.0;
	Real lgam_sum_param = lgamma(pmswitch::math::sum(param, filter));
	Real sum_lgam_param = 0.0;
	Real sum_par_lt     = 0.0;
	{
		for(Int i = 0; i < param.size(); i++)if(filter[i]){ sum_lgam_param -= lgamma(param[i]);}
		for(Int i = 0; i < param.size(); i++)if(filter[i]){ sum_par_lt     += (param[i]-1.0) * logTheta[i]; }
	}
	ans += lgam_sum_param + sum_lgam_param + sum_par_lt;
	return ans;
}


/*
		// log_theta[0] =  Eq[ log_Var[0] ]
    static inline double calELogDir(std::vector<double> log_theta, std::vector<double> param) {
        LOG(logDEBUG) << " start calELogDir " << std::endl;
        //E[log p(theta|param)] : Dir
        //- sum(lgamma(param)) + lgamma(sum(param)) + sum((param-1) * log_theta)
        double sum_lgam_param = 0.0;
        for (int i = 0; i < param.size(); i++) {
            sum_lgam_param -= lgamma(param[i]);
            LOG(logDEBUG) << " done lgamma " << std::endl;
        }
        double lgam_sum_param = lgamma(sumVec(param));
        double sum_par_lt = 0.0;
        for (int i = 0; i < param.size(); i++) {
            sum_par_lt += (param[i] - 1) * log_theta[i];
        }
        if (log_theta.size() != param.size()) {
            LOG(logDEBUG) << "--- calELogDir ---" << std::endl;
            printVec(log_theta);
            LOG(logDEBUG) << std::endl;
            printVec(param);
            LOG(logDEBUG) << std::endl;
            LOG(logDEBUG) << sum_lgam_param << " " << lgam_sum_param << " " << sum_par_lt << std::endl;
            throw std::string("calELogDir");
        }
        return sum_lgam_param + lgam_sum_param + sum_par_lt;
    }
*/


/*
################
template<typename Int = long long, typename Real = double>
class Inference{
public:
	Inference(FeatureData<Int, Real> _fvData, DBData<Int, Real> _dbData, PriorParameters<Int, Real> prior);
	~Inference();

	InferenceData<Int, Real> vb(Int T, bool updateDB = true);
	InferenceData<Int, Real> onlineVb();
private:
	FeatureData<Int, Real> fvData;
	DBData<Int, Real> dbData;
	PriorParameters<Int, Real> prior;
};
################
*/

/*
	public functions
*/
template<typename Int, typename Real>
pmswitch::Inference<Int, Real>::Inference(FeatureData<Int, Real> _fvData, DBData<Int, Real> _dbData, PriorParameters<Int, Real> _prior)
		:fvData(_fvData), dbData(_dbData), prior(_prior){
}

template<typename Int, typename Real>
pmswitch::InferenceData<Int, Real> pmswitch::Inference<Int, Real>::vb(Int T, bool updateDB ){
	using namespace pmswitch;
	Int N = dbData.N;
	Int I = fvData.I;
	FixedSizeMultiVector<Int> Js = fvData.Js;
	Int maxJ = fvData.maxJ;
	Int L = dbData.L;
	assert(dbData.L == fvData.L);
	Int maxMl = dbData.maxMl;
	FixedSizeMultiVector<Int> ml = fvData.Ml;
	FixedSizeMultiVector<Int> X  = fvData.X; // i,maxJ,l,maxMl

	FixedSizeMultiVector<Real, Int> EZ = math::initEZ<Int, Real>(I, Js, maxJ, T); // I, maxJ, T
	FixedSizeMultiVector<Real, Int> ES = math::initES<Int, Real>(I, Js, maxJ);    // I, maxJ, 2
	FixedSizeMultiVector<Real, Int> EY = math::initEY<Int, Real>(I, Js, maxJ, N); // I, maxJ, N

	FixedSizeMultiVector<Real, Int> alpha = prior.alpha; // I, T, 2
	FixedSizeMultiVector<Real, Int> beta  = prior.beta;  // I, 2
	FixedSizeMultiVector<Real, Int> gamma = prior.gamma; // I, N
	FixedSizeMultiVector<Real, Int> eta   = prior.eta;   // T, Ml.size(), maxMl

	FixedSizeMultiVector<Real, Int> dig_v     = math::calDirExp(alpha, 1); // I, T, 2
	FixedSizeMultiVector<Real, Int> dig_theta = math::calDirExp(beta,  0); // I, 2
	FixedSizeMultiVector<Real, Int> dig_pi    = math::calDirExp(gamma, 0); // I, N
	FixedSizeMultiVector<Real, Int> dig_f     = math::calDirExp(eta,   1); // T, L, maxMl
	FixedSizeMultiVector<Real, Int> g         = dbData.g; // N, l, maxMl

	FixedSizeMultiVector<Real, Int> dig_v_sig(0.0, I, T); // I, T

	FixedSizeMultiVector<bool, Int> EZFilter(true, I, maxJ, T);   math::makeFilter(EZFilter,  I, Js, FixedSizeMultiVector<Int, Int>(T, maxJ) );
	FixedSizeMultiVector<bool, Int> ESFilter(true, I, maxJ, 2);   math::makeFilter(ESFilter,  I, Js, FixedSizeMultiVector<Int, Int>(2, maxJ) );
	FixedSizeMultiVector<bool, Int> EYFilter(true, I, maxJ, N);   math::makeFilter(EYFilter,  I, Js, FixedSizeMultiVector<Int, Int>(N, maxJ) );
	FixedSizeMultiVector<bool, Int> etaFilter(true, T, L, maxMl); math::makeFilter(etaFilter, T, FixedSizeMultiVector<Int, Int>(T, L), ml);
	FixedSizeMultiVector<bool, Int> gFilter(true, I, maxJ, N);    math::makeFilter(gFilter,   N, FixedSizeMultiVector<Int, Int>(N, L), ml);

	FixedSizeMultiVector<Real, Int> lnEZ = math::applied(EZ, std::log, EZFilter) ; // I, maxJ, T
	FixedSizeMultiVector<Real, Int> lnES = math::applied(ES, std::log, ESFilter) ; // I, maxJ, 2
	FixedSizeMultiVector<Real, Int> lnEY = math::applied(EY, std::log, EYFilter) ; // I, maxJ, N

	Real minLqConvergence = 1e-6;
    Real beforeLq = -1.0 * 1e50;
    Real nextLq = -1.0 * 1e50;
    bool first = true;
    Int update = 0;

    Real INF = 1e300;

    while( first || std::abs(beforeLq - nextLq) >=  minLqConvergence  ){
	    if(first) first = false;
	    if( beforeLq > nextLq + minLqConvergence ){ std::cerr << "decrease lower bound @ vb" << std::endl; }
	    assert(beforeLq > nextLq + minLqConvergence);
	    // E-step
	    {// EZ
		    for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++)for(Int k = 0; k < T; k++) EZ(i,j,k) = 0.0;

		    for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T; k++){
	    		for(Int kk = 0; kk < k-1; kk++){ EZ(i,j,k) += dig_v(i,kk,1);}
	    		if(k < T) {EZ(i,j,k) += dig_v(i,k,0);}
	    		for(Int l = 0; l < L; l++)for(Int m = 0; m < ml(l); m++){
	    			EZ(i,j,k) += ES(i,j) * X(i,j,l,m) * dig_f(k, l, m);
	    		}
	    		math::filter(EZ, EZFilter, -INF);
	    		subByMax(EZ, 1);
	    		math::apply(EZ, std::exp);
	    		math::filter(EZ, EZFilter, 0.0);
	    		math::norm(EZ, 1);
	    		lnEZ = math::apply(EZ, std::log, EZFilter);
		    }
	    }
	    {// ES
	    	for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++){ES(i,j,0) = 0.0; ES(i,j,1) = 0.0;}

	    	for(Int i = 0; i < I; i++) for(Int j =0; j < maxJ; ){
	    		ES(i,j,0) += dig_theta(i, 0);
	    		for(Int l = 0; l < L; l++)for(Int k = 0; k < T; k++)for(Int m = 0; m < ml(l); m++){
	    			ES(i,j,0) += EZ(i,j,k) * X(i,j,l,m) * dig_f(k, l, m);
	    		}

	    		ES(i,j,1) += dig_theta(i,1);
	    		for(Int l = 0; l < L; l++)for(Int n = 0; n < N; n++)for(int m = 0; m < ml(l); m++){
	    			ES(i,j,1) += EY(i,j,n) * X(i,j,l,m) + std::log(g(n,l,m));
	    		}
	    		math::filter(ES, ESFilter, -INF);
	    		subByMax(ES, 1);
	    		math::apply(ES, std::exp);
	    		math::filter(ES, ESFilter, 0.0);
	    		math::norm(ES, 1);
	    		lnES = math::apply(ES, std::log, ESFilter);
	    	}
	    }
		{// EY // I, maxJ, N
	    	for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++)for(Int n = 0; n < N; n++){EY(i,j,n) = 0.0;}

	    	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int n = 0; n < N; n++){
	    		EY(i,j,n) += dig_pi(i,n);
	    		for(Int l = 0; l < L; l++)for(Int m = 0; m < ml(l); m++){
	    			EY(i,j,n) += X(i,j,l,m) * ES(i,j,1) * std::log(g(n,l,m));
	    		}
	    		math::filter(EY, EYFilter, -INF);
	    		subByMax(EY, 1);
	    		math::apply(EY, std::exp);
	    		math::filter(EY, EYFilter, 0.0);
	    		math::norm(EY, 1);
	    		lnEY = math::apply(EY, std::log, EYFilter);
	    	}
	    }
	    // M-step
	    {// FixedSizeMultiVector<Real, Int> eta   = prior.eta;   // T, Ml.size(), maxMl
	   		FixedSizeMultiVector<Real, Int> eta   = prior.eta;
	   		for(Int k = 0; k < T; k++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml(l); m++){
	   			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
		   			eta(k,l,m) += X(i,j,m,l) * EZ(i,j,k) * ES(i,j);
	   			}
	   		}
	   		dig_f = math::calDirExp(eta, 1, etaFilter);
	    }
	    {// FixedSizeMultiVector<Real, Int> alpha = prior.alpha; // I, T, 2
			FixedSizeMultiVector<Real, Int> alpha = prior.alpha; // I, T, 2
			for(Int i = 0; i < I; i++)for(Int k = 0; k < T; k++){
				for(Int j = 0; j < Js(i); j++){
					alpha(i,k,0)  += EZ(i,j,k);
					for(Int t = k+1; t < T; t++) alpha(i,k,1)  += EZ(i,j,t);
				}
			}
			dig_v = math::calDirExp(alpha, 1); // I, T, 2
	    }
	    {// FixedSizeMultiVector<Real, Int> beta  = prior.beta;  // I, 2
	    	FixedSizeMultiVector<Real, Int> beta  = prior.beta;  // I, 2
	    	for(Int i = 0; i < I; i++){
	    		for(Int j = 0; j < Js(i); j++){
		    		beta(i, 0) += ES(i,j,0); beta(i, 1) += ES(i,j,1);
	    		}
	    	}
	    	dig_theta = math::calDirExp(beta,  0);
	    }
		{// FixedSizeMultiVector<Real, Int> gamma = prior.gamma; // I, N
			FixedSizeMultiVector<Real, Int> gamma = prior.gamma; // I, N
			for(Int i = 0; i < I; i++)for(Int n = 0; n < N; n++){
				for(Int j = 0; j < Js(i); j++){
					gamma(i,n) += EY(i,j,n);
				}
	    	}
    	 	dig_pi    = math::calDirExp(gamma, 0); // I, N
		}

		if(updateDB){
			for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < maxMl; m++) g(n,l,m) = 0.0;

			for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml(l); m++){
				for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
					g(n,l,m) *= X(i,j,m,l) * EY(i,j,n) * ES(i,j,1);
				}
			}
    		math::filter(g, gFilter, 0.0);
    		math::norm(g, 2);
		}
	// FixedSizeMultiVector<Int> X  = fvData.X; // i,maxJ,l,maxMl
	// FixedSizeMultiVector<Real, Int> EZ = math::initEZ<Int, Real>(I, Js, maxJ, T); // I, maxJ, T
	// FixedSizeMultiVector<Real, Int> ES = math::initES<Int, Real>(I, Js, maxJ);    // I, maxJ, 2
	// FixedSizeMultiVector<Real, Int> EY = math::initEY<Int, Real>(I, Js, maxJ, N); // I, maxJ, N

	// FixedSizeMultiVector<Real, Int> alpha = prior.alpha; // I, T, 2
	// FixedSizeMultiVector<Real, Int> beta  = prior.beta;  // I, 2
	// FixedSizeMultiVector<Real, Int> gamma = prior.gamma; // I, N
	// FixedSizeMultiVector<Real, Int> eta   = prior.eta;   // T, Ml.size(), maxMl

	// FixedSizeMultiVector<Real, Int> dig_v     = math::calDirExp(alpha, 1); // I, T, 2
	// FixedSizeMultiVector<Real, Int> dig_theta = math::calDirExp(beta,  0); // I, 2
	// FixedSizeMultiVector<Real, Int> dig_pi    = math::calDirExp(gamma, 0); // I, N
	// FixedSizeMultiVector<Real, Int> dig_f     = math::calDirExp(eta,   1); // T, L, maxMl

	// FixedSizeMultiVector<Real, Int> g         = dbData.g; // N, l, maxMl
		{// nextLq = 0.0
			nextLq = 0.0;
			nextLq += math::calELogDir(dig_f,     prior.eta,    1, etaFilter);
			nextLq += math::calELogDir(dig_v,     prior.alpha,  1);
			nextLq += math::calELogDir(dig_theta, prior.beta,   0);
			nextLq += math::calELogDir(dig_pi,    prior.gamma,  0);

			{//dig_v_sig
				for(Int i = 0; i < I; i++)for(Int k = 0; k < T; k++){ dig_v_sig(i,k) = 0.0; }
				for(Int i = 0; i < I; i++)for(Int k = 0; k < T; k++){
					if(k < T-1){ dig_v_sig(i, k) += dig_v(i, k, 0);}
					for(Int t = 0; t < k-1; t++ ){ dig_v_sig(i, k) += dig_v(i, t, 1);}
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
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int l = 0; l < L; l++){
				for(Int k = 0; k < T; k++)for(Int m = 0; m < ml(l); m++){
					nextLq += EZ(i,j,k) * ES(i,j,0) * X(i,j,l,m) * dig_f(k,l,m);
				}
			}

			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int l = 0; l < L; l++){
				for(Int n = 0; n < N; n++)for(Int m = 0; m < ml(l); m++){
					nextLq += EY(i,j,n) * ES(i,j,1) * X(i,j,l,m) * std::log(g(n,l,m));
				}
			}
			nextLq -= math::calELogDir(dig_f,     eta,    1, etaFilter); // T, L, maxMl
			nextLq -= math::calELogDir(dig_v,     alpha,  1); 			 // I, T, 2
			nextLq -= math::calELogDir(dig_theta, beta,   0); 			 // I, 2
			nextLq -= math::calELogDir(dig_pi,    gamma,  0); 			 // I, N

			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				for(Int k = 0; k < T; k++) nextLq -= EZ(i,j,k) * lnEZ(i,j,k);
			}
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				nextLq -= (ES(i,j,0) * lnES(i,j,0) + ES(i,j,1) * lnES(i,j,1) );
			}
			for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
				for(Int n = 0; n < N; n++) nextLq -= EY(i,j,n) * lnEY(i,j,n);
			}
		}

	    beforeLq = nextLq;
	    nextLq = 0.0;
	    update++;
	}

	// InferenceData(FixedSizeMultiVector<double, Int> _EZ, FixedSizeMultiVector<double, Int> _ES, FixedSizeMultiVector<double, Int> _EY,
	// 			  FixedSizeMultiVector<double, Int> _alpha, FixedSizeMultiVector<double, Int> _eta,
	// 			  FixedSizeMultiVector<double, Int> _beta, FixedSizeMultiVector<double, Int> _gamma,
	// 			  FixedSizeMultiVector<double, Int> _g);
	pmswitch::InferenceData<Int, Real> ans(EZ, ES, EY, alpha, eta, beta, gamma, g);
	return ans;
}


/*
必要なデータ形式
	variable :: type //coment
	T     :: INT // truncation number
	N     :: INT // db size
	I     :: INT // sample size
	J_{i} :: INT // position size
	L     :: INT // dimension of feature vector
	M_{l} :: INT // # of possible state at l-th col of feature vector

	x_{i,j,l, }    :: 1 out of m_{l} vector // x_{i,j,l,m}

	//	k in [T]
	//	n in [N]
	//	i in [I]
	//		j in [J_{i}]
	//  l in [L]
	//		m in [M_{l}]


	EZ_{i,j, }     :: T-dimensional non-negative simplex // EZ_{i,j,k   }
	ES_{i, }	   :: 2-dimensional non-negative simplex // ES_{i,{0,1} }
	EY_{i, }       :: N-dimensional non-negative simplex // EY_{i,n}

###EZ_{i,j, }###
	v_{i, }        :: T-dimensional real in [0,1] //v_{i, k}
	dig_v_{i,}     :: T*2-dimensional real
	\alpha{i, }    :: T*2-dimensional real

	f_{k,l, }        :: m_{l}-dimensional non-negative simplex // f_{k,l}
	dig_f_{k,l, }    :: m_{l}-dimensional real vector          // dig_f_{k, l}
	\eta_{k,l, }	 :: m_{l}-dimensional real vector 		   // eta_{k, l}
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

/*
	private functions
*/

/*
################
template<typename Int = long long, typename Real = double>
class InferenceCreator{
public:
	InferenceCreator();
	pmswitch::Inference<Int, Real> createInference(std::string pathFV, std::string pathDB,
												   Int T, Real _alpha, Real beta0, Real beta1,
												   Real _gamma, Real _eta );
};
################
*/
/*
	public functions
*/
template<typename Int, typename Real>
pmswitch::InferenceCreator<Int, Real>::InferenceCreator(){

}

template<typename Int, typename Real>
pmswitch::Inference<Int, Real> pmswitch::InferenceCreator<Int, Real>::createInference(std::string pathFV, std::string pathDB,
													   Int T, Real _alpha, Real beta0, Real beta1,
													   Real _gamma, Real _eta ){
	using namespace pmswitch;
	InputFileParser<Int, Real> parser;
	InferenceCreator inferenceCreator;

	FeatureData<Int, Real> fvData = parser.parseFeatureFile(pathFV);
	DBData<Int, Real> dbData = parser.parseFeatureFile(pathDB);
	Inference<Int, Real> inference = inferenceCreator.createInference(fvData, dbData, T, _alpha, beta0, beta1, _gamma, _eta);
	return inference;
}

/*
	private functions
*/
#endif
