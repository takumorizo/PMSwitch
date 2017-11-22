#ifndef MathUtil_h
#define MathUtil_h

#include "FixedSizeMultiVector.hpp"

#include <vector>
#include <iostream>
#include <assert.h>
#include <string>
#include <boost/math/special_functions/digamma.hpp>
#include <math.h>
#include <cmath>
#include <random>

namespace pmswitch{
	namespace math{
		template<typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t baseIdx, pmswitch::size_t dim);

		template<typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real> &vec);

		template<typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t baseIdx, pmswitch::size_t dim,
				 const pmswitch::FixedSizeMultiVector<bool> &filter);

		template<typename Real = double>
		Real sum(const pmswitch::FixedSizeMultiVector<Real> &vec,
				 const pmswitch::FixedSizeMultiVector<bool> &filter);

		template<typename Real = double>
		Real dot(const pmswitch::FixedSizeMultiVector<Real> &a, const pmswitch::FixedSizeMultiVector<Real> &b);

		template<typename Real = double>
		Real dot(const pmswitch::FixedSizeMultiVector<Real> &a, const pmswitch::FixedSizeMultiVector<Real> &b,
				 const pmswitch::FixedSizeMultiVector<Real> &filter);

		template<typename Real = double>
		void subMax(pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t dim, Real INF = 1e300);

		template<typename Real = double>
		void subMin(pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t dim, Real INF = 1e300);

		template<typename Real = double>
		void norm(pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t dim);

		template<typename Real = double>
		void norm(pmswitch::FixedSizeMultiVector<Real> &vec);

		template<typename Real = double>
		void apply(pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real) );

		template<typename Real = double>
		pmswitch::FixedSizeMultiVector<Real> applied(const pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real) );

		template<typename Real = double>
		void apply(pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real),
				   const pmswitch::FixedSizeMultiVector<bool> &filter);

		template<typename Real = double>
		pmswitch::FixedSizeMultiVector<Real> applied(const pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real),
				   									   const pmswitch::FixedSizeMultiVector<bool> &filter);

		template<typename Real = double>
		void filter(pmswitch::FixedSizeMultiVector<Real> &vec, const pmswitch::FixedSizeMultiVector<bool> &filterVec, Real defaultVal);

		template<typename Int, class... Args>
		void makeFilter(pmswitch::FixedSizeMultiVector<bool> &vec, Int range, const Args&... args);

		namespace impl{
			template<typename Int, class... Args>
			void makeFilterImpl(pmswitch::FixedSizeMultiVector<bool> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int> &range, const Args&... args);

			template<typename Int, class... Args>
			void makeFilterImpl(pmswitch::FixedSizeMultiVector<bool> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int> &range);
		}

		template<typename Real = double>
		pmswitch::FixedSizeMultiVector<Real> calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha, pmswitch::size_t dim);

		template<typename Real = double>
		pmswitch::FixedSizeMultiVector<Real> calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha);

		template<typename Real = double>
		pmswitch::FixedSizeMultiVector<Real> calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha, pmswitch::size_t dim,
														   const pmswitch::FixedSizeMultiVector<bool> &filter);

		template<typename Real = double>
		pmswitch::FixedSizeMultiVector<Real> calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha,
														   const pmswitch::FixedSizeMultiVector<bool> &filter);

		template<typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param, pmswitch::size_t dim);

		template<typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param);

		template<typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param, pmswitch::size_t dim,
						const pmswitch::FixedSizeMultiVector<bool> &filter);

		template<typename Real = double>
		Real calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param,
						const pmswitch::FixedSizeMultiVector<bool> &filter);

		namespace random{
			template<typename URNG,typename Real>
			Real sampleGamma(URNG& engine, Real a, Real b);

			template<typename URNG, typename Real, typename Int = long long>
			std::vector<Real> sampleGamma(URNG& engine, Real a, Real b, Int num);

			template<typename URNG, typename Real>
			std::vector<Real> sampleDirichlet(URNG& engine, std::vector<Real> alpha);

			template<typename URNG, typename Real, typename Int = long long>
			pmswitch::FixedSizeMultiVector<Real> sampleDirichlet(URNG& engine, std::vector<Real> alpha, Int num);
		}
	}
}


template<typename URNG,typename Real>
Real pmswitch::math::random::sampleGamma(URNG& engine, Real a, Real b){
 	std::gamma_distribution<Real> distGamma(a,b);
 	Real p = distGamma(engine);
 	return p;
}

template<typename URNG, typename Real, typename Int>
std::vector<Real> pmswitch::math::random::sampleGamma(URNG& engine, Real a, Real b, Int num){
	std::vector<Real> ans;
	std::gamma_distribution<Real> distGamma(a,b);
	for(Int i = 0; i < num; i++) ans.push_back(distGamma(engine));
 	return ans;
}

template<typename URNG, typename Real>
std::vector<Real> pmswitch::math::random::sampleDirichlet(URNG& engine, std::vector<Real> alpha){
	std::vector<Real> ans;
	Real total = 0.0;
	for(pmswitch::size_t i = 0; i < alpha.size(); i++){
		std::gamma_distribution<Real> distGamma(alpha[i],1);
		ans.push_back(distGamma(engine));
		total += ans[i];
	}
	for(pmswitch::size_t i = 0; i < alpha.size(); i++) ans[i] /= total;
 	return ans;
}

template<typename URNG, typename Real, typename Int>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::random::sampleDirichlet(URNG& engine, std::vector<Real> alpha, Int num){
	pmswitch::FixedSizeMultiVector<Real> ans(0, num, alpha.size());
	for(Int n = 0; n < num; n++)for(Int i = 0; i < alpha.size(); i++){
		std::gamma_distribution<Real> distGamma(alpha[i],1);
		ans(n,i) = distGamma(engine);
	}
	pmswitch::math::norm<Int, Real>(ans, 0);
 	return ans;
}


template<typename Real>
void pmswitch::math::subMax(pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t dim, Real INF){
	pmswitch::size_t outN = vec.size() / vec.volume(dim);
	pmswitch::size_t inN  = vec.volume(dim);
	Real localMax = -INF;
	pmswitch::size_t baseIdx = 0;
	for(pmswitch::size_t i = 0; i < outN; i++){
		baseIdx = i * inN;
		localMax = -INF;
		for(pmswitch::size_t j = 0; j < inN; j++){ localMax = std::max(vec[baseIdx + j], localMax);}
		for(pmswitch::size_t j = 0; j < inN; j++){ vec[baseIdx + j] -= localMax;}
	}
}

template<typename Real>
void pmswitch::math::subMin(pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t dim, Real INF){
	pmswitch::size_t outN = vec.size() / vec.volume(dim);
	pmswitch::size_t inN  = vec.volume(dim);
	Real localMin = INF;
	pmswitch::size_t baseIdx = 0;
	for(pmswitch::size_t i = 0; i < outN; i++){
		baseIdx = i * inN;
		localMin = INF;
		for(pmswitch::size_t j = 0; j < inN; j++){ localMin = std::min(vec[baseIdx + j], localMin) ;}
		for(pmswitch::size_t j = 0; j < inN; j++){ vec[baseIdx + j] -= localMin;}
	}
}

template<typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t baseIdx, pmswitch::size_t dim){
	pmswitch::size_t inN = vec.volume(dim);
	Real sum = 0.0;
	for(pmswitch::size_t j = 0; j < inN; j++){ sum += vec[baseIdx + j];}
	return sum;
}

template<typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real> &vec){
	Real sum = 0.0;
	for(pmswitch::size_t i = 0; i < vec.size(); i++){ sum += vec[i];}
	return sum;
}

template<typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t baseIdx, pmswitch::size_t dim,
							const pmswitch::FixedSizeMultiVector<bool> &filter){
	assert(vec.size() == filter.size() && vec.dim() == filter.dim());
	#ifdef NDEBUG
	#else
		for(pmswitch::size_t d = 0; d < vec.dim(); d++){ assert(vec.size(d) == filter.size(d)); }
	#endif

	pmswitch::size_t inN = vec.volume(dim);
	Real sum = 0.0;
	for(pmswitch::size_t j = 0; j < inN; j++)if(filter[baseIdx + j]){ sum += vec[baseIdx + j];}
	return sum;
}

template<typename Real>
Real pmswitch::math::sum(const pmswitch::FixedSizeMultiVector<Real> &vec,
						 const pmswitch::FixedSizeMultiVector<bool> &filter){
	assert(vec.size() == filter.size() && vec.dim() == filter.dim());
	#ifdef NDEBUG
	#else
		for(pmswitch::size_t d = 0; d < vec.dim(); d++){ assert(vec.size(d) == filter.size(d)); }
	#endif

	Real sum = 0.0;
	for(pmswitch::size_t i = 0; i < vec.size(); i++){ sum += (pmswitch::size_t)(filter[i]) * vec[i];}
	return sum;
}


template<typename Real>
Real pmswitch::math::dot(const pmswitch::FixedSizeMultiVector<Real> &a, const pmswitch::FixedSizeMultiVector<Real> &b){
	assert(a.size() == b.size() && a.dim() == b.dim());
	#ifdef NDEBUG
	#else
		for(pmswitch::size_t d = 0; d < a.dim(); d++){ assert(a.size(d) == b.size(d)); }
	#endif
	Real dot = 0.0;
	for(pmswitch::size_t i = 0; i < a.size(); i++){ dot += a(i) * b(i);}
	return dot;
}

template<typename Real>
Real pmswitch::math::dot(const pmswitch::FixedSizeMultiVector<Real> &a, const pmswitch::FixedSizeMultiVector<Real> &b,
						 const pmswitch::FixedSizeMultiVector<Real> &filter){
	assert(a.size() == b.size()      && a.dim() == b.dim());
	assert(a.size() == filter.size() && a.dim() == filter.dim());
	#ifdef NDEBUG
	#else
		for(pmswitch::size_t d = 0; d < a.dim(); d++){ assert(a.size(d) == b.size(d)); }
		for(pmswitch::size_t d = 0; d < a.dim(); d++){ assert(a.size(d) == filter.size(d)); }
	#endif
	Real dot = 0.0;
	for(pmswitch::size_t i = 0; i < a.size(); i++){ dot += (pmswitch::size_t)(filter[i]) * a[i] * b[i];}
	return dot;
}

template<typename Real>
void pmswitch::math::norm(pmswitch::FixedSizeMultiVector<Real> &vec, pmswitch::size_t dim){
	pmswitch::size_t outN = vec.size() / vec.volume(dim);
	pmswitch::size_t inN  = vec.volume(dim);
	Real sum = 0.0;
	pmswitch::size_t baseIdx = 0;
	for(pmswitch::size_t i = 0; i < outN; i++){
		baseIdx = i * inN;
		sum = 0.0;
		for(pmswitch::size_t j = 0; j < inN; j++){ sum += vec[baseIdx + j];}
		if(sum != 0)for(pmswitch::size_t j = 0; j < inN; j++){ vec[baseIdx + j] /= sum;}
	}
}

template<typename Real>
void pmswitch::math::norm(pmswitch::FixedSizeMultiVector<Real> &vec){
	Real sum = 0.0;
	for(pmswitch::size_t i = 0; i < vec.size(); i++){ sum += vec[i];}
	if(sum != 0)for(pmswitch::size_t i = 0; i < vec.size(); i++){ vec[i] /= sum;}
}

template<typename Real>
void pmswitch::math::apply(pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real) ){
	for(pmswitch::size_t i = 0; i < vec.size(); i++){ vec[i] = func(vec[i]);}
}

template<typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::applied(const pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real) ){
	pmswitch::FixedSizeMultiVector<Real> ans = vec;
	apply(ans, func);
	return ans;
}

template<typename Real>
void pmswitch::math::apply(pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real),
						   const FixedSizeMultiVector<bool> &filter){
	{// size, dim check in debug
		assert(vec.size() == filter.size() && vec.dim() == filter.dim());
		#ifdef NDEBUG
		#else
			for(pmswitch::size_t d = 0; d < vec.dim(); d++){ assert(vec.size(d) == filter.size(d)); }
		#endif
	}
	for(pmswitch::size_t i = 0; i < vec.size(); i++)if(filter[i]){ vec[i] = func(vec[i]);}
}

template<typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::applied(const pmswitch::FixedSizeMultiVector<Real> &vec, Real (*func)(Real),
							  									  const FixedSizeMultiVector<bool> &filter){
	pmswitch::FixedSizeMultiVector<Real> ans = vec;
	apply(ans, func, filter);
	return ans;
}

template<typename Real>
void pmswitch::math::filter(pmswitch::FixedSizeMultiVector<Real> &vec, const FixedSizeMultiVector<bool> &filterVec, Real defaultVal){
	assert(vec.size() == filterVec.size());
	for(pmswitch::size_t i = 0; i < vec.size(); i++)if( not filterVec[i] ){ vec[i] = defaultVal; }
}


template<typename Int, class... Args>
void pmswitch::math::makeFilter(pmswitch::FixedSizeMultiVector<bool> &vec, Int range, const Args&... args){
	for(Int i = 0; i < vec.size(0); i++){
		pmswitch::math::impl::makeFilterImpl(vec, (Int)1, (Int)(i*vec.volume((pmswitch::size_t)0)), (Int)i, (i < range), args...);
	}
}

template<typename Int, class... Args>
void pmswitch::math::impl::makeFilterImpl(pmswitch::FixedSizeMultiVector<bool> &vec, Int n, Int idx, Int prevIdx,  bool val,
	  									  const pmswitch::FixedSizeMultiVector<Int> &range, const Args&... args){
	for(Int i = 0; i < vec.size(n); i++){
		pmswitch::math::impl::makeFilterImpl(vec, n+1, (Int)(idx+i*vec.volume((pmswitch::size_t)n)), i, val && (i < range(prevIdx)), args...);
	}
}

template<typename Int, class... Args>
void pmswitch::math::impl::makeFilterImpl(pmswitch::FixedSizeMultiVector<bool> &vec, Int n, Int idx, Int prevIdx, bool val,
										  const pmswitch::FixedSizeMultiVector<Int> &range){
	for(Int i = 0; i < vec.size(n); i++){
		vec[(pmswitch::size_t)(idx+i*vec.volume(n))] = ( val && ( i < range(prevIdx) ) );
	}
}

template<typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha, pmswitch::size_t dim){

	pmswitch::FixedSizeMultiVector<Real> ansDig = alpha;
	pmswitch::size_t baseIdx = 0;
	Real sum = 0.0, digSum = 0.0;
	pmswitch::size_t outN = alpha.size() / alpha.volume(dim);
	pmswitch::size_t inN	= alpha.volume(dim);
	for(pmswitch::size_t i = 0; i < outN; i++){
		sum = 0.0;
		baseIdx = i * inN;
		for(pmswitch::size_t j = 0; j < inN; j++){ sum += alpha[baseIdx + j];}
		digSum =  boost::math::digamma(sum);
		for(pmswitch::size_t j = 0; j < inN; j++){
			ansDig[baseIdx + j] = boost::math::digamma(ansDig[baseIdx + j]) - digSum;
		}
	}
	return ansDig;
}

template<typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha){
	pmswitch::FixedSizeMultiVector<Real> ansDig = alpha;
	Real sum = 0.0, digSum = 0.0;
	sum = 0.0;
	for(pmswitch::size_t j = 0; j < alpha.size(); j++){ sum += alpha[j];}
	digSum =  boost::math::digamma(sum);
	for(pmswitch::size_t j = 0; j < alpha.size(); j++){ ansDig[j] = boost::math::digamma(ansDig[j]) - digSum;}
	return ansDig;
}

template<typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha, pmswitch::size_t dim,
															   const pmswitch::FixedSizeMultiVector<bool> &filter){

	pmswitch::FixedSizeMultiVector<Real> ansDig = alpha;
	pmswitch::size_t baseIdx = 0;
	Real sum = 0.0, digSum = 0.0;
	pmswitch::size_t outN = alpha.size() / alpha.volume(dim);
	pmswitch::size_t inN	= alpha.volume(dim);
	for(pmswitch::size_t i = 0; i < outN; i++){
		sum = 0.0;
		baseIdx = i * inN;
		for(pmswitch::size_t j = 0; j < inN; j++)if(filter[baseIdx + j]){ sum += alpha[baseIdx + j]; }
		digSum =  boost::math::digamma(sum);
		for(pmswitch::size_t j = 0; j < inN; j++)if(filter[baseIdx + j]){
			ansDig[baseIdx + j] =  boost::math::digamma(ansDig[baseIdx + j]) - digSum;
		}
	}
	return ansDig;
}

template<typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::math::calDirExp(const pmswitch::FixedSizeMultiVector<Real> &alpha, const pmswitch::FixedSizeMultiVector<bool> &filter){
	pmswitch::FixedSizeMultiVector<Real> ansDig = alpha;
	Real sum = 0.0, digSum = 0.0;
	sum = 0.0;
	for(pmswitch::size_t j = 0; j < alpha.size(); j++)if(filter[j]){ sum += alpha[j];}
	digSum =  boost::math::digamma(sum);
	for(pmswitch::size_t j = 0; j < alpha.size(); j++)if(filter[j]){ ansDig[j] = boost::math::digamma(ansDig[j]) - digSum;}
	return ansDig;
}

template<typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param, pmswitch::size_t dim){
	pmswitch::size_t outN = param.size() / param.volume(dim);
	pmswitch::size_t inN  = param.volume(dim);
	{// size, dim check in debug
		assert(logTheta.size() == param.size() && logTheta.dim() == param.dim());
		#ifdef NDEBUG
		#else
			for(pmswitch::size_t d = 0; d < logTheta.dim(); d++){ assert(logTheta.size(d) == param.size(d)); }
		#endif
	}
	Real ans = 0.0;
	pmswitch::size_t baseIdx = 0;
	for(pmswitch::size_t i = 0; i < outN; i++){
		baseIdx = i * param.volume(dim);
		Real local_lgam_sum_param = lgamma(pmswitch::math::sum(param, baseIdx, dim));
		Real local_sum_lgam_param = 0.0;
		Real local_sum_par_lt     = 0.0;
		{
			for(pmswitch::size_t j = 0; j < inN; j++){ local_sum_lgam_param -= lgamma(param[baseIdx + j]);}
			for(pmswitch::size_t j = 0; j < inN; j++){ local_sum_par_lt += (param[baseIdx + j]-1.0) * logTheta[baseIdx + j]; }
		}
		ans += local_lgam_sum_param + local_sum_lgam_param + local_sum_par_lt;
	}
	return ans;
}

template<typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param){
	assert(logTheta.size() == param.size());
	Real ans = 0.0;
	Real lgam_sum_param = lgamma(pmswitch::math::sum(param));
	Real sum_lgam_param = 0.0;
	Real sum_par_lt     = 0.0;
	{
		for(pmswitch::size_t i = 0; i < param.size(); i++){ sum_lgam_param -= lgamma(param[i]);}
		for(pmswitch::size_t i = 0; i < param.size(); i++){ sum_par_lt += (param[i]-1.0) * logTheta[i]; }
	}
	ans += lgam_sum_param + sum_lgam_param + sum_par_lt;
	return ans;
}

template<typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param, pmswitch::size_t dim,
								const pmswitch::FixedSizeMultiVector<bool> &filter){
	pmswitch::size_t outN = param.size() / param.volume(dim);
	pmswitch::size_t inN  = param.volume(dim);
	{// size, dim check in debug
		assert(logTheta.size() == param.size() && logTheta.dim() == param.dim());
		assert(filter.size() == param.size() && filter.dim() == param.dim());
		#ifdef NDEBUG
		#else
			for(pmswitch::size_t d = 0; d < logTheta.dim(); d++){ assert(logTheta.size(d) == param.size(d)); }
			for(pmswitch::size_t d = 0; d < filter.dim();   d++){ assert(filter.size(d)   == param.size(d)); }
		#endif
	}
	Real ans = 0.0;
	pmswitch::size_t baseIdx = 0;
	for(pmswitch::size_t i = 0; i < outN; i++){
		baseIdx = i * param.volume(dim);
		{/*
			std::cerr << "baseIdx: " << baseIdx << std::endl;
			for(pmswitch::size_t j = 0; j < inN; j++){
				std::cerr << "j: " << j << ", param[baseIdx+j]: " << param[baseIdx+j] << ", filter[baseIdx+j]: " << filter[baseIdx+j] << std::endl;
			}
			std::cerr << pmswitch::math::sum(param, baseIdx, dim, filter) << std::endl;
		*/
			// Real tmp = pmswitch::math::sum(param, baseIdx, dim, filter);
			// std::cerr << baseIdx << " : sum : " << tmp << std::endl;
		}
		Real local_lgam_sum_param = lgamma(pmswitch::math::sum(param, baseIdx, dim, filter));
		Real local_sum_lgam_param = 0.0;
		Real local_sum_par_lt     = 0.0;
		{
			for(pmswitch::size_t j = 0; j < inN; j++)if(filter[baseIdx + j]){ local_sum_lgam_param -= lgamma(param[baseIdx + j]); }
			for(pmswitch::size_t j = 0; j < inN; j++)if(filter[baseIdx + j]){ local_sum_par_lt     +=  (param[baseIdx + j]-1.0) * logTheta[baseIdx + j]; }
		}
		ans += local_lgam_sum_param + local_sum_lgam_param + local_sum_par_lt;
	}
	return ans;
}

template<typename Real>
Real pmswitch::math::calELogDir(const pmswitch::FixedSizeMultiVector<Real> &logTheta, const pmswitch::FixedSizeMultiVector<Real> &param,
								const pmswitch::FixedSizeMultiVector<bool> &filter){
	assert(logTheta.size() == param.size());
	Real ans = 0.0;
	Real lgam_sum_param = lgamma(pmswitch::math::sum(param, filter));
	Real sum_lgam_param = 0.0;
	Real sum_par_lt     = 0.0;
	{
		for(pmswitch::size_t i = 0; i < param.size(); i++)if(filter[i]){ sum_lgam_param -= lgamma(param[i]);}
		for(pmswitch::size_t i = 0; i < param.size(); i++)if(filter[i]){ sum_par_lt     += (param[i]-1.0) * logTheta[i]; }
	}
	ans += lgam_sum_param + sum_lgam_param + sum_par_lt;
	return ans;
}

#endif
