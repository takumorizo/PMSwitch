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

namespace pmswitch{
	namespace math{
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
		void subMin(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim, Real INF = 1e300);

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

		namespace impl{
			template<typename Int = long long, class... Args>
			void makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int, Int> &range, const Args&... args);

			template<typename Int = long long, class... Args>
			void makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int, Int> &range);
		}

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

template<typename Int, typename Real>
void pmswitch::math::subMax(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim, Real INF){
	Int outN = vec.size() / vec.volume(dim);
	Int inN  = vec.volume(dim);
	Real localMax = -INF;
	Int baseIdx = 0;
	for(Int i = 0; i < outN; i++){
		baseIdx = i * inN;
		localMax = -INF;
		for(Int j = 0; j < inN; j++){ localMax = std::max(vec[baseIdx + j], localMax);}
		for(Int j = 0; j < inN; j++){ vec[baseIdx + j] -= localMax;}
	}
}

template<typename Int, typename Real>
void pmswitch::math::subMin(pmswitch::FixedSizeMultiVector<Real,Int> &vec, Int dim, Real INF){
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
		if(sum != 0)for(Int j = 0; j < inN; j++){ vec[baseIdx + j] /= sum;}
	}
}

template<typename Int, typename Real>
void pmswitch::math::norm(pmswitch::FixedSizeMultiVector<Real,Int> &vec){
	Real sum = 0.0;
	for(Int i = 0; i < vec.size(); i++){ sum += vec[i];}
	if(sum != 0)for(Int i = 0; i < vec.size(); i++){ vec[i] /= sum;}
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
		pmswitch::math::impl::makeFilterImpl(vec, (long long)1, i*vec.volume(0), i, (i < range), args...);
	}
}

template<typename Int, class... Args>
void pmswitch::math::impl::makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx,  bool val,
  									   const pmswitch::FixedSizeMultiVector<Int, Int> &range, const Args&... args){
	for(Int i = 0; i < vec.size(n); i++){
		pmswitch::math::impl::makeFilterImpl(vec, n+1, idx+i*vec.volume(n), i, val && (i < range(prevIdx)), args...);
	}
}

template<typename Int, class... Args>
void pmswitch::math::impl::makeFilterImpl(pmswitch::FixedSizeMultiVector<bool, Int> &vec, Int n, Int idx, Int prevIdx, bool val, const pmswitch::FixedSizeMultiVector<Int, Int> &range){
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
		*/
			// Real tmp = pmswitch::math::sum(param, baseIdx, dim, filter);
			// std::cerr << baseIdx << " : sum : " << tmp << std::endl;
		}
		Real local_lgam_sum_param = lgamma(pmswitch::math::sum(param, baseIdx, dim, filter));
		Real local_sum_lgam_param = 0.0;
		Real local_sum_par_lt     = 0.0;
		{
			for(Int j = 0; j < inN; j++)if(filter[baseIdx + j]){ local_sum_lgam_param -= lgamma(param[baseIdx + j]); }
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

#endif
