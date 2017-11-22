#ifndef Inference_h
#define Inference_h

#include "ExceptionUtil.hpp"
#include "DBData.hpp"
#include "FeatureData.hpp"
#include "InferenceData.hpp"
#include "PriorParameters.hpp"
#include "InputFileParser.hpp"
#include "MathUtil.hpp"
#include "Parameters.hpp"


#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <string>
#include <boost/math/special_functions/digamma.hpp>
#include <math.h>
#include <cmath>
#include <random>

namespace pmswitch{
	template<typename Int = long long, typename Real = double>
	class Inference{
	public:
		Inference(FeatureData<Int, Real> _fvData, DBData<Int, Real> _dbData, PriorParameters<Int, Real> prior, Int _T);
		InferenceData<Int, Real> vb(bool updateDB = true, const Real epsilon = 1e-9, std::string outputErrLatentDir = "", std::string inputLatentDir = "") const;
		InferenceData<Int, Real> vbFull(bool updateDB = true, const Real epsilon = 1e-9, std::string outputErrLatentDir = "", std::string inputLatentDir = "") const;
		InferenceData<Int, Real> onlineVb();
	private:
		FeatureData<Int, Real> fvData;
		DBData<Int, Real> dbData;
		PriorParameters<Int, Real> prior;
		Int T;

		FixedSizeMultiVector<Real> initEZ(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int T, std::string inputPath = "") const;
		FixedSizeMultiVector<Real> initES(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, 		std::string inputPath = "") const;
		FixedSizeMultiVector<Real> initEY(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int N, std::string inputPath = "") const;

		void updateEZ(FixedSizeMultiVector<Real> &EZ, FixedSizeMultiVector<Real> &lnEZ, const FixedSizeMultiVector<bool> &EZFilter,
					  const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml, const Int &N,
			          const FixedSizeMultiVector<Real> &dig_f, const FixedSizeMultiVector<Real> &ES,
			          const FixedSizeMultiVector<Int>  &X,     const FixedSizeMultiVector<Real> &dig_v) const;
		void updateES(FixedSizeMultiVector<Real> &ES, FixedSizeMultiVector<Real> &lnES, const FixedSizeMultiVector<bool> &ESFilter,
					  const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml, const Int &N,
					  const FixedSizeMultiVector<Real> &dig_f, const FixedSizeMultiVector<Real> &EZ,
					  const FixedSizeMultiVector<Real> &lng,   const FixedSizeMultiVector<Real> &EY,
					  const FixedSizeMultiVector<Int>  &X,     const FixedSizeMultiVector<Real> &dig_theta) const;

		void updateEY(FixedSizeMultiVector<Real> &EY, FixedSizeMultiVector<Real> &lnEY, const FixedSizeMultiVector<bool> &EYFilter,
					  const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &L, const FixedSizeMultiVector<Int> &ml, const Int &N,
					  const FixedSizeMultiVector<Real> &lng,   const FixedSizeMultiVector<Real> &ES,
					  const FixedSizeMultiVector<Int>  &X,     const FixedSizeMultiVector<Real> &dig_pi) const;

		void udpateF(FixedSizeMultiVector<Real> &eta, const FixedSizeMultiVector<bool> &etaFilter, FixedSizeMultiVector<Real> &dig_f,
					const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml,
					const PriorParameters<Int, Real> &prior,   const FixedSizeMultiVector<Int>  &X,
					const FixedSizeMultiVector<Real> &EZ, const FixedSizeMultiVector<Real> &ES) const;

		void updateTheta(FixedSizeMultiVector<Real> &beta, FixedSizeMultiVector<Real> &dig_theta,
						 const Int &I, const FixedSizeMultiVector<Int> &Js,
						 const PriorParameters<Int, Real> &prior, const FixedSizeMultiVector<Real> &ES) const;

		void updateThetaWithBetaDistributed(FixedSizeMultiVector<Real> &beta, FixedSizeMultiVector<Real> &dig_theta,
											const Int &I, const FixedSizeMultiVector<Int> &Js,
											const FixedSizeMultiVector<Real> &avgBeta, const FixedSizeMultiVector<Real> &ES) const;

		void updatePi(FixedSizeMultiVector<Real> &gamma, FixedSizeMultiVector<Real> &dig_pi,
					  const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &N,
					  const PriorParameters<Int, Real> &prior, const FixedSizeMultiVector<Real> &EY) const;

		void updateV(FixedSizeMultiVector<Real> &alpha, FixedSizeMultiVector<Real> &dig_v,
			         const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &T,
			         const PriorParameters<Int, Real> &prior, const FixedSizeMultiVector<Real> &EZ) const;

		void updateVWithAlphaDistributed(FixedSizeMultiVector<Real> &alpha, FixedSizeMultiVector<Real> &dig_v,
								         const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &T,
								         const FixedSizeMultiVector<Real> &avgAlpha, const FixedSizeMultiVector<Real> &EZ) const;

		void updateAlpha() const; 		                //updateAlpha

		void updateThetaWithBetaDistributed() const;

		void updateBeta() const;

		void updateG(FixedSizeMultiVector<Real> &g, FixedSizeMultiVector<Real> &lng, const FixedSizeMultiVector<bool> &gFilter,
					 const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml, Int maxMl, const Int &N,
					 const FixedSizeMultiVector<Real> &EY, const FixedSizeMultiVector<Real> &ES,
					 const FixedSizeMultiVector<Int> &X, Real eps = 1e-300) const;

		FixedSizeMultiVector<Real> makeAvgAlpha(const FixedSizeMultiVector<Real> &s, Int I, Int T) const;
		FixedSizeMultiVector<Real> makeAvgBeta(const FixedSizeMultiVector<Real> &u, Int I) const;

		void makeAvgAlpha(FixedSizeMultiVector<Real> &ans, const FixedSizeMultiVector<Real> &s, Int I, Int T) const;
		void makeAvgBeta(FixedSizeMultiVector<Real> &ans,  const FixedSizeMultiVector<Real> &u, Int I) const;

		Real lnProbGamma(const FixedSizeMultiVector<Real> &param, const FixedSizeMultiVector<Real> &trueParam) const;
	};

	template<typename Int = long long, typename Real = double>
	class InferenceCreator{
	public:
		InferenceCreator();
	  	static pmswitch::Inference<Int, Real> createInference(std::string pathFV, std::string pathDB,
													   Int T, Real _alpha, Real _beta0, Real _beta1,
													   Real _gamma, Real _eta,
													   Real s0, Real s1,
													   Real u0, Real u1,
													   std::string alphaPath = "", std::string betaPath = "",
													   std::string gammaPath = "", std::string etaPath  = "",
													   std::string sPath     = "", std::string uPath     = "",
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
Real pmswitch::Inference<Int, Real>::lnProbGamma(const FixedSizeMultiVector<Real> &param, const FixedSizeMultiVector<Real> &trueParam) const{
	// E[ -ln x - s1 * x + s0 * ln_x - lnGamma(s0) + s0 * ln_s1 ]
	// E[ x ] = s0 / s1
	// E[ ln_x ] = diggamma(s0) - log(s1)
	Real E_ln_x = boost::math::digamma(trueParam(0)) - std::log(trueParam(1));
	Real E_x    = trueParam(0) / trueParam(1);
	return  -1.0*E_ln_x - param(1)*E_x + param(0)*E_ln_x - std::lgamma(param(0)) + param(0)*std::log(param(1));
}

template<typename Int, typename Real>
pmswitch::InferenceData<Int, Real> pmswitch::Inference<Int, Real>::vbFull(bool updateDB, const Real epsilon,
																		  std::string outputErrLatentDir, std::string inputLatentDir) const{
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

	FixedSizeMultiVector<Real> firstEZ = initEZ(I, Js, maxJ, T, inputLatentDir + "/EZ.txt"); // I, maxJ, T
	FixedSizeMultiVector<Real> firstES = initES(I, Js, maxJ   , inputLatentDir + "/ES.txt");    // I, maxJ, 2
	FixedSizeMultiVector<Real> firstEY = initEY(I, Js, maxJ, N, inputLatentDir + "/EY.txt"); // I, maxJ, N

	FixedSizeMultiVector<Real> EZ = firstEZ;
	FixedSizeMultiVector<Real> ES = firstES;
	FixedSizeMultiVector<Real> EY = firstEY;

	FixedSizeMultiVector<Real> S     = prior.s;     // 2
	FixedSizeMultiVector<Real> U     = prior.u;     // 2
	FixedSizeMultiVector<Real> avg_alpha = makeAvgAlpha(S, I, T); // I, T-1, 2
	FixedSizeMultiVector<Real> avg_beta  = makeAvgBeta(U, I); // I, 2

	FixedSizeMultiVector<Real> gamma = prior.gamma; // I, N
	FixedSizeMultiVector<Real> eta   = prior.eta;   // T, Ml.size(), maxMl

	// FixedSizeMultiVector<Real> prior_alpha = makeAvgAlpha(S, I, T);
	// FixedSizeMultiVector<Real> prior_beta  = makeAvgBeta(U, I);
	FixedSizeMultiVector<Real> alpha       = avg_alpha;
	FixedSizeMultiVector<Real> beta        = avg_beta;
	FixedSizeMultiVector<Real> dig_v       = math::calDirExp<Real>(avg_alpha, 1); // I, T-1, 2
	FixedSizeMultiVector<Real> dig_theta   = math::calDirExp<Real>(avg_beta,  0); // I, 2
	FixedSizeMultiVector<Real> dig_pi      = math::calDirExp<Real>(gamma, 0); // I, N
	FixedSizeMultiVector<Real> dig_f       = math::calDirExp<Real>(eta,   1); // T, L, maxMl
	FixedSizeMultiVector<Real> g           = dbData.g; // N, l, maxMl
	FixedSizeMultiVector<Real> dig_v_sig(0.0, I, T); // I, T

	FixedSizeMultiVector<bool> EZFilter(true, I, maxJ, T);   math::makeFilter(EZFilter,  I, Js, FixedSizeMultiVector<Int>(T, maxJ) );
	FixedSizeMultiVector<bool> ESFilter(true, I, maxJ, 2);   math::makeFilter(ESFilter,  I, Js, FixedSizeMultiVector<Int>(2, maxJ) );
	FixedSizeMultiVector<bool> EYFilter(true, I, maxJ, N);   math::makeFilter(EYFilter,  I, Js, FixedSizeMultiVector<Int>(N, maxJ) );
	FixedSizeMultiVector<bool> etaFilter(true, T, L, maxMl); math::makeFilter(etaFilter, T, FixedSizeMultiVector<Int>(L, T), ml);
	FixedSizeMultiVector<bool> gFilter(true, N, L, maxMl);   math::makeFilter(gFilter,   N, FixedSizeMultiVector<Int>(L, N), ml);

	FixedSizeMultiVector<Real> lnEZ = math::applied<Real>(EZ, std::log, EZFilter) ; // I, maxJ, T
	FixedSizeMultiVector<Real> lnES = math::applied<Real>(ES, std::log, ESFilter) ; // I, maxJ, 2
	FixedSizeMultiVector<Real> lnEY = math::applied<Real>(EY, std::log, EYFilter) ; // I, maxJ, N
	FixedSizeMultiVector<Real> lng  = math::applied<Real>(g,  std::log,  gFilter) ; // N, L, maxM

	// {   // init g
	// 	Real eps = 1e-300;
	// 	math::filter<Int, Real>(g, gFilter, 0);
	// 	math::norm<Int, Real>(g, 1);
	// 	for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++)if(g(n,l,m) < eps) g(n,l,m) = eps;
	// }

    Real beforeLq = -1.0 * 1e50;
    Real nextLq = -1.0 * 1e50;
    bool first = true;
    Int update = 0;

    Real INF = 1e300;

    while( first || std::abs(nextLq - beforeLq)/std::abs(beforeLq) >=  epsilon){
	    beforeLq = nextLq;
	    nextLq = 0.0;

	    if(updateDB){
			updateG(g, lng, gFilter, I, Js, maxJ, T, L, ml, maxMl, N, EY, ES, X);
		}

	    // E-step
    	updateEZ(EZ, lnEZ, EZFilter, I, Js, maxJ, T, L, ml, N, dig_f, ES, X, dig_v);
     	updateES(ES, lnES, ESFilter, I, Js, maxJ, T, L, ml, N, dig_f, EZ, lng, EY, X, dig_theta);
		updateEY(EY, lnEY, EYFilter, I, Js, maxJ, L, ml, N, lng, ES, X, dig_pi);

	    // M-step
		udpateF(eta, etaFilter, dig_f, I, Js, T, L, ml, prior, X, EZ, ES);
		updatePi(gamma, dig_pi, I, Js, N, prior, EY);
		updateVWithAlphaDistributed(alpha, dig_v, I, Js, T, avg_alpha, EZ);
		updateThetaWithBetaDistributed(beta, dig_theta, I, Js, avg_beta, ES);
		// updateS, updateU
	    {//updateS
	    	S(0) = prior.s(0) + ( I * (T-1) );
	    	S(1) = prior.s(1);
	    	for(Int i = 0; i < I; i++)for(Int k = 0; k < T-1; k++) S(1) -= dig_v(i,k,1);
    		makeAvgAlpha(avg_alpha, S, I, T); // updateAvgAlpha
	    }
	    {//updateU
	    	U(0) = prior.u(0) + I ;
	    	U(1) = prior.u(1);
	    	for(Int i = 0; i < I; i++) U(1) -= dig_theta(i,1);
		    makeAvgBeta(avg_beta, U, I); //updateAvgBeta
	    }

		{// nextLq = 0.0
			nextLq = 0.0;
			nextLq += lnProbGamma(prior.s, S);
			nextLq += lnProbGamma(prior.u, U);
			nextLq += math::calELogDir<Real>(dig_f,     prior.eta,    1, etaFilter); // T, L, maxMl
			nextLq += math::calELogDir<Real>(dig_v,     avg_alpha,    1); // I, T-1, 2
			nextLq += math::calELogDir<Real>(dig_theta, avg_beta,     0); // I, 2
			nextLq += math::calELogDir<Real>(dig_pi,    prior.gamma,  0); // I, N

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

			nextLq -= lnProbGamma(S, S);
			nextLq -= lnProbGamma(U, U);
			nextLq -= math::calELogDir<Real>(dig_f,     eta,    1, etaFilter);  // T, L, maxMl
			nextLq -= math::calELogDir<Real>(dig_v,     alpha,  1); 			 // I, T-1, 2
			nextLq -= math::calELogDir<Real>(dig_theta, beta,   0); 			 // I, 2
			nextLq -= math::calELogDir<Real>(dig_pi,    gamma,  0); 			 // I, N

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
	    	std::cerr << "!decrease lower bound @ vb, update : " << update << std::endl;
	    	std::cerr << "vb, update : " << update << ", beforeLq :" << beforeLq << " nextLq : "  << nextLq << ", nextLq-beforeLq/beforeLq :" << (nextLq - beforeLq)/std::abs(beforeLq) << std::endl;

	    	std::ofstream ofs(outputErrLatentDir + "/EZ.txt");
		   	if( !ofs.fail()) {firstEZ.print(outputErrLatentDir + "/EZ.txt");}
			else{ std::cerr << "failed to open EZ.txt \n";}
		   	ofs.close();

	    	ofs.open(outputErrLatentDir + "/ES.txt");
		   	if( !ofs.fail()) {firstES.print(outputErrLatentDir + "/ES.txt");}
		   	else{ std::cerr << "failed to open ES.txt \n";}
		   	ofs.close();

	    	ofs.open(outputErrLatentDir + "/EY.txt");
		   	if( !ofs.fail()) {firstEY.print(outputErrLatentDir + "/EY.txt");}
		   	else{ std::cerr << "failed to open EY.txt \n";}

		   	ofs.open(outputErrLatentDir + "/EZnew.txt");
		   	if( !ofs.fail()) {EZ.print(outputErrLatentDir + "/EZnew.txt");}
		   	else{ std::cerr << "failed to open EZnew.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/ESnew.txt");
		   	if( !ofs.fail()) {ES.print(outputErrLatentDir + "/ESnew.txt");}
		   	else{ std::cerr << "failed to open ESnew.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/EYnew.txt");
		   	if( !ofs.fail()) {EY.print(outputErrLatentDir + "/EYnew.txt");}
		   	else{ std::cerr << "failed to open EYnew.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/alpha.txt");
		   	if( !ofs.fail()) {alpha.print(outputErrLatentDir + "/alpha.txt");}
		   	else{ std::cerr << "failed to open alpha.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/beta.txt");
		   	if( !ofs.fail()) {beta.print(outputErrLatentDir + "/beta.txt");}
		   	else{ std::cerr << "failed to open beta.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/g.txt");
		   	if( !ofs.fail()) {g.print(outputErrLatentDir + "/g.txt");}
		   	else{ std::cerr << "failed to open g.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/s.txt");
		   	if( !ofs.fail()) {S.print(outputErrLatentDir + "/s.txt");}
		   	else{ std::cerr << "failed to open s.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/u.txt");
		   	if( !ofs.fail()) {U.print(outputErrLatentDir + "/u.txt");}
		   	else{ std::cerr << "failed to open u.txt \n";}
		   	ofs.close();

		   	ofs.close();
		    assert( (beforeLq - nextLq)/std::abs(beforeLq) < epsilon  );
	    }
		if(first) first = false;
		update++;
	}
	pmswitch::InferenceData<Int, Real> ans(EZ, ES, EY, alpha, eta, beta, gamma, S, U, g, nextLq);
	// std::cerr << "init := " << update << ", g := " << std::endl;
	// dbData.g.print();
	return ans;
}

template<typename Int, typename Real>
pmswitch::InferenceData<Int, Real> pmswitch::Inference<Int, Real>::vb(bool updateDB, const Real epsilon,
																	  std::string outputErrLatentDir, std::string inputLatentDir) const{
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

	FixedSizeMultiVector<Real> firstEZ = initEZ(I, Js, maxJ, T, inputLatentDir + "/EZ.txt"); // I, maxJ, T
	FixedSizeMultiVector<Real> firstES = initES(I, Js, maxJ   , inputLatentDir + "/ES.txt");    // I, maxJ, 2
	FixedSizeMultiVector<Real> firstEY = initEY(I, Js, maxJ, N, inputLatentDir + "/EY.txt"); // I, maxJ, N

	FixedSizeMultiVector<Real> EZ = firstEZ;
	FixedSizeMultiVector<Real> ES = firstES;
	FixedSizeMultiVector<Real> EY = firstEY;

	FixedSizeMultiVector<Real> alpha = prior.alpha; // I, T-1, 2
	FixedSizeMultiVector<Real> beta  = prior.beta;  // I, 2
	FixedSizeMultiVector<Real> gamma = prior.gamma; // I, N
	FixedSizeMultiVector<Real> eta   = prior.eta;   // T, Ml.size(), maxMl

	FixedSizeMultiVector<Real> dig_v     = math::calDirExp<Real>(alpha, 1); // I, T-1, 2
	FixedSizeMultiVector<Real> dig_theta = math::calDirExp<Real>(beta,  0); // I, 2
	FixedSizeMultiVector<Real> dig_pi    = math::calDirExp<Real>(gamma, 0); // I, N
	FixedSizeMultiVector<Real> dig_f     = math::calDirExp<Real>(eta,   1); // T, L, maxMl
	FixedSizeMultiVector<Real> g         = dbData.g; // N, l, maxMl
	FixedSizeMultiVector<Real> dig_v_sig(0.0, I, T); // I, T

	FixedSizeMultiVector<bool> EZFilter(true, I, maxJ, T);   math::makeFilter(EZFilter,  I, Js, FixedSizeMultiVector<Int>(T, maxJ) );
	FixedSizeMultiVector<bool> ESFilter(true, I, maxJ, 2);   math::makeFilter(ESFilter,  I, Js, FixedSizeMultiVector<Int>(2, maxJ) );
	FixedSizeMultiVector<bool> EYFilter(true, I, maxJ, N);   math::makeFilter(EYFilter,  I, Js, FixedSizeMultiVector<Int>(N, maxJ) );
	FixedSizeMultiVector<bool> etaFilter(true, T, L, maxMl); math::makeFilter(etaFilter, T, FixedSizeMultiVector<Int>(L, T), ml);
	FixedSizeMultiVector<bool> gFilter(true, N, L, maxMl);   math::makeFilter(gFilter,   N, FixedSizeMultiVector<Int>(L, N), ml);

	FixedSizeMultiVector<Real> lnEZ = math::applied<Real>(EZ, std::log, EZFilter) ; // I, maxJ, T
	FixedSizeMultiVector<Real> lnES = math::applied<Real>(ES, std::log, ESFilter) ; // I, maxJ, 2
	FixedSizeMultiVector<Real> lnEY = math::applied<Real>(EY, std::log, EYFilter) ; // I, maxJ, N
	FixedSizeMultiVector<Real> lng  = math::applied<Real>(g,  std::log,  gFilter) ; // I, maxJ, N

	// {   // init g
	// 	Real eps = 1e-300;
	// 	math::filter<Int, Real>(g, gFilter, 0);
	// 	math::norm<Int, Real>(g, 1);
	// 	for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++)if(g(n,l,m) < eps) g(n,l,m) = eps;
	// }

    Real beforeLq = -1.0 * 1e30;
    Real nextLq = -1.0 * 1e30;
    bool first = true;
    Int update = 0;

    Real INF = 1e300;

    while( first || std::abs(nextLq - beforeLq)/std::abs(beforeLq) >=  epsilon){
	    beforeLq = nextLq;
	    nextLq = 0.0;

	    // M-step
		udpateF(eta, etaFilter, dig_f, I, Js, T, L, ml, prior, X, EZ, ES);
		updatePi(gamma, dig_pi, I, Js, N, prior, EY);
		updateV(alpha, dig_v, I, Js, T, prior,EZ);
		updateTheta(beta, dig_theta, I, Js, prior, ES);

		if(updateDB){ // update < 275 -> pass ,update < 276 -> not pass
			updateG(g, lng, gFilter, I, Js, maxJ, T, L, ml, maxMl, N, EY, ES, X);
			// std::cerr << "update := " << update << ", g := " << std::endl;
			// g.print();
		}

	    // E-step
    	updateEZ(EZ, lnEZ, EZFilter, I, Js, maxJ, T, L, ml, N, dig_f, ES, X, dig_v);
     	updateES(ES, lnES, ESFilter, I, Js, maxJ, T, L, ml, N, dig_f, EZ, lng, EY, X, dig_theta);
		updateEY(EY, lnEY, EYFilter, I, Js, maxJ, L, ml, N, lng, ES, X, dig_pi);

		{// nextLq = 0.0
			nextLq = 0.0;
			nextLq += math::calELogDir<Real>(dig_f,     prior.eta,    1, etaFilter); // T, L, maxMl
			nextLq += math::calELogDir<Real>(dig_v,     prior.alpha,  1); // I, T-1, 2
			nextLq += math::calELogDir<Real>(dig_theta, prior.beta,   0); // I, 2
			nextLq += math::calELogDir<Real>(dig_pi,    prior.gamma,  0); // I, N

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
					// nextLq     += tmpXLogY[ (EY(i,j,n) * ES(i,j,1) * X(i,j,l,m) != 0.0)  ];
			}

			nextLq -= math::calELogDir<Real>(dig_f,     eta,    1, etaFilter);  // T, L, maxMl
			nextLq -= math::calELogDir<Real>(dig_v,     alpha,  1); 			 // I, T-1, 2
			nextLq -= math::calELogDir<Real>(dig_theta, beta,   0); 			 // I, 2
			nextLq -= math::calELogDir<Real>(dig_pi,    gamma,  0); 			 // I, N

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

	    	std::ofstream ofsEZ(outputErrLatentDir + "/EZ.txt");
		   	if( !ofsEZ.fail()) {firstEZ.print(outputErrLatentDir + "/EZ.txt");}
		   	ofsEZ.close();

	    	std::ofstream ofsES(outputErrLatentDir + "/ES.txt");
		   	if( !ofsES.fail()) {firstES.print(outputErrLatentDir + "/ES.txt");}
		   	ofsES.close();

	    	std::ofstream ofsEY(outputErrLatentDir + "/EY.txt");
		   	if( !ofsEY.fail()) {firstEY.print(outputErrLatentDir + "/EY.txt");}
		   	ofsEY.close();

		   	std::ofstream ofs(outputErrLatentDir + "/EZnew.txt");
		   	if( !ofs.fail()) {EZ.print(outputErrLatentDir + "/EZnew.txt");}
		   	else{ std::cerr << "failed to open EZnew.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/ESnew.txt");
		   	if( !ofs.fail()) {ES.print(outputErrLatentDir + "/ESnew.txt");}
		   	else{ std::cerr << "failed to open ESnew.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/EYnew.txt");
		   	if( !ofs.fail()) {EY.print(outputErrLatentDir + "/EYnew.txt");}
		   	else{ std::cerr << "failed to open EYnew.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/alpha.txt");
		   	if( !ofs.fail()) {alpha.print(outputErrLatentDir + "/alpha.txt");}
		   	else{ std::cerr << "failed to open alpha.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/beta.txt");
		   	if( !ofs.fail()) {beta.print(outputErrLatentDir + "/beta.txt");}
		   	else{ std::cerr << "failed to open beta.txt \n";}
		   	ofs.close();

		   	ofs.open(outputErrLatentDir + "/g.txt");
		   	if( !ofs.fail()) {g.print(outputErrLatentDir + "/g.txt");}
		   	else{ std::cerr << "failed to open g.txt \n";}
		   	ofs.close();

		    assert( (beforeLq - nextLq)/std::abs(beforeLq) < epsilon  );
	    }
		if(first) first = false;
		update++;
	}
	pmswitch::InferenceData<Int, Real> ans(EZ, ES, EY, alpha, eta, beta, gamma, g, nextLq);
	// std::cerr << "init := " << update << ", g := " << std::endl;
	// dbData.g.print();
	return ans;
}

// EZ_{i,j, }     :: T-dimensional non-negative simplex // EZ_{i,j,k   }
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::initEZ(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int T, std::string inputPath) const{
	std::ifstream ifs(inputPath);
	if(ifs.fail()) {
		ifs.close();
		pmswitch::FixedSizeMultiVector<Real> EZ(0.0, I, maxJ, T);
		std::vector<Real> alpha; for(Int k = 0; k < T; k++) alpha.push_back((Real)1.0);
		std::random_device rd;
		std::mt19937 gen(rd());
		for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
			std::vector<Real> rand = pmswitch::math::random::sampleDirichlet<std::mt19937, Real>(gen, alpha);
			for(Int k = 0; k < T; k++) EZ(i,j,k) = rand[k];
		}
		return EZ;
	}else{
		std::cerr << "parse EZ from file " << std::endl;
		ifs.close();
	    FixedSizeMultiVector<Real> EZ = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(inputPath, FixedSizeMultiVectorCreator<Real>::strToReal);
	    if(EZ.dim() != 3){ pmswitch::die("invalid dimension @ pmswitch::Inference<Int, Real>::initEZ"); }
	    if(EZ.size(0) != I || EZ.size(1) != maxJ || EZ.size(2) != T){pmswitch::die("invalid dimension size @ pmswitch::Inference<Int, Real>::initEZ");}
	    return EZ;
	}
}

// ES_{i, }	   :: 2-dimensional non-negative simplex // ES_{i,{0,1} }
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::initES(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, std::string inputPath) const{
	std::ifstream ifs(inputPath);
	if(ifs.fail()) {
		ifs.close();
		pmswitch::FixedSizeMultiVector<Real> ES(0.0, I, maxJ, 2);
		std::vector<Real> alpha = {1.0, 1.0};
		std::random_device rd;
		std::mt19937 gen(rd());
		for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
			std::vector<Real> rand = pmswitch::math::random::sampleDirichlet<std::mt19937, Real>(gen, alpha);
			ES(i,j,0) = rand[0]; ES(i,j,1) = rand[1];
		}
		return ES;
	}else{
		std::cerr << "parse ES from file " << std::endl;
		ifs.close();
	    FixedSizeMultiVector<Real> ES = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(inputPath, FixedSizeMultiVectorCreator<Real>::strToReal);
	    if(ES.dim() != 3){ pmswitch::die("invalid dimension @ pmswitch::Inference<Int, Real>::initES"); }
	    if(ES.size(0) != I || ES.size(1) != maxJ || ES.size(2) != 2){pmswitch::die("invalid dimension size @ pmswitch::Inference<Int, Real>::initES");}
	    return ES;
	}
}

// EY_{i, }       :: N-dimensional non-negative simplex // EY_{i,n}
template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::initEY(Int I, FixedSizeMultiVector<Int> Js, Int maxJ, Int N, std::string inputPath) const{
	std::ifstream ifs(inputPath);
	if(ifs.fail()) {
		ifs.close();
		pmswitch::FixedSizeMultiVector<Real> EY(0.0, I, maxJ, N);
		std::vector<Real> alpha; for(Int n = 0; n < N; n++) alpha.push_back((Real)1.0);
		std::random_device rd;
		std::mt19937 gen(rd());
		for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
			std::vector<Real> rand = pmswitch::math::random::sampleDirichlet<std::mt19937, Real>(gen, alpha);
			for(Int n = 0; n < N; n++){ EY(i,j,n) = rand[n];}
		}
		return EY;
	}else{
		std::cerr << "parse EY from file " << std::endl;
		ifs.close();
	    FixedSizeMultiVector<Real> EY = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(inputPath, FixedSizeMultiVectorCreator<Real>::strToReal);
	    if(EY.dim() != 3){ pmswitch::die("invalid dimension @ pmswitch::Inference<Int, Real>::initEY"); }
	    if(EY.size(0) != I || EY.size(1) != maxJ || EY.size(2) != N){pmswitch::die("invalid dimension size @ pmswitch::Inference<Int, Real>::initEY");}
	    return EY;
	}
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateEZ(FixedSizeMultiVector<Real> &EZ, FixedSizeMultiVector<Real> &lnEZ, const FixedSizeMultiVector<bool> &EZFilter,
			  const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml, const Int &N,
	          const FixedSizeMultiVector<Real> &dig_f, const FixedSizeMultiVector<Real> &ES,
	          const FixedSizeMultiVector<Int>  &X,     const FixedSizeMultiVector<Real> &dig_v) const {
    for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++)for(Int k = 0; k < T; k++) EZ(i,j,k) = 0.0;

    for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T; k++){
		for(Int kk = 0; kk <= k-1; kk++){ EZ(i,j,k) += dig_v(i,kk,1);}
		if(k < T-1) {EZ(i,j,k) += dig_v(i,k,0);}
		for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
			EZ(i,j,k) += ES(i,j,0) * X(i,j,l,m) * dig_f(k, l, m);
		}
    }
	math::subMax<Real>(EZ, 1);
	math::apply<Real>(EZ, std::exp);
	math::filter<Real>(EZ, EZFilter, 0.0);
	math::norm<Real>(EZ, 1);
	lnEZ = math::applied<Real>(EZ, std::log, EZFilter);
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateES(FixedSizeMultiVector<Real> &ES, FixedSizeMultiVector<Real> &lnES, const FixedSizeMultiVector<bool> &ESFilter,
			 const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml, const Int &N,
			  const FixedSizeMultiVector<Real> &dig_f, const FixedSizeMultiVector<Real> &EZ,
			  const FixedSizeMultiVector<Real> &lng,   const FixedSizeMultiVector<Real> &EY,
			  const FixedSizeMultiVector<Int>  &X,     const FixedSizeMultiVector<Real> &dig_theta) const {
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
	math::subMax<Real>(ES, 1);
	math::apply<Real>(ES, std::exp);
	math::filter<Real>(ES, ESFilter, 0.0);
	math::norm<Real>(ES, 1);
	lnES = math::applied<Real>(ES, std::log, ESFilter);
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateEY(FixedSizeMultiVector<Real> &EY, FixedSizeMultiVector<Real> &lnEY, const FixedSizeMultiVector<bool> &EYFilter,
			  const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &L, const FixedSizeMultiVector<Int> &ml, const Int &N,
			  const FixedSizeMultiVector<Real> &lng,   const FixedSizeMultiVector<Real> &ES,
			  const FixedSizeMultiVector<Int>  &X,     const FixedSizeMultiVector<Real> &dig_pi) const {
	std::vector<Real> tmpXLogY = {0, 0};
	for(Int i = 0; i < I; i++)for(Int j = 0; j < maxJ; j++)for(Int n = 0; n < N; n++){EY(i,j,n) = 0.0;}

	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int n = 0; n < N; n++){
		EY(i,j,n) += dig_pi(i,n);
		for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
			tmpXLogY[1] = X(i,j,l,m) * ES(i,j,1) * lng(n,l,m);
			EY(i,j,n) += tmpXLogY[!std::isnan(tmpXLogY[1])];
		}
	}
	math::subMax<Real>(EY, 1);
	math::apply<Real>(EY, std::exp);
	math::filter<Real>(EY, EYFilter, 0.0);
	math::norm<Real>(EY, 1);
	lnEY = math::applied<Real>(EY, std::log, EYFilter);
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::udpateF(FixedSizeMultiVector<Real> &eta, const FixedSizeMultiVector<bool> &etaFilter, FixedSizeMultiVector<Real> &dig_f,
			const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml,
			const PriorParameters<Int, Real> &prior,   const FixedSizeMultiVector<Int>  &X,
			const FixedSizeMultiVector<Real> &EZ, const FixedSizeMultiVector<Real> &ES) const {
		eta = prior.eta;
		for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T; k++){
			for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++){
				eta(k,l,m) += X(i,j,l,m) * EZ(i,j,k) * ES(i,j,0);
			}
		}
		dig_f = math::calDirExp<Real>(eta, 1, etaFilter);
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateTheta(FixedSizeMultiVector<Real> &beta, FixedSizeMultiVector<Real> &dig_theta,
				 const Int &I, const FixedSizeMultiVector<Int> &Js,
				 const PriorParameters<Int, Real> &prior, const FixedSizeMultiVector<Real> &ES) const {
	beta  = prior.beta;  // I, 2
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
		beta(i, 0) += ES(i,j,0); beta(i, 1) += ES(i,j,1);
	}
	dig_theta = math::calDirExp<Real>(beta,  0);
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateThetaWithBetaDistributed(FixedSizeMultiVector<Real> &beta, FixedSizeMultiVector<Real> &dig_theta,
																	const Int &I, const FixedSizeMultiVector<Int> &Js,
																	const FixedSizeMultiVector<Real> &avgBeta,
																	const FixedSizeMultiVector<Real> &ES) const {
	beta  = avgBeta;  // I, 2
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++){
		beta(i, 0) += ES(i,j,0); beta(i, 1) += ES(i,j,1);
	}
	dig_theta = math::calDirExp<Real>(beta,  0);
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updatePi(FixedSizeMultiVector<Real> &gamma, FixedSizeMultiVector<Real> &dig_pi,
			  const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &N,
			  const PriorParameters<Int, Real> &prior, const FixedSizeMultiVector<Real> &EY) const {
	gamma = prior.gamma; // I, N
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int n = 0; n < N; n++){
			gamma(i,n) += EY(i,j,n);
	}
 	dig_pi    = math::calDirExp<Real>(gamma, 0); // I, N
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateV(FixedSizeMultiVector<Real> &alpha, FixedSizeMultiVector<Real> &dig_v,
	                       const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &T,
	                       const PriorParameters<Int, Real> &prior, const FixedSizeMultiVector<Real> &EZ) const {
	alpha = prior.alpha; // I, T-1, 2
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T-1; k++){
		alpha(i,k,0) += EZ(i,j,k);
		for(Int t = k+1; t < T; t++){ alpha(i,k,1) += EZ(i,j,t);}
	}
	dig_v = math::calDirExp<Real>(alpha, 1); // I, T-1, 2
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateVWithAlphaDistributed(FixedSizeMultiVector<Real> &alpha, FixedSizeMultiVector<Real> &dig_v,
	                       const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &T,
	                       const FixedSizeMultiVector<Real> &avgAlpha, const FixedSizeMultiVector<Real> &EZ) const {
	alpha = avgAlpha; // I, T-1, 2
	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int k = 0; k < T-1; k++){
		alpha(i,k,0) += EZ(i,j,k);
		for(Int t = k+1; t < T; t++){ alpha(i,k,1) += EZ(i,j,t);}
	}
	dig_v = math::calDirExp<Real>(alpha, 1); // I, T-1, 2
}


template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::updateG(FixedSizeMultiVector<Real> &g, FixedSizeMultiVector<Real> &lng, const FixedSizeMultiVector<bool> &gFilter,
			 const Int &I, const FixedSizeMultiVector<Int> &Js, const Int &maxJ, const Int &T, const Int &L, const FixedSizeMultiVector<Int> &ml, Int maxMl, const Int &N,
			 const FixedSizeMultiVector<Real> &EY, const FixedSizeMultiVector<Real> &ES,
			  const FixedSizeMultiVector<Int> &X, Real eps) const { // update < 275 -> pass ,update < 276 -> not pass
	for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < maxMl; m++){ g(n,l,m) = 0.0;}

	for(Int i = 0; i < I; i++)for(Int j = 0; j < Js(i); j++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++)for(Int n = 0; n < N; n++){
		g(n,l,m) += X(i,j,l,m) * EY(i,j,n) * ES(i,j,1);
	}
	math::filter<Real>(g, gFilter, 0);
	math::norm<Real>(g, 1);
	for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++)for(Int m = 0; m < ml[l]; m++)if(g(n,l,m) < eps) g(n,l,m) = eps;
	lng  = math::applied<Real>(g,  std::log,  gFilter) ; // I, maxJ, N
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::makeAvgAlpha(const FixedSizeMultiVector<Real> &s, Int I, Int T) const{
	using namespace pmswitch;
	FixedSizeMultiVector<Real> ans((Real)0, I, T-1, 2);
	for(Int i = 0; i < I; i++)for(Int k = 0; k < T-1; k++){
		ans(i,k,0) = (Real) 1.0;
		ans(i,k,1) = (Real) ( s(0) / s(1) );
	}
	return ans;
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::makeAvgAlpha(FixedSizeMultiVector<Real> &ans, const FixedSizeMultiVector<Real> &s, Int I, Int T) const{
	for(Int i = 0; i < I; i++)for(Int k = 0; k < T-1; k++){
		ans(i,k,0) = (Real) 1.0;
		ans(i,k,1) = (Real) ( s(0) / s(1) );
	}
}


template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::Inference<Int, Real>::makeAvgBeta(const FixedSizeMultiVector<Real> &u, Int I) const{
	using namespace pmswitch;
	FixedSizeMultiVector<Real> ans((Real)0, I, 2);
	for(Int i = 0; i < I; i++){
		ans(i,0) = (Real) 1.0;
		ans(i,1) = (Real) ( u(0) / u(1) );
	}
	return ans;
}

template<typename Int, typename Real>
void pmswitch::Inference<Int, Real>::makeAvgBeta(FixedSizeMultiVector<Real> &ans, const FixedSizeMultiVector<Real> &u, Int I) const{
	for(Int i = 0; i < I; i++){
		ans(i,0) = (Real) 1.0;
		ans(i,1) = (Real) ( u(0) / u(1) );
	}
}

		// FixedSizeMultiVector<Real, Int> makeAvgBeta(const FixedSizeMultiVector<Real, Int> &u, Int I, Int T) const;

		// void makeAvgAlpha(FixedSizeMultiVector<Real, Int> &ans, const FixedSizeMultiVector<Real, Int> &s, Int I, Int T);
		// void makeAvgBeta(FixedSizeMultiVector<Real, Int> &ans,  const FixedSizeMultiVector<Real, Int> &u, Int I, Int T);

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
													   Real _s0, Real _s1,
													   Real _u0, Real _u1,
													   std::string alphaPath, std::string betaPath,
													   std::string gammaPath, std::string etaPath ,
													   std::string sPath,     std::string uPath,
													   Real (*strToE)(std::string)){
	using namespace pmswitch;
	InputFileParser<Int, Real> parser;

	FeatureData<Int, Real> fvData = parser.parseFeatureFile(pathFV);
	DBData<Int, Real> dbData = parser.parseDBFile(pathDB);
	// PriorParametersCreator<Int, Real> priorCreator;
	PriorParameters<Int, Real> prior = PriorParametersCreator<Int, Real>::createPriorParameters(fvData,  dbData, _T, _alpha, _beta0, _beta1, _gamma, _eta, _s0, _s1, _u0, _u1,
 																								alphaPath,  betaPath, gammaPath,  etaPath, sPath, uPath, strToE);
	Inference<Int, Real> inference(fvData, dbData, prior, _T);
	return inference;
}


template<typename Int, typename Real>
pmswitch::Inference<Int, Real> pmswitch::InferenceCreator<Int, Real>::createInference(pmswitch::Parameters<Int, Real> param){
	return pmswitch::InferenceCreator<Int, Real>::createInference(param.pathFV, param.pathDB, param.T, param.alpha, param.beta0, param.beta1,
																  param.gamma,  param.eta, param.s0, param.s1, param.u0, param.u1,
																  param.alphaPath, param.betaPath, param.gammaPath, param.etaPath, param.sPath, param.uPath,
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
