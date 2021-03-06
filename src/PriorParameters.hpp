#ifndef PriorParameters_h
#define PriorParameters_h

#include "ExceptionUtil.hpp"
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
						pmswitch::FixedSizeMultiVector<Real> eta,
						pmswitch::FixedSizeMultiVector<Real> s,
						pmswitch::FixedSizeMultiVector<Real> u);
		void print() const;
		void print(std::string path) const;
		const pmswitch::FixedSizeMultiVector<Real> alpha;
		const pmswitch::FixedSizeMultiVector<Real> beta;
		const pmswitch::FixedSizeMultiVector<Real> gamma;
		const pmswitch::FixedSizeMultiVector<Real> eta;
		const pmswitch::FixedSizeMultiVector<Real> s;
		const pmswitch::FixedSizeMultiVector<Real> u;
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
		static pmswitch::PriorParameters<Int, Real> createPriorParameters(
														   pmswitch::FeatureData<Int,Real> fvData,
														   pmswitch::DBData<Int,Real> dbData,
														   Int T,
														   Real _alpha,
														   Real _beta0, Real _beta1,
														   Real _gamma, Real _eta,
														   Real _s0,    Real _s1,
														   Real _u0,    Real _u1,
														   std::string alphaPath = "", std::string betaPath = "",
														   std::string gammaPath = "", std::string etaPath = "",
														   std::string sPath = "",     std::string uPath = "",
														   Real (*strToE)(std::string) = NULL);
	private:
		static pmswitch::FixedSizeMultiVector<Real> makeAlpha(Int I, Int T, Real _alpha);
		template<typename String> static pmswitch::FixedSizeMultiVector<Real> makeAlpha(Int I, Int T, std::string path, Real (*strToE)(String));

		static pmswitch::FixedSizeMultiVector<Real> makeBeta(Int I, Real _beta0, Real _beta1);
		template<typename String> static pmswitch::FixedSizeMultiVector<Real> makeBeta(Int I, std::string path, Real (*strToE)(String));

		static pmswitch::FixedSizeMultiVector<Real> makeGamma(Int I, Int N, Real _gamma);
		template<typename String> static pmswitch::FixedSizeMultiVector<Real> makeGamma(Int I, Int N, std::string path, Real (*strToE)(String) );

		static pmswitch::FixedSizeMultiVector<Real> makeEta(Int T, Int L, pmswitch::FixedSizeMultiVector<Int> Ml, Real _eta);
		template<typename String> static pmswitch::FixedSizeMultiVector<Real> makeEta(Int T, Int L, pmswitch::FixedSizeMultiVector<Int> Ml, std::string path, Real (*strToE)(String) );

		static pmswitch::FixedSizeMultiVector<Real> makeS(Real _s0, Real _s1);
		template<typename String> static pmswitch::FixedSizeMultiVector<Real> makeS(std::string path, Real (*strToE)(String));

		static pmswitch::FixedSizeMultiVector<Real> makeU(Real _s0, Real _s1);
		template<typename String> static pmswitch::FixedSizeMultiVector<Real> makeU(std::string path, Real (*strToE)(String));

	};
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeAlpha(Int I, Int T, Real _alpha){
	pmswitch::FixedSizeMultiVector<Real> alpha(_alpha, I, T-1, 2);
	for(Int i = 0; i < I; i++)for(Int t = 0; t < T-1; t++){
		alpha(i, t, 0) = 1.0; alpha(i, t, 1) = _alpha;
	}
	return alpha;
}

template<typename Int, typename Real>
template<typename String>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeAlpha(Int I, Int T, std::string path, Real (*strToE)(String)){
	using namespace pmswitch;
	pmswitch::FixedSizeMultiVector<Real> alpha = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(path,  strToE);
	if( alpha.dim() != 3  ){ die("alpha.dim() differs @ PriorParametersCreator<Int, Real>::makeAlpha(Int I, Int T, std::string path, Real (*strToE)(String))"); }
	if( alpha.size(0) != I || alpha.size(1) != T || alpha.size(2) != 2){
		die("alpha.size(each) differs @ PriorParametersCreator<Int, Real>::makeAlpha(Int I, Int T, std::string path, Real (*strToE)(String))");
	}
	return alpha;
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeBeta(Int I, Real _beta0, Real _beta1){
	using namespace pmswitch;
	pmswitch::FixedSizeMultiVector<Real> beta(0, I, 2);
	for(Int i = 0; i < I; i++){ beta(i,0) = _beta0; beta(i,1) = _beta1; }
	return beta;
}

template<typename Int, typename Real>
template<typename String>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeBeta(Int I, std::string path, Real (*strToE)(String) ){
	using namespace pmswitch;
	pmswitch::FixedSizeMultiVector<Real> beta = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(path,  strToE);
	if( beta.dim() != 2  ){ die("beta.dim() differs @ PriorParametersCreator<Int, Real>::makebeta(Int I, Int T, std::string path, Real (*strToE)(String))"); }
	if( beta.size(0) != I || beta.size(1) != 2){
		die("beta.size(each) differs @ PriorParametersCreator<Int, Real>::makebeta(Int I, Int T, std::string path, Real (*strToE)(String))");
	}
	return beta;
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeGamma(Int I, Int N, Real _gamma){
	pmswitch::FixedSizeMultiVector<Real> gamma(_gamma, I, N);
	return gamma;
}

template<typename Int, typename Real>
template<typename String>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeGamma(Int I, Int N, std::string path, Real (*strToE)(String)){
	pmswitch::FixedSizeMultiVector<Real> gamma = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(path,  strToE);
	if( gamma.dim()   != 2  ){ die("gamma.dim() differs @ PriorParametersCreator<Int, Real>::makeGamma(Int I, Int T, std::string path, Real (*strToE)(String))"); }
	if( gamma.size(0) != I || gamma.size(1) != N){
		die("gamma.size(each) differs @ PriorParametersCreator<Int, Real>::makeGamma(Int I, Int T, std::string path, Real (*strToE)(String))");
	}
	return gamma;
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeEta(Int T, Int L, pmswitch::FixedSizeMultiVector<Int> Ml, Real _eta){
	Int maxMl = 0;
	for(Int l = 0; l < Ml.size(); l++) {
		maxMl = std::max(maxMl, Ml[l]);
		if(Ml[l] < 0) {
			pmswitch::die("Ml is invalid value. @ PriorParametersCreator<Int, Real>::makeEta(Int T, Int L, pmswitch::FixedSizeMultiVector<Int> Ml, std::string path, Real (*strToE)(String))");
		}
	}
	pmswitch::FixedSizeMultiVector<Real> eta(_eta, T, Ml.size(), maxMl);
	for(Int t = 0; t < T; t++)for(Int l = 0; l < Ml.size(); l++)for(Int m = 0; m < Ml(l); m++){
		eta(t,l,m) = _eta;
	}
	return eta;
}

template<typename Int, typename Real>
template<typename String>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeEta(Int T, Int L, pmswitch::FixedSizeMultiVector<Int> Ml, std::string path, Real (*strToE)(String)){
	Int maxMl = 0;
	for(Int l = 0; l < Ml.size(); l++) {
		maxMl = std::max(maxMl, Ml[l]);
		if(Ml[l] < 0) {
			pmswitch::die("Ml is invalid value. @ PriorParametersCreator<Int, Real>::makeEta(Int T, Int L, pmswitch::FixedSizeMultiVector<Int> Ml, std::string path, Real (*strToE)(String))");
		}
	}
	pmswitch::FixedSizeMultiVector<Real> eta = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(path,  strToE);
	if( eta.dim()   != 3  ){ die("eta.dim() differs @ PriorParametersCreator<Int, Real>::makeeta(Int I, Int T, std::string path, Real (*strToE)(String))"); }
	if( eta.size(0) != T || eta.size(1) != Ml.size() || eta.size(2) != maxMl){
		die("eta.size(each) differs @ PriorParametersCreator<Int, Real>::makebeta(Int I, Int T, std::string path, Real (*strToE)(String))");
	}
	return eta;
}

template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeS(Real _s0, Real _s1){
	pmswitch::FixedSizeMultiVector<Real> s(0, 2);
	s(0) = _s0; s(1) = _s1;
	return s;
}

template<typename Int, typename Real>
template<typename String>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeS(std::string path, Real (*strToE)(String)){
	pmswitch::FixedSizeMultiVector<Real> s = FixedSizeMultiVectorCreator<Real>::createFixedSizeMultiVector(path,  strToE);
	if( s.dim()   != 1  ){ die("s.dim() differs @ pmswitch::PriorParametersCreator<Int, Real>::makeS(std::string path, Real (*strToE)(String))"); }
	if( s.size(0) != 2  ){
		die("s.size(each) differs @ pmswitch::PriorParametersCreator<Int, Real>::makeS(std::string path, Real (*strToE)(String))");
	}
	return s;
}


template<typename Int, typename Real>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeU(Real _u0, Real _u1){
	return makeS(_u0,_u1);
}

template<typename Int, typename Real>
template<typename String>
pmswitch::FixedSizeMultiVector<Real> pmswitch::PriorParametersCreator<Int, Real>::makeU(std::string path, Real (*strToE)(String)){
	return makeS(path, strToE);
}


// static pmswitch::FixedSizeMultiVector<Real, Int> makeS(Real _s0, Real _s1);
// template<typename String> static pmswitch::FixedSizeMultiVector<Real, Int> makeS(std::string path, Real (*strToE)(String));


template<typename Int, typename Real>
pmswitch::PriorParameters<Int, Real>::PriorParameters(pmswitch::FixedSizeMultiVector<Real> _alpha,
													  pmswitch::FixedSizeMultiVector<Real> _beta,
													  pmswitch::FixedSizeMultiVector<Real> _gamma,
													  pmswitch::FixedSizeMultiVector<Real> _eta,
													  pmswitch::FixedSizeMultiVector<Real> _s,
													  pmswitch::FixedSizeMultiVector<Real> _u) : alpha(_alpha), beta(_beta), gamma(_gamma), eta(_eta), s(_s), u(_u){

}

template<typename Int, typename Real>
void pmswitch::PriorParameters<Int, Real>::print() const{
	std::cerr<<"==================== PriorParameters.alpha ====================" << std::endl;
	alpha.print(std::cerr);
	std::cerr<<"==================== PriorParameters.beta ====================" << std::endl;
	beta.print(std::cerr);
	std::cerr<<"==================== PriorParameters.gamma ====================" << std::endl;
	gamma.print(std::cerr);
	std::cerr<<"==================== PriorParameters.eta ====================" << std::endl;
	eta.print(std::cerr);
	std::cerr<<"==================== PriorParameters.s ====================" << std::endl;
	s.print(std::cerr);
	std::cerr<<"==================== PriorParameters.u ====================" << std::endl;
	u.print(std::cerr);
}

template<typename Int, typename Real>
void pmswitch::PriorParameters<Int, Real>::print(std::string path) const{
	alpha.print(path + "/alpha.txt");
	beta.print(path + "/beta.txt");
	eta.print(path + "/eta.txt");
	gamma.print(path + "/gamma.txt");
	s.print(path + "/s.txt");
	u.print(path + "/u.txt");
}

template<typename Int, typename Real>
pmswitch::PriorParametersCreator<Int, Real>::PriorParametersCreator(){

}

template<typename Int, typename Real>
pmswitch::PriorParameters<Int, Real> pmswitch::PriorParametersCreator<Int, Real>::createPriorParameters(
   														   pmswitch::FeatureData<Int,Real> fvData,
														   pmswitch::DBData<Int,Real> dbData,
														   Int T,
														   Real _alpha,
														   Real _beta0, Real _beta1,
														   Real _gamma, Real _eta,
														   Real _s0,    Real _s1,
														   Real _u0,    Real _u1,
														   std::string alphaPath, std::string betaPath,
														   std::string gammaPath, std::string etaPath,
														   std::string sPath,     std::string uPath,
														   Real (*strToE)(std::string) ){
	Int I = fvData.I;
	Int N = dbData.N;
	const pmswitch::FixedSizeMultiVector<Int> Ml = fvData.Ml;

	Int maxMl = 0; for(Int l = 0; l < Ml.size(); l++) maxMl = std::max(maxMl, Ml[l]);
	pmswitch::FixedSizeMultiVector<Real> alpha; //(_alpha, I, T-1, 2);
	pmswitch::FixedSizeMultiVector<Real> beta;  //(0.0, I, 2);
	pmswitch::FixedSizeMultiVector<Real> gamma; //(_gamma, I, N);
	pmswitch::FixedSizeMultiVector<Real> eta;   //(_eta, T, Ml.size(), maxMl);
	pmswitch::FixedSizeMultiVector<Real> s;     //(2)
	pmswitch::FixedSizeMultiVector<Real> u;     //(2)


	std::ifstream ifs;
	ifs.open(alphaPath);
   	if(!ifs.fail()) { ifs.close(); alpha = makeAlpha(I, T, alphaPath, strToE);}
   	else { 			  ifs.close(); alpha = makeAlpha(I, T, _alpha);}

   	ifs.open(betaPath);
   	if(!ifs.fail()) { ifs.close(); beta  = makeBeta(I, betaPath, strToE);}
   	else { 			  ifs.close(); beta  = makeBeta(I, _beta0, _beta1);}

   	ifs.open(gammaPath);
   	if(!ifs.fail()) { ifs.close(); gamma = makeGamma(I, N, gammaPath, strToE);}
   	else { 			  ifs.close(); gamma = makeGamma(I, N, _gamma);}

   	ifs.open(etaPath);
   	if(!ifs.fail()) { ifs.close(); eta = makeEta(T, Ml.size(), Ml, etaPath, strToE);}
   	else { 			  ifs.close(); eta = makeEta(T, Ml.size(), Ml, _eta);}

   	ifs.open(sPath);
   	if(!ifs.fail()) { ifs.close(); s = makeS(sPath, strToE);}
   	else { 			  ifs.close(); s = makeS(_s0, _s1);}

   	ifs.open(uPath);
   	if(!ifs.fail()) { ifs.close(); u = makeU(uPath, strToE);}
   	else { 			  ifs.close(); u = makeU(_u0, _u1);}


	pmswitch::PriorParameters<Int, Real> prior(alpha, beta, gamma, eta, s, u);
	return prior;
}


#endif
