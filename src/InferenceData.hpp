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
		InferenceData(FixedSizeMultiVector<Real> _EZ, FixedSizeMultiVector<Real> _ES, FixedSizeMultiVector<Real> _EY,
					  FixedSizeMultiVector<Real> _alpha, FixedSizeMultiVector<Real> _eta,
					  FixedSizeMultiVector<Real> _beta, FixedSizeMultiVector<Real> _gamma,
					  FixedSizeMultiVector<Real> _s, FixedSizeMultiVector<Real> _u,
					  FixedSizeMultiVector<Real> _g, Real _Lq = 0.0);

		InferenceData(FixedSizeMultiVector<Real> _EZ, FixedSizeMultiVector<Real> _ES, FixedSizeMultiVector<Real> _EY,
					  FixedSizeMultiVector<Real> _alpha, FixedSizeMultiVector<Real> _eta,
					  FixedSizeMultiVector<Real> _beta, FixedSizeMultiVector<Real> _gamma,
					  FixedSizeMultiVector<Real> _g, Real _Lq = 0.0);

		void printHyperParameters(std::string baseDir) 		const;
		void printLatents(std::string baseDir)		   		const;
		void printDB(std::string baseDir)              		const;
		void printAllData(std::string baseDir)         		const;
		void printGammaHyperParameters(std::string baseDir) const;

		FixedSizeMultiVector<Real> getAlpha() const {return alpha;}
		FixedSizeMultiVector<Real> getBeta()  const {return beta;}
		FixedSizeMultiVector<Real> getGamma() const {return gamma;}
		FixedSizeMultiVector<Real> getEta()   const {return eta;}
		FixedSizeMultiVector<Real> getG()     const {return g;}
		FixedSizeMultiVector<Real> getS()     const {return g;}
		FixedSizeMultiVector<Real> getU()     const {return g;}

		FixedSizeMultiVector<Real> getEZ()    const {return EZ;}
		FixedSizeMultiVector<Real> getES()    const {return ES;}
		FixedSizeMultiVector<Real> getEY()    const {return EY;}

		Real getLq() const {return Lq;}

	private:
		FixedSizeMultiVector<Real> EZ;
		FixedSizeMultiVector<Real> ES;
		FixedSizeMultiVector<Real> EY;

		FixedSizeMultiVector<Real> alpha;
		FixedSizeMultiVector<Real> eta;
		FixedSizeMultiVector<Real> beta;
		FixedSizeMultiVector<Real> gamma;

		FixedSizeMultiVector<Real> s;
		FixedSizeMultiVector<Real> u;

		FixedSizeMultiVector<Real> g;
		Real Lq;
	};
}

template<typename Int , typename Real >
pmswitch::InferenceData<Int, Real>::InferenceData(FixedSizeMultiVector<Real> _EZ, FixedSizeMultiVector<Real> _ES, FixedSizeMultiVector<Real> _EY,
					  FixedSizeMultiVector<Real> _alpha, FixedSizeMultiVector<Real> _eta,
					  FixedSizeMultiVector<Real> _beta, FixedSizeMultiVector<Real> _gamma,
					  FixedSizeMultiVector<Real> _g, Real _Lq) : EZ(_EZ), ES(_ES), EY(_EY), alpha(_alpha), eta(_eta), beta(_beta), gamma(_gamma), g(_g), Lq(_Lq) {
}

template<typename Int , typename Real >
pmswitch::InferenceData<Int, Real>::InferenceData(FixedSizeMultiVector<Real> _EZ, FixedSizeMultiVector<Real> _ES, FixedSizeMultiVector<Real> _EY,
					  FixedSizeMultiVector<Real> _alpha, FixedSizeMultiVector<Real> _eta,
					  FixedSizeMultiVector<Real> _beta, FixedSizeMultiVector<Real> _gamma,
					  FixedSizeMultiVector<Real> _s, FixedSizeMultiVector<Real> _u,
					  FixedSizeMultiVector<Real> _g,
					  Real _Lq) : EZ(_EZ), ES(_ES), EY(_EY), alpha(_alpha), eta(_eta), beta(_beta), gamma(_gamma), g(_g), s(_s), u(_u), Lq(_Lq) {
}

template<typename Int , typename Real >
void pmswitch::InferenceData<Int, Real>::printHyperParameters(std::string baseDir) const {
	std::string alphaPath = baseDir + "/alpha.txt";
	std::string betaPath  = baseDir + "/beta.txt";
	std::string gammaPath = baseDir + "/gamma.txt";
	std::string etaPath   = baseDir + "/eta.txt";
	std::string sPath     = baseDir + "/s.txt";
	std::string uPath     = baseDir + "/u.txt";

	if(alpha.size() > 0) alpha.print(alphaPath);
  	if( beta.size() > 0)  beta.print( betaPath);
	if(gamma.size() > 0)  gamma.print(gammaPath);
	if(  eta.size() > 0)    eta.print(  etaPath);
	if(    s.size() > 0)     s.print(    sPath);
	if(    u.size() > 0)     u.print(    uPath);
}

template<typename Int , typename Real >
void pmswitch::InferenceData<Int, Real>::printLatents(std::string baseDir) const {
	std::string EZPath = baseDir + "/EZ.txt";
	std::string ESPath = baseDir + "/ES.txt";
	std::string EYPath = baseDir + "/EY.txt";

	if(EZ.size() > 0) EZ.print(EZPath);
	if(ES.size() > 0) ES.print(ESPath);
	if(EY.size() > 0) EY.print(EYPath);
}

template<typename Int , typename Real >
void pmswitch::InferenceData<Int, Real>::printDB(std::string baseDir) const {
	std::string gPath = baseDir + "/g.txt";
	if(g.size() > 0) g.print(gPath);
}

template<typename Int , typename Real >
void pmswitch::InferenceData<Int, Real>::printAllData(std::string baseDir) const {
	printHyperParameters(baseDir);
	printLatents(baseDir);
	printDB(baseDir);
}



/*
	public functions
*/





/*
	private functions
*/


#endif
