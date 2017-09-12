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
					  FixedSizeMultiVector<Real, Int> _g, Real _Lq = 0.0);

		void printHyperParameters(std::string baseDir) const;
		void printLatents(std::string baseDir)		   const;
		void printDB(std::string baseDir)              const;
		void printAllData(std::string baseDir)         const;

		FixedSizeMultiVector<Real, Int> getAlpha() const {return alpha;}
		FixedSizeMultiVector<Real, Int> getBeta()  const {return beta;}
		FixedSizeMultiVector<Real, Int> getGamma() const {return gamma;}
		FixedSizeMultiVector<Real, Int> getEta()   const {return eta;}
		FixedSizeMultiVector<Real, Int> getG()     const {return g;}

		FixedSizeMultiVector<Real, Int> getEZ()    const {return EZ;}
		FixedSizeMultiVector<Real, Int> getES()    const {return ES;}
		FixedSizeMultiVector<Real, Int> getEY()    const {return EY;}

		Real getLq() const {return Lq;}

	private:
		FixedSizeMultiVector<Real, Int> EZ;
		FixedSizeMultiVector<Real, Int> ES;
		FixedSizeMultiVector<Real, Int> EY;

		FixedSizeMultiVector<Real, Int> alpha;
		FixedSizeMultiVector<Real, Int> eta;
		FixedSizeMultiVector<Real, Int> beta;
		FixedSizeMultiVector<Real, Int> gamma;

		FixedSizeMultiVector<Real, Int> g;
		Real Lq;
	};
}

template<typename Int , typename Real >
pmswitch::InferenceData<Int, Real>::InferenceData(FixedSizeMultiVector<Real, Int> _EZ, FixedSizeMultiVector<Real, Int> _ES, FixedSizeMultiVector<Real, Int> _EY,
					  FixedSizeMultiVector<Real, Int> _alpha, FixedSizeMultiVector<Real, Int> _eta,
					  FixedSizeMultiVector<Real, Int> _beta, FixedSizeMultiVector<Real, Int> _gamma,
					  FixedSizeMultiVector<Real, Int> _g, Real _Lq) : EZ(_EZ), ES(_ES), EY(_EY), alpha(_alpha), eta(_eta), beta(_beta), gamma(_gamma), g(_g), Lq(_Lq) {
}

template<typename Int , typename Real >
void pmswitch::InferenceData<Int, Real>::printHyperParameters(std::string baseDir) const {
	std::string alphaPath = baseDir + "/alpha.txt";
	std::string betaPath  = baseDir + "/beta.txt";
	std::string gammaPath = baseDir + "/gamma.txt";
	std::string etaPath   = baseDir + "/eta.txt";
	alpha.print(alphaPath);
	beta.print(betaPath);
	gamma.print(gammaPath);
	eta.print(etaPath);
}

template<typename Int , typename Real >
void pmswitch::InferenceData<Int, Real>::printLatents(std::string baseDir) const {
	std::string EZPath = baseDir + "/EZ.txt";
	std::string ESPath = baseDir + "/ES.txt";
	std::string EYPath = baseDir + "/EY.txt";
	EZ.print(EZPath);
	ES.print(ESPath);
	EY.print(EYPath);
}

template<typename Int , typename Real >
void pmswitch::InferenceData<Int, Real>::printDB(std::string baseDir) const {
	std::string gPath = baseDir + "/g.txt";
	g.print(gPath);
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
