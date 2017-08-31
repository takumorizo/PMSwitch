#ifndef DBData_h
#define DBData_h


#include "FixedSizeMultiVector.hpp"
/*
	N     :: Int // db size
	L     :: Int // dimension of feature vector
	M_{l} :: Int // # of possible state at l-th col of feature vector
	maxMl :: Int
	g_{n, } :: m_{l}-dimensional non-negative simplex //g_{n,l} : m_{l} sized simplex
*/
namespace pmswitch{
	template<typename Int = int, typename Real = double>
	class DBData{
	public:
		DBData(Int N, Int L, Int maxMl, FixedSizeMultiVector<Int> Ml, FixedSizeMultiVector<Real> g);

		const Int N;
		const Int L;
		const Int maxMl;
		const FixedSizeMultiVector<Int> Ml;
		const FixedSizeMultiVector<Real> g;
	};
}



template<typename Int, typename Real>
pmswitch::DBData<Int, Real>::DBData(Int N,
	 					  Int L, Int maxMl, FixedSizeMultiVector<Int> Ml,
	 					  FixedSizeMultiVector<Real> g) : N(N), L(L), maxMl(maxMl), Ml(Ml), g(g){
}



#endif