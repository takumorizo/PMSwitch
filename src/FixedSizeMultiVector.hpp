#ifndef FixedSizeMultiVector_h
#define FixedSizeMultiVector_h
#include <vector>
#include <iostream>
#include <assert.h>

namespace pmswitch{
	template<typename E, typename Int = long long >
	class FixedSizeMultiVector{
	public:
		FixedSizeMultiVector();
		template<class... Dims >
		FixedSizeMultiVector(E e, Dims... dims);
		template <class... Args> Int at(Args... args) const;
		template <class... Args> typename std::vector<E>::const_reference operator()(Args... args) const;
		template <class... Args> typename std::vector<E>::reference operator()(Args... args);
		typename std::vector<E>::const_reference operator[](Int idx) const;
		typename std::vector<E>::reference operator[](Int idx);

		Int size() const;
		Int size(Int dim) const;
		Int volume(Int dim) const;
		Int dim() const;
		bool operator==(const FixedSizeMultiVector<E, Int>& rhs) const;

		void print() const;

	private:
		/*
			Dims: d1, d2, d3
			sizeVec: {d2*d3, d3, 1}
			dimVec: {d1, d2, d3}
		*/
		std::vector<E> vec;
		std::vector<Int> sizeVec;
		std::vector<Int> dimVec;

		template <typename F>
		void setImpl(std::vector<F> &list);

		template <typename F, class... Args>
		void setImpl(std::vector<F> &list, F e, Args... args);

		template <typename F, class... Args>
		void set(std::vector<F> &list, Args... args);

		Int atImpl(Int Result, const std::vector<Int> &list, Int at) const;

		template <class... Args >
		Int atImpl(Int Result, const std::vector<Int> &list, Int at, Int e, Args... args) const;
	};

/*
	template<typename E, typename Int = long long, bool Forward = true >
	class Range {
	public:
	    class iterator {
	    public:
	        iterator( Int _num );
	        iterator& operator++(){ num = Forward ? num + 1: num - 1; return *this; }
	        iterator operator++(int) {iterator retval = *this; ++(*this); return retval;}
	        bool operator==(iterator other) const {return num == other.num;}
	        bool operator!=(iterator other) const {return !(*this == other);}
	        Int operator*() {return num;}
	        // iterator traits
	        using difference_type = Int;
	        using value_type = Int;
	        using pointer = const Int*;
	        using reference = const Int&;
	        using iterator_category = std::forward_iterator_tag;
	    private:
		    Int num  = 0;
	    };

	    iterator begin() {return from;}
	    iterator end() {return Forward? to+1 : to-1;}
	    Range( FixedSizeMultiVector<E, Int> v, std::vector<Int> _from, std::vector<Int> _to);
	private:
		Int from = 0;
		Int to = 0;
	};
*/
}

/*
	public functions
*/

template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E,Int>::FixedSizeMultiVector(){
}


template<typename E, typename Int >
template<class... Dims >
pmswitch::FixedSizeMultiVector<E,Int>::FixedSizeMultiVector(E e, Dims... dims){
	{//dimVec
		std::vector<Int> tmp;
		FixedSizeMultiVector<E,Int>::template set<Int, Dims...>(tmp, dims...);
		dimVec = tmp;
		// std::cerr << "dimVec" << std::endl;
		// for(auto a : dimVec){ std::cerr << a << std::endl; }
	}
	{//sizeVec 4,3,2 -> 4*3*2, 3*2, 2 -> 4*3*2, 3*2, 2, 1 -> 3*2, 2, 1
		FixedSizeMultiVector<E,Int>::template set<Int, Dims...>(sizeVec, dims...);
		sizeVec.push_back((Int)1);
		for(Int i = sizeVec.size()-1; i-1 >=0; i-- ){ sizeVec[i-1] *= sizeVec[i]; }
		sizeVec.erase(sizeVec.begin());
		// std::cerr << "sizeVec " << std::endl;
		// for(auto a : sizeVec){ std::cerr << a << std::endl; }
	}
	{//vec
		Int total = 1;
		for(auto d : dimVec){ total *= d;}
		vec = std::vector<E>(total, e);
	}
}

template<typename E, typename Int >
template <class... Args>
Int pmswitch::FixedSizeMultiVector<E,Int>::at(Args... args) const {
	// assert(len(args) == len(sizeVec))
	return atImpl(0, sizeVec, 0, args...);
}


template<typename E, typename Int >
template <class... Args>
typename std::vector<E>::const_reference pmswitch::FixedSizeMultiVector<E,Int>::operator()(Args... args) const {
	return vec[at(args...)];
	// #ifdef NDEBUG
	// 	return vec[at(args...)];
	// #else
	// 	return vec.at(at(args...));
	// #endif
}

template<typename E, typename Int >
template <class... Args>
typename std::vector<E>::reference pmswitch::FixedSizeMultiVector<E,Int>::operator()(Args... args) {
	return vec[at(args...)];
	// #ifdef NDEBUG
	// 	return vec[at(args...)];
	// #else
	// 	return vec.at(at(args...));
	// #endif
}

template<typename E, typename Int >
typename std::vector<E>::reference pmswitch::FixedSizeMultiVector<E,Int>::operator[](Int idx) {
	return vec[idx];
}


// template<typename E, typename Int >
// E& pmswitch::FixedSizeMultiVector<E,Int>::operator[](Int idx) {
// 	return vec[idx];
// }

template<typename E, typename Int >
typename std::vector<E>::const_reference pmswitch::FixedSizeMultiVector<E,Int>::operator[](Int idx) const {
	return vec[idx];
}

template<typename E, typename Int >
Int pmswitch::FixedSizeMultiVector<E,Int>::size() const {
	return (Int)(vec.size());
}

template<typename E, typename Int >
Int pmswitch::FixedSizeMultiVector<E,Int>::size(Int dim) const {
	// return dimVec[dim];
	#ifdef NDEBUG
		// std::cerr << "[] called" << std::endl;
		return dimVec[dim];
	#else
		// std::cerr << ".at called" << std::endl;
		return dimVec.at(dim);
	#endif
}

template<typename E, typename Int >
Int pmswitch::FixedSizeMultiVector<E,Int>::dim() const {
	return (Int)(dimVec.size());
}

template<typename E, typename Int >
Int pmswitch::FixedSizeMultiVector<E,Int>::volume(Int dim) const {
	#ifdef NDEBUG
		// std::cerr << "[] called" << std::endl;
		return (Int)(sizeVec[dim]);
	#else
		// std::cerr << ".at called" << std::endl;
		return (Int)(sizeVec.at(dim));
	#endif
}

template<typename E, typename Int >
void pmswitch::FixedSizeMultiVector<E,Int>::print() const {
	std::cerr << "dimVec" << std::endl;
	for(auto a : dimVec){ std::cerr << a << std::endl; }

	std::cerr << "sizeVec " << std::endl;
	for(auto a : sizeVec){ std::cerr << a << std::endl; }

	std::cerr << "vec " << std::endl;
	for(auto a : vec){ std::cerr << a << std::endl; }
}

template<typename E, typename Int >
bool pmswitch::FixedSizeMultiVector<E,Int>::operator==(const FixedSizeMultiVector<E, Int>& rhs) const{
	bool equalSizeVec = (sizeVec.size() == rhs.sizeVec.size() ) and ( std::equal(sizeVec.cbegin(), sizeVec.cend(), rhs.sizeVec.cbegin()) );
	bool equalDimVec  = ( dimVec.size() == rhs.dimVec.size()  ) and ( std::equal(dimVec.cbegin(),  dimVec.cend(),  rhs.dimVec.cbegin()) );
	if( not( equalSizeVec and equalDimVec )) return false;

	bool equalVec = (vec.size() == rhs.vec.size()) and ( std::equal(vec.cbegin(), vec.cend(), rhs.vec.cbegin()) );
	return (equalSizeVec and equalDimVec and equalVec);
}


/*
	private functions
*/
template<typename E, typename Int >
template <typename F>
void pmswitch::FixedSizeMultiVector<E,Int>::setImpl(std::vector<F> &list){
}

template<typename E, typename Int >
template <typename F, class... Args>
void pmswitch::FixedSizeMultiVector<E,Int>::setImpl(std::vector<F> &list, F e, Args... args){
	list.push_back(e);
	setImpl<F>(list, args...);
}

template<typename E, typename Int >
template <typename F, class... Args>
void pmswitch::FixedSizeMultiVector<E,Int>::set(std::vector<F> &list, Args... args){
	setImpl<F>(list, args...);
}

template<typename E, typename Int >
Int pmswitch::FixedSizeMultiVector<E,Int>::atImpl(Int Result, const std::vector<Int> &list, Int at) const {
	return Result;
}

template<typename E, typename Int >
template <class... Args >
Int pmswitch::FixedSizeMultiVector<E,Int>::atImpl(Int Result, const std::vector<Int> &list, Int at, Int e, Args... args) const{
	return atImpl(Result + e*list[at], list, at+1, args...);
}


#endif
