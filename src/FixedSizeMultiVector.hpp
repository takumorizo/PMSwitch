#ifndef FixedSizeMultiVector_h
#define FixedSizeMultiVector_h

#include "ExceptionUtil.hpp"
#include "StringUtil.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>
#include <limits>

namespace pmswitch{
	template<typename E, typename Int = long long >
	class FixedSizeMultiVector{
	public:
		FixedSizeMultiVector();
		FixedSizeMultiVector(std::vector<Int> dims);
		FixedSizeMultiVector(E e, std::vector<Int> dims);
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

		void print(std::ostream& out = std::cerr, std::string (*eToStr)(E) = NULL) const;
		void print(std::string path, std::string (*eToStr)(E) = NULL) const;

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

		bool isReal() const;
		bool isInt()  const;
	};

	template<typename E, typename Int = long long >
	class FixedSizeMultiVectorCreator{
	public:
		FixedSizeMultiVectorCreator();
		static FixedSizeMultiVector<E, Int> createFixedSizeMultiVector(std::string path, E (*strToE)(std::string));
		static FixedSizeMultiVector<E, Int> createFixedSizeMultiVector(std::string path, E defaultValue,  E (*strToE)(std::string));

		static FixedSizeMultiVector<E, Int> createFixedSizeMultiVector(std::string path, E (*strToE)(const std::string&));
		static FixedSizeMultiVector<E, Int> createFixedSizeMultiVector(std::string path, E defaultValue,  E (*strToE)(const std::string&));

	private:
		static void parseHeaderFile(std::string path, Int& size, Int& dim, std::vector<Int> &dimVec, std::vector<Int> &volumeVec);
		static void updateVector(std::string path, FixedSizeMultiVector<E, Int> &ans, E (*strToE)(std::string), bool checkNum = true);
		static void updateVector(std::string path, FixedSizeMultiVector<E, Int> &ans, E (*strToE)(const std::string&), bool checkNum = true);

	};
}
/*
	public functions
*/

template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E,Int>::FixedSizeMultiVector(){
}

template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E,Int>::FixedSizeMultiVector(std::vector<Int> dims){
	{//dimVec
		dimVec = dims;
	}
	{//sizeVec 4,3,2 -> 4*3*2, 3*2, 2 -> 4*3*2, 3*2, 2, 1 -> 3*2, 2, 1
		sizeVec = dims;
		sizeVec.push_back((Int)1);
		for(Int i = sizeVec.size()-1; i-1 >=0; i-- ){ sizeVec[i-1] *= sizeVec[i]; }
		sizeVec.erase(sizeVec.begin());
	}
	{//vec
		Int total = 1;
		for(auto d : dimVec){ total *= d;}
		vec = std::vector<E>(total);
	}
}

template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E,Int>::FixedSizeMultiVector(E e, std::vector<Int> dims){
	{//dimVec
		dimVec = dims;
	}
	{//sizeVec 4,3,2 -> 4*3*2, 3*2, 2 -> 4*3*2, 3*2, 2, 1 -> 3*2, 2, 1
		sizeVec = dims;
		sizeVec.push_back((Int)1);
		for(Int i = sizeVec.size()-1; i-1 >=0; i-- ){ sizeVec[i-1] *= sizeVec[i]; }
		sizeVec.erase(sizeVec.begin());
	}
	{//vec
		Int total = 1;
		for(auto d : dimVec){ total *= d;}
		vec = std::vector<E>(total, e);
	}
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
	// if(dimVec.size() != sizeof...(args)){
	// 	std::cerr << "dimVec.size(): " << dimVec.size() << ", " << sizeof...(args) << std::endl;
	// 	for(auto d : dimVec)std::cerr << "d: " << d << std::endl;
	// }
	assert(dimVec.size() == sizeof...(args));
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
bool pmswitch::FixedSizeMultiVector<E,Int>::isReal() const {
	return
		(typeid(E) == typeid(float)) ||
		(typeid(E) == typeid(double)) ||
		(typeid(E) == typeid(long double));
}


template<typename E, typename Int >
bool pmswitch::FixedSizeMultiVector<E,Int>::isInt() const {
	return
		(typeid(E) == typeid(char)) ||
		(typeid(E) == typeid(signed char)) ||
		(typeid(E) == typeid(unsigned char)) ||
		(typeid(E) == typeid(wchar_t)) ||
		(typeid(E) == typeid(char16_t)) ||
		(typeid(E) == typeid(char32_t)) ||
		(typeid(E) == typeid(short)) ||
		(typeid(E) == typeid(unsigned short)) ||
		(typeid(E) == typeid(int)) ||
		(typeid(E) == typeid(unsigned int)) ||
		(typeid(E) == typeid(long)) ||
		(typeid(E) == typeid(unsigned long)) ||
		(typeid(E) == typeid(long long)) ||
		(typeid(E) == typeid(unsigned long long));
}

template<typename E, typename Int >
void pmswitch::FixedSizeMultiVector<E,Int>::print(std::ostream& out, std::string (*eToStr)(E)) const {
	const std::streamsize ss = out.precision();
   	out << "#size\t" << vec.size()    << std::endl;
   	out << "#dim\t"  << dimVec.size() << std::endl;

   	out << "#dimVec\t";
   	if(dimVec.size() > 0){
	   	for(Int i = 0; i < dimVec.size()-1; i++){ out << dimVec[i] << ",";}
		out << dimVec[dimVec.size()-1] << std::endl;;
   	}else { out << std::endl; }

   	out << "#volumeVec\t";
	if(sizeVec.size() > 0){
	   	for(Int i = 0; i < sizeVec.size()-1; i++){ out << sizeVec[i] << ",";}
		out << sizeVec[sizeVec.size()-1] << std::endl;;
   	}else { out << std::endl; }

	for(Int i = 0; i < vec.size(); i++){
		out << i / volume(0);
		Int at = i % volume(0);
		for(Int d = 1; d < dimVec.size(); d++){
			out << "," << at / volume(d);
			at = i % volume(d);
		}
		out << "\t";
		if(    eToStr != NULL){
			out << eToStr(vec[i]) << std::endl;
		}else if(eToStr == NULL){
			if(isReal()){ out << std::setprecision( std::numeric_limits<long double>::max_digits10);}
			out <<        vec[i]  << std::endl;
			out << std::setprecision( ss );
		}
	}
}

template<typename E, typename Int >
void pmswitch::FixedSizeMultiVector<E,Int>::print(std::string path, std::string (*eToStr)(E)) const {
	using namespace pmswitch;
	std::ofstream ofs(path);
   	std::string buffer; std::vector<std::string> record;
   	if(ofs.fail()) {die("cannot open output file");}

   	ofs << "#size\t" << vec.size()    << std::endl;
   	ofs << "#dim\t"  << dimVec.size() << std::endl;

   	ofs << "#dimVec\t";
   	if(dimVec.size() > 0){
	   	for(Int i = 0; i < dimVec.size()-1; i++){ ofs << dimVec[i] << ",";}
		ofs << dimVec[dimVec.size()-1] << std::endl;;
   	}else { ofs << std::endl; }

   	ofs << "#volumeVec\t";
	if(sizeVec.size() > 0){
	   	for(Int i = 0; i < sizeVec.size()-1; i++){ ofs << sizeVec[i] << ",";}
		ofs << sizeVec[sizeVec.size()-1] << std::endl;;
   	}else { ofs << std::endl; }

	for(Int i = 0; i < vec.size(); i++){
		ofs << i / volume(0);
		Int at = i % volume(0);
		for(Int d = 1; d < dimVec.size(); d++){
			ofs << "," << at / volume(d);
			at = i % volume(d);
		}
		ofs << "\t";
		if(    eToStr != NULL){
			ofs << eToStr(vec[i]) << std::endl;
		}else if(eToStr == NULL){
			if(isReal()){ ofs << std::setprecision( std::numeric_limits<long double>::max_digits10);}
			ofs <<        vec[i]  << std::endl;
		}
	}
   	ofs.close();
}


template<typename E, typename Int >
void pmswitch::FixedSizeMultiVectorCreator<E, Int>::parseHeaderFile(std::string path, Int& size, Int& dim, std::vector<Int> &dimVec, std::vector<Int> &volumeVec){
	using namespace pmswitch;
	std::ifstream ifs(path);
   	std::string buffer; std::vector<std::string> record;
   	if(ifs.fail()) {die("cannot open output file");}

   	dimVec.clear(); volumeVec.clear();

   	unsigned int cheked = 0;
	for(int i = 0; i < 4; i++){
   		std::getline (ifs,buffer);
   		record = string::split(buffer, '\t');
   		if(record.size() != 2) {
   			die("invalid format header files @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::parseHeaderFile");
   		}
   		if(      record[0] == "#size"){
   			size = (Int)std::stoull(record[1]);
   			cheked |= (1 << 0);
   		}else if(record[0] == "#dim"){
   			dim  =(Int)std::stoull(record[1]);
   			cheked |= (1 << 1);
   		}else if(record[0] == "#dimVec"){
   			record = string::split(record[1], ',');
   			for(auto strD : record){ dimVec.push_back( (Int)std::stoull(strD) ); }
			cheked |= (1 << 2);
   		}else if(record[0] == "#volumeVec"){
   			record = string::split(record[1], ',');
   			for(auto strV : record){ volumeVec.push_back( (Int)std::stoull(strV) ); }
			cheked |= (1 << 3);
   		}else{
   			die("unknown header cols @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::parseHeaderFile");
   		}
   	}
   	ifs.close();

   	if(dimVec.size() != dim || volumeVec.size() != dim){
   		std::cerr << dim  << ", " << dimVec.size() << ", " << volumeVec.size() << std::endl;
   		die("dim not consistent @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::parseHeaderFile");
   	}
   	if(cheked != 15) die("insufficient header @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::parseHeaderFile");
}


template<typename E, typename Int >
void pmswitch::FixedSizeMultiVectorCreator<E, Int>::updateVector(std::string path, FixedSizeMultiVector<E, Int> &ans, E (*strToE)(std::string), bool checkNum){
	std::ifstream ifs(path);
   	std::string buffer; std::vector<std::string> record;
   	if(ifs.fail()) {die("cannot open output file");}

   	Int inputLine = 0;
	for(int i = 0; i < 4; i++){std::getline (ifs,buffer);} // pass header
   	while(std::getline (ifs,buffer)){
   		record = string::split(buffer, '\t');
   		if(record.size() != 2) {
   			die("invalid content format @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::updateVector");
   		}
   		Int at = 0;
   		E value = strToE(record[1]);
   		record = string::split(record[0], ',');
   		for(Int i = 0; i < record.size(); i++){
   			Int d = (Int)std::stoull(record[i]);
   			at += d * ans.volume(i);
   		}
   		ans[at] = value;
   		inputLine++;
   	}
   	ifs.close();

   	if(checkNum and (inputLine != ans.size())){ die("line number is not same as size @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::updateVector");}
}

template<typename E, typename Int >
void pmswitch::FixedSizeMultiVectorCreator<E, Int>::updateVector(std::string path, FixedSizeMultiVector<E, Int> &ans, E (*strToE)(const std::string&), bool checkNum){
	std::ifstream ifs(path);
   	std::string buffer; std::vector<std::string> record;
   	if(ifs.fail()) {die("cannot open output file");}

   	Int inputLine = 0;
	for(int i = 0; i < 4; i++){std::getline (ifs,buffer);} // pass header
   	while(std::getline (ifs,buffer)){
   		record = string::split(buffer, '\t');
   		if(record.size() != 2) {
   			die("invalid content format @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::updateVector");
   		}
   		Int at = 0;
   		E value = strToE(record[1]);
   		record = string::split(record[0], ',');
   		for(Int i = 0; i < record.size(); i++){
   			Int d = (Int)std::stoull(record[i]);
   			at += d * ans.volume(i);
   		}
   		ans[at] = value;
   		inputLine++;
   	}
   	ifs.close();

   	if(checkNum and (inputLine != ans.size())){ die("line number is not same as size @ pmswitch::FixedSizeMultiVectorCreator<E, Int>::updateVector");}
}


template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E, Int> pmswitch::FixedSizeMultiVectorCreator<E, Int>::createFixedSizeMultiVector(std::string path, E (*strToE)(std::string)){
	using namespace pmswitch;

   	Int size; Int dim;
   	std::vector<Int> dimVec;
	std::vector<Int> volumeVec;
	{
		parseHeaderFile(path, size, dim, dimVec, volumeVec);
	}

   	FixedSizeMultiVector<E, Int> ans(dimVec);
   	{
   		updateVector(path, ans, strToE, true);
   	}
   	return ans;
}

template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E, Int> pmswitch::FixedSizeMultiVectorCreator<E, Int>::createFixedSizeMultiVector(std::string path, E defaultValue, E (*strToE)(std::string)){
	using namespace pmswitch;

   	Int size; Int dim;
   	std::vector<Int> dimVec;
	std::vector<Int> volumeVec;
	{
		parseHeaderFile(path, size, dim, dimVec, volumeVec);
	}

   	FixedSizeMultiVector<E, Int> ans(defaultValue, dimVec);
   	{
   		updateVector(path, ans, strToE, false);
   	}
   	return ans;
}

template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E, Int> pmswitch::FixedSizeMultiVectorCreator<E, Int>::createFixedSizeMultiVector(std::string path, E (*strToE)(const std::string&)){
	using namespace pmswitch;

   	Int size; Int dim;
   	std::vector<Int> dimVec;
	std::vector<Int> volumeVec;
	{
		parseHeaderFile(path, size, dim, dimVec, volumeVec);
	}

   	FixedSizeMultiVector<E, Int> ans(dimVec);
   	{
   		updateVector(path, ans, strToE, true);
   	}
   	return ans;
}

template<typename E, typename Int >
pmswitch::FixedSizeMultiVector<E, Int> pmswitch::FixedSizeMultiVectorCreator<E, Int>::createFixedSizeMultiVector(std::string path, E defaultValue, E (*strToE)(const std::string&)){
	using namespace pmswitch;

   	Int size; Int dim;
   	std::vector<Int> dimVec;
	std::vector<Int> volumeVec;
	{
		parseHeaderFile(path, size, dim, dimVec, volumeVec);
	}

   	FixedSizeMultiVector<E, Int> ans(defaultValue, dimVec);
   	{
   		updateVector(path, ans, strToE, false);
   	}
   	return ans;
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
