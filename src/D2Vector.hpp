#ifndef D2Vector_h
#define D2Vector_h
#include <vector>
#include <iostream>


template<typename E, typename Int = int>
class D2Vector{
public:
	D2Vector(Int D, Int s, E e){
		vec = std::vector<E>(s,e);
		d1 = D;
	};
	Int at2(Int x, Int y){ return x * d1 + y; };

	const E& operator[](Int i) const;
	E& operator[](Int i);

	const E& operator()(Int i) const;
	E& operator()(Int i);

	const E& operator()(Int x, Int y) const;
	E& operator()(Int x, Int y);

	void push_back(E e);
	void push_back(E e, Int dim);
	Int size();
	Int size(Int dim);

private:
	std::vector<E> vec;
	Int d1;
};

template<typename E, typename Int>
const E& D2Vector<E,Int>::operator[](Int i) const {
	return vec[i];
}

template<typename E, typename Int>
E& D2Vector<E,Int>::operator[](Int i) {
	return vec[i];
}

template<typename E, typename Int>
const E& D2Vector<E,Int>::operator()(Int i) const {
	return vec.at(i);
}

template<typename E, typename Int>
E& D2Vector<E,Int>::operator()(Int i) {
	return vec.at(i);
}

template<typename E, typename Int>
const E& D2Vector<E,Int>::operator()(Int x, Int y) const {
	return vec[ x*d1 + y ];
}

template<typename E, typename Int>
E& D2Vector<E,Int>::operator()(Int x, Int y){
	return vec[ x*d1 + y ];
}

template<typename E, typename Int>
void D2Vector<E,Int>::push_back(E e){
	vec.push_back(e);
}

template<typename E, typename Int>
Int D2Vector<E,Int>::size(){
	return (Int)vec.size();
}

template<typename E, typename Int>
Int D2Vector<E,Int>::size(Int dim){
	if(dim == 0){
		return std::min(vec.size(), d1);
	}else if(dim == 1){
		return vec.size() % d1;
	}else{
		std::cout << "out of dimensional push_back" << std::endl;
		exit(1);
		return -1;
	}
}

template<typename E, typename Int>
void D2Vector<E,Int>::push_back(E e, Int dim){
	if(dim == 0){
		vec.push_back(e);
	}else if(dim == 1){
		for(int i = 0; i < d1; i++) vec.push_back(e);
	}else{
		std::cout << "out of dimensional push_back" << std::endl;
		exit(1);
	}
}

#endif
