// #ifndef MultiVector_h
// #define MultiVector_h
// #include <vector>
// #include <iostream>

// // template <typename E>
// // void set(std::vector<E> &list) {
// // }

// // template <typename E, class... Args>
// // void set(std::vector<E> &list, E e, Args... args){
// //   list.push_back(e);
// //   set(list, args...);
// // }

// // template <typename Int>
// // Int atImpl(Int Result, const std::vector<Int> &list, Int at) {
// // 	return Result;
// // }

// // template <typename Int, class... Args>
// // Int atImpl(Int Result, const std::vector<Int> &list, Int at, Int e, Args... args){
// // 	return atImpl(Result + e*list[at], list, at+1, args...);
// // }

// // x*d1       + y                  -> 1次元 :  d1,       2次元 : 1
// // x*d1*d2    + y*d2 + z           -> 1次元 :  d1*d2,    2次元 : d2,    3次元 : 1
// // x*d1*d2*d3 + y*d2*d3 + d3*z + w -> 1次元 :  d1*d2*d3, 2次元 : d2*d3, 3次元 : d3, 4次元 : 1

// template<typename E, typename Int = long long, class... Dims >
// class MultiVector{
// public:
// 	MultiVector(Int D, Int s, E e, Dims... dims){
// 		//dimVec
// 		{
// 			set<Int, Dims...>(dimVec, dims...);
// 			std::cerr << " dimvec after set " << std::endl;
// 			for(auto a : dimVec){ std::cerr << a << std::endl; }
// 			dimVec.push_back((Int)1);
// 			for(Int i = dimVec.size()-1; i--; i >=0){ dimVec[i] *= dimVec[i+1]; }
// 			dimVec.erase(dimVec.begin());
// 		}
// 	};
// 	template <class... Args> Int at(Args... args) const;
// 	template <class... Args> const E& operator()(Args... args) const;
// 	template <class... Args> E& operator()(Args... args);

// 	const E& operator()(Int i) const;
// 	E& operator()(Int i);

// 	void push_back(E e);
// 	Int size();
// 	template <class... Args > Int size(Int dim, Args... args); // TODO
// private:
// 	std::vector<E> vec;
// 	std::vector<Int> dimVec;

// 	template <typename F>
// 	static void set(std::vector<F> &list){
// 	}

// 	template <typename F, class... Args>
// 	static void set(std::vector<F> &list, F e, Args... args){
// 		  list.push_back(e);
// 		  set(list, args...);
// 	}

// 	static Int atImpl(Int Result, const std::vector<Int> &list, Int at) {
// 		return Result;
// 	}

// 	template <class... Args >
// 	static Int atImpl(Int Result, const std::vector<Int> &list, Int at, Int e, Args... args){
// 		//assert(sizeof...(args) == len(list) - at  )
// 		return atImpl(Result + e*list[at], list, at+1, args...);
// 	}

// 	Int sizeImpl(Int res1, Int res2, Int rnum, Int dim){
// 		rnum++;
// 		return  vec.size()/dimVec[rnum] + 1;
// 	}

// 	Int sizeImpl(Int res1, Int res2, Int rnum, Int dim, Int pos){
// 		res1 += dimVec[rnum] * (pos+1);
// 		res2 += dimVec[rnum] *  pos   ;
// 		rnum ++;
// 		return (min( vec.size(), res1 ) - res2 )/dimVec[rnum] + 1;
// 	}

// 	template <class... Args >
// 	Int sizeImpl(Int res1, Int res2, Int rnum, Int dim, Int pos, Args... args){
// 		res1 += dimVec[rnum] *  pos;
// 		res2 += dimVec[rnum] *  pos;
// 		return sizeImpl(res1, res2, rnum, dim, args...);
// 	}
// 	// Int size(Int dim, Args... args){
// 	// 	return sizeImpl(0, 0, 0, dim, args...);
// 	// }

// 	/*
// 		Dims: d1, d2, d3
// 		dimVec: {d2*d3, d3, 1}

// 		sizes @(x,y) : {size/d1+ 1, size/d1+ 1}
// 		parse vector in this way
// 		for(x in size(0))for(y in size(1,x))for(z in size(2,x,y)){
// 			multiVec(x,y,z)
// 		}
// 		x : 0:{  vec.size()/(d1*d2) + 1 }
// 		y : 0:{ (min( vec.size(), d1*d2*(x+1)) - d1*d2*x)/d2 + 1 }
// 		z : 0:{ (min( vec.size(), d1*d2*x + d2*(y+1) ) - d1*d2*x + d2*y )/1 + 1 }

// 		x : 0:{  vec.size()/(d1*d2) + 1 }
// 		y : 0:{ (min( vec.size(), d1*d2*(x+1)) - d1*d2*x)/d2 + 1 }
// 		z : 0:{ (min( vec.size(), d1*d2*x + d2*(y+1) ) - d1*d2*x - d2*y )/1 + 1 }

// 		Int sizeImpl(Int res1, Int res2, Int rnum, Int dim){
// 			rnum++;
// 			return  vec.size()/dimList[rnum] + 1;
// 		}

// 		Int sizeImpl(Int res1, Int res2, Int rnum, Int dim, Int pos){
// 			res1 += dimList[rnum] * (pos+1);
// 			res2 += dimList[rnum] *  pos   ;
// 			rnum ++;
// 			return (min( vec.size(), res1 ) - res2 )/dimList[rnum] + 1
// 		}

// 		Int sizeImpl(Int res1, Int res2, Int rnum, Int dim, Int pos, Args... args){
// 			res1 += dimList[rnum] *  pos;
// 			res2 += dimList[rnum] *  pos;
// 			return sizeImpl(Int res1, Int res2, Int rnum, Int dim, Args... args);
// 		}
// 		Int size(Int dim, Args... args){
// 			return sizeImpl(0, 0, 0, dim, args...);
// 		}
// 	*/

// };

// template<typename E, typename Int, class... Dims >
// template <class... Args>
// Int MultiVector<E,Int,Dims... >::at(Args... args) const {
// 	return atImpl<args... >(0, dimVec, 0, args...);
// }


// template<typename E, typename Int, class... Dims >
// template <class... Args>
// const E& MultiVector<E,Int,Dims... >::operator()(Args... args) const {
// 	return vec[at(args...)];
// }

// template<typename E, typename Int, class... Dims >
// template <class... Args>
// E& MultiVector<E,Int,Dims... >::operator()(Args... args) {
// 	return vec[at(args...)];
// }

// template<typename E, typename Int, class... Dims >
// const E& MultiVector<E,Int,Dims... >::operator()(Int i) const {
// 	return vec[i];
// }

// template<typename E, typename Int, class... Dims >
// E& MultiVector<E,Int,Dims... >::operator()(Int i) {
// 	return vec[i];
// }

// template<typename E, typename Int, class... Dims >
// void MultiVector<E,Int,Dims... >::push_back(E e){
// 	vec.push_back(e);
// }

// template<typename E, typename Int, class... Dims >
// Int MultiVector<E,Int,Dims... >::size(){
// 	return (Int)(vec.size());
// }

// template<typename E, typename Int, class... Dims >
// template <class... Args>
// Int MultiVector<E,Int,Dims... >::size(Int dim, Args... args){
// 	return sizeImpl<args... >(0, 0, 0, dim, args...);
// }


// #endif
