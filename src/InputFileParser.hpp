#ifndef InputFileParser_h
#define InputFileParser_h

#include "FeatureData.hpp"
#include "DBData.hpp"
#include "ExceptionUtil.hpp"

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace pmswitch{
	template<typename Int = long long, typename Real = double >
	class InputFileParser{
	public:
		InputFileParser(){};
		pmswitch::FeatureData<Int, Real> parseFeatureFile(std::string pathFV);
		pmswitch::DBData<Int, Real> parseDBFile(std::string pathDB);
	private:
		std::vector<std::string> split(const std::string &str, char delimiter) const;

		void prepareParseFV(std::string pathFV,
							Int &maxMl, Int &L, FixedSizeMultiVector<Int> &ml,
							Int &I, Int &maxJ, FixedSizeMultiVector<Int> &Js, std::vector<Int> &idxToPosSize,
							std::unordered_map<Int, std::string> &idxToSampleName,
							std::unordered_map<std::string, Int> &sampleNameToIdx) const;
		pmswitch::FixedSizeMultiVector<Int> parseFVX(std::string pathFV, Int I, Int maxJ, Int L, Int maxMl,
			          const std::unordered_map<std::string, Int> &sampleNameToIdx) const;

		void getHeadInfoFV(const std::vector<std::string> &record,
						   Int &maxMl, Int &L, FixedSizeMultiVector<Int> &ml ) const;
		void updateSampleNameIndexRelationFV(const std::string &sampleName,
										   Int& I,
										   std::vector<Int> &idxToPosSize,
										   std::unordered_map<Int, std::string> &idxToSampleName,
										   std::unordered_map<std::string, Int> &sampleNameToIdx) const;
		void getJsInfoFV(const std::vector<Int>& idxToPosSize, Int I,
						 FixedSizeMultiVector<Int> &Js, Int &maxJ) const;
	};
}


/*
.db ファイルフォーマット
	N L m1 m2, ... mL
	n 1 g_{n,1,1}, ..., g_{n,1,m1}
	n . .
	n . .
	n L g_{n,L,1}, ..., g_{n,L,mL}

	N     :: Int // db size
	L     :: Int // dimension of feature vector
	M_{l} :: Int // # of possible state at l-th col of feature vector
	maxMl :: Int
	g_{n, } :: m_{l}-dimensional non-negative simplex //g_{n,l} : m_{l} sized simplex

.db ファイルフォーマット
	N L m1 m2, ... mL
	n 1 g_{n,1,1}, ..., g_{n,1,m1}
	n . .
	n . .
	n L g_{n,L,1}, ..., g_{n,L,mL}

template<typename Int = int, typename Real = double>
class DBData{
public:
	DBData(Int N, Int L, Int maxMl, FixedSizeMultiVector<Int> Ml, FixedSizeMultiVector<Real> g);

	const Int N;
	const Int L;
	const FixedSizeMultiVector<Int> Ml;
	const Int maxMl;
	const FixedSizeMultiVector<Real> g;
};
*/

/*
public
*/
template<typename Int, typename Real >
pmswitch::DBData<Int, Real > pmswitch::InputFileParser<Int, Real >::parseDBFile(std::string pathDB){
	using namespace pmswitch;
	std::ifstream ifs(pathDB);
   	std::string buffer; std::vector<std::string> record;
   	if(ifs.fail()) {die("cannot open .db file");}

	//header
	Int N, L, maxMl=0;
	FixedSizeMultiVector<Int>  ml;
	{
		std::getline (ifs,buffer); record = split(buffer, '\t');
		N = std::stoi(record[0]); L = std::stoi(record[1]);
		if(record.size() != 2 + L) { die("failed to read header of db file.");}
		FixedSizeMultiVector<Int>  tmp(0,L);
		for(Int l = 0; l < L; l++){ tmp(l) = std::stoi(record[2+l]); }
		ml = tmp;
		for(Int l = 0; l < L; l++){ maxMl = std::max(ml(l), maxMl); }
		buffer = ""; record.clear();
	}
	//content
	FixedSizeMultiVector<Real> g(0.0, N, L, maxMl);
	{
		for(Int n = 0; n < N; n++)for(Int l = 0; l < L; l++){
			std::getline (ifs,buffer); record = split(buffer, '\t');
			if(record.size() != 2 + ml(l)) { die("failed to read content of db file.");}
			for(Int ll = 0; ll < ml(l); ll++){ g(n,l,ll) = std::stod(record[2+ll]); }
			buffer = ""; record.clear();
		}
	}
	ifs.close();

	pmswitch::DBData<Int, Real > ans(N, L, maxMl, ml, g);
	return ans;
}

/*
.fv ファイルフォーマット
	L m1 ... mL //header
	sample, x1, ..., xL //eachSample

template<typename Int = long long, typename Real = double>
class FeatureData{
public:
	FeatureData(Int I, Int maxJ, std::vector<Int> Js,
				Int L, Int maxMl, FixedSizeMultiVector<Int> Ml,
				FixedSizeMultiVector<Int> X,
				std::unordered_map<Int, std::string> sampleIdxToSampleName);
	const Int I;
	const FixedSizeMultiVector<Int> Js;
	const Int maxJ;
	const Int maxMl;
	const FixedSizeMultiVector<Int> X;
	const Int L;
	const FixedSizeMultiVector<Int> Ml;

	const std::unordered_map<Int, std::string> sampleIdxToSampleName;
	// std::unordered_map<Int, Int> idxToPosSize;
};
*/
template<typename Int, typename Real >
pmswitch::FeatureData<Int, Real > pmswitch::InputFileParser<Int, Real >::parseFeatureFile(std::string pathFV){
	using namespace pmswitch;

	Int maxMl=0, L;
	FixedSizeMultiVector<Int> ml;
   	Int I=0, maxJ=0; FixedSizeMultiVector<Int> Js; std::vector<Int> idxToPosSize;
	std::unordered_map<Int, std::string> idxToSampleName;
	std::unordered_map<std::string, Int> sampleNameToIdx;
	{
		prepareParseFV(pathFV, maxMl, L, ml,
					   I, maxJ, Js, idxToPosSize,
					   idxToSampleName,
					   sampleNameToIdx);
	}

	FixedSizeMultiVector<Int> X(0,I,maxJ,L,maxMl);
	{
		X = parseFVX(pathFV, I, maxJ,  L, maxMl, sampleNameToIdx);
	}

	pmswitch::FeatureData<Int, Real > ans(I, maxJ, Js, L, maxMl, ml, X, idxToSampleName);
	return ans;
}


/*
private
*/
template<typename Int, typename Real >
std::vector<std::string> pmswitch::InputFileParser<Int, Real >::split(const std::string &str, char delimiter) const {
	std::vector<std::string> ans;
	std::istringstream sstr(str);
	std::string token;
	while (getline(sstr, token, delimiter)) ans.push_back(token);
	return ans;
}

template<typename Int, typename Real >
void pmswitch::InputFileParser<Int, Real >::prepareParseFV(std::string pathFV,
					Int &maxMl, Int &L, FixedSizeMultiVector<Int> &ml,
					Int &I, Int &maxJ, FixedSizeMultiVector<Int> &Js, std::vector<Int> &idxToPosSize,
					std::unordered_map<Int, std::string> &idxToSampleName,
					std::unordered_map<std::string, Int> &sampleNameToIdx) const{
	using namespace pmswitch;
	std::ifstream ifs(pathFV);
	if(ifs.fail()) {die("cannot open .fv file");}
   	std::string buffer; std::vector<std::string> record;
	// header
	maxMl=0;
	{// Int maxMl=0, L; FixedSizeMultiVector<Int> ml;
		std::getline (ifs,buffer);
		record = split(buffer, '\t');
		getHeadInfoFV(record, maxMl, L, ml );
		buffer = ""; record.clear();
	}
	// Supplemental infor for reading X;
	I=0; maxJ=0;
	{	// Int I=0, maxJ=0; FixedSizeMultiVector<Int> Js; std::vector<Int> idxToPosSize;
		// std::unordered_map<Int, std::string> idxToSampleName;
		// std::unordered_map<std::string, Int> sampleNameToIdx;
		while(std::getline (ifs,buffer)){
			record = split(buffer, '\t'); std::string sampleName = record[0];
			updateSampleNameIndexRelationFV(sampleName, I, idxToPosSize, idxToSampleName, sampleNameToIdx);
			buffer = ""; record.clear();
		}
		getJsInfoFV(idxToPosSize, I, Js, maxJ);
	}
	ifs.close();
	// std::cerr << "done prepare parse FV file. " << std::endl;
}


template<typename Int, typename Real >
pmswitch::FixedSizeMultiVector<Int> pmswitch::InputFileParser<Int, Real >::parseFVX(std::string pathFV, Int I, Int maxJ, Int L, Int maxMl,
													 const std::unordered_map<std::string, Int> &sampleNameToIdx) const{

	std::string buffer; std::vector<std::string> record;
	std::ifstream ifs(pathFV);
	if(ifs.fail()) {die("cannot open .fv file");}

	std::getline(ifs,buffer); buffer = "";
	FixedSizeMultiVector<Int> X(0,I,maxJ,L,maxMl);
	{
		FixedSizeMultiVector<Int> nowPos(0, I);
		while(std::getline(ifs, buffer)){
			std::vector<Int> Xs(L, 0);
			record = split(buffer, '\t');
			std::string sampleName = record[0];
			for(Int l = 0; l < L; l++){ Xs[l] = std::stoi(record[l+1]); }
			Int idx = sampleNameToIdx.at(sampleName);
			Int pos = nowPos(idx);
			for(Int l = 0; l < L; l++){
				X(idx, pos, l, Xs[l]) = 1;
			}
			nowPos(idx) += 1;

			buffer = ""; record.clear();
		}
	}
	ifs.close();
	return X;
}

template<typename Int, typename Real >
void pmswitch::InputFileParser<Int, Real >::getHeadInfoFV(const std::vector<std::string> &record,
						   Int &maxMl, Int &L, FixedSizeMultiVector<Int> &ml ) const{
	L = std::stoi(record[0]);
	if(record.size() != L + 1){die("failed to read header of FV file.");}
	FixedSizeMultiVector<Int> tmp(0, L);
	for(Int l = 0; l < L; l++){ tmp(l) = std::stoi(record[l+1]); }
	ml = tmp;
	for(Int l = 0; l < L; l++){ maxMl = std::max(maxMl, ml(l));}
}

template<typename Int, typename Real >
void pmswitch::InputFileParser<Int, Real >::updateSampleNameIndexRelationFV(
											const std::string &sampleName,
											Int& I,
											std::vector<Int> &idxToPosSize,
											std::unordered_map<Int, std::string> &idxToSampleName,
											std::unordered_map<std::string, Int> &sampleNameToIdx) const {
	if(sampleNameToIdx.count(sampleName) == 0){
		sampleNameToIdx[sampleName] = I;
		idxToSampleName[I] = sampleName;
		I++;
		idxToPosSize.push_back(1);
	}else{
		Int idx = sampleNameToIdx[sampleName];
		idxToPosSize[idx] += 1;
	}
}

template<typename Int, typename Real >
void pmswitch::InputFileParser<Int, Real >::getJsInfoFV(const std::vector<Int>& idxToPosSize, Int I,
														FixedSizeMultiVector<Int> &Js, Int &maxJ) const{
	FixedSizeMultiVector<Int> tmpJs(0,I);
	for(Int i = 0; i < I; i++){ tmpJs(i) = idxToPosSize[i];}
	Js = tmpJs;

	for(Int i = 0; i < I; i++){ maxJ = *std::max_element(idxToPosSize.begin(), idxToPosSize.end());}
}


#endif