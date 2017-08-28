#include "Inference.hpp"
#include "InputFileParser.hpp"
#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

std::string data_path;

// The fixture for testing class.
class PMSwithTest : public ::testing::Test {
protected:
    // You can remove any or all of the following functions if its body
    // is empty.
    PMSwithTest() {
        // You can do set-up work for each test here.
    }
    virtual ~PMSwithTest() {
        // You can do clean-up work that doesn't throw exceptions here.
    }
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:
    virtual void SetUp() {
        // Code here will be called immediately after the constructor (right
        // before each test).
    }
    virtual void TearDown() {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }
    // Objects declared here can be used by all tests in the test case for Foo.
};

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
TEST_F(PMSwithTest, InputParserDB) {
    pmswitch::InputFileParser<long long, double> parser;
    long long N = 10;
    long long L = 4;
    std::vector<long long> ml = {2,3,4,5};
    pmswitch::FixedSizeMultiVector<double> gOriginal(0,10,4,5);
    for(long long n = 0; n < N; n++)for(long long l = 0; l < L; l++)for(long long ll = 0; ll < ml[l]; ll++){
        gOriginal(n,l,ll) = 0.1;
    }
    std::string testDBPath = data_path + "/test.db";

    std::ofstream ofs(testDBPath);
    ofs << N << '\t' << L;
    for(long long l = 0; l < L; l++) ofs << '\t' << ml[l];
    ofs << std::endl;

    for(long long n = 0; n < N; n++){
        ofs << n;
        for(long long l = 0; l < L; l++){
            ofs << '\t' << l;
            for(long long ll = 0; ll < ml[l]; ll++){ ofs << '\t' << gOriginal(n,l,ll); }
            ofs << std::endl;
        }
    }
    // std::cout << "stoped output files" << std::endl;
    ofs.close();
    pmswitch::DBData<long long, double> dbData = parser.parseDBFile(testDBPath);
    // {
    //     gOriginal.print();
    //     dbData.g.print();
    // }

    EXPECT_TRUE(N == dbData.N);
    EXPECT_TRUE(L == dbData.L);

    EXPECT_TRUE(gOriginal == dbData.g);
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

TEST_F(PMSwithTest, InputParserFV) {
    using namespace pmswitch;
    InputFileParser<long long, double> parser;
    long long sampleSize = 5;
    std::vector<long long> posSizeEach = {2,2,3,3,4};
    long long maxJ = *std::max_element(posSizeEach.begin(), posSizeEach.end());
    long long L = 4;
    std::vector<long long> ml = {2,3,4,5};
    long long maxMl = *std::max_element(ml.begin(), ml.end());
    FixedSizeMultiVector<long long> xOriginal(0, sampleSize, maxJ, L, maxMl);
    FixedSizeMultiVector<long long> xOutput(0, sampleSize, maxJ, L);
    srand(0);
    for(long long n = 0; n < sampleSize; n++)for(long long p = 0; p < posSizeEach[n]; p++)for(long long l = 0; l < L; l++){
        long long value = rand() % ml[l];
        xOriginal(n, p, l, value ) = 1; xOutput(n,p,l) = value;
    }

    std::string testFVPath = data_path + "/test.fv";
    std::ofstream ofs(testFVPath);
    // headers
    ofs << L;
    for(long long l = 0; l < L; l++) { ofs << '\t' << ml[l];}
    ofs << std::endl;

    for(long long n = 0; n < sampleSize; n++)for(long long p = 0; p < 1; p++){
        ofs << "sample:" << n;
        for(long long l = 0; l < L; l++){
            ofs << '\t' << xOutput(n,p,l);
        }
        ofs << std::endl;
    }

    for(long long n = 0; n < sampleSize; n++)for(long long p = 1; p < posSizeEach[n]; p++){
        ofs << "sample:" << n;
        for(long long l = 0; l < L; l++){
            ofs << '\t' << xOutput(n,p,l);
        }
        ofs << std::endl;
    }
    ofs.close();
    // std::cerr << "finish making test FV file. " << std::endl;
    FeatureData<> fvData = parser.parseFeatureFile(testFVPath);
    // std::cerr << "finish read FV file. " << std::endl;

    // header file same?
    // L m1 ... mL //header
    // std::cerr << "fv header test " << std::endl;
    EXPECT_TRUE(fvData.I == sampleSize);
    EXPECT_TRUE(fvData.L == L);
    for(long long l = 0; l < L; l++){
        EXPECT_TRUE(ml[l] == fvData.Ml(l));
    }
    // main content sub information same?
    // n : position size
    // std::cerr << "fv prepare test " << std::endl;
    for(long long n = 0; n < sampleSize; n++){
        EXPECT_TRUE(fvData.Js(n) == posSizeEach[n]);
    }
    for(long long l = 0; l < L; l++){
        EXPECT_TRUE(fvData.maxJ == maxJ);
    }
    // main content same?
    // n : position size
    // std::cerr << "fv main test " << std::endl;
    for(long long n = 0; n < sampleSize; n++)for(long long p = 0; p < posSizeEach[n]; p++)for(long long l = 0; l < L; l++){
        for(long long v = 0; v < ml[l]; v++){
            EXPECT_TRUE( xOriginal(n,p,l,v) == fvData.X(n,p,l,v) );
        }
    }
}

TEST_F(PMSwithTest, math_apply) {
    using namespace pmswitch;
    double val = 1.2;
    FixedSizeMultiVector<double> vec(val, 2, 3, 4);
    FixedSizeMultiVector<double> applied = math::applied(vec, std::exp);
    for (int i = 0; i < 2; ++i)for(int j = 0; j < 3; j++)for (int k = 0; k < 4; ++k){
        EXPECT_EQ( std::exp(val), applied(i,j,k));
    }
    math::apply(vec, std::exp);
    for (int i = 0; i < 2; ++i)for(int j = 0; j < 3; j++)for (int k = 0; k < 4; ++k){
        EXPECT_EQ( std::exp(val), vec(i,j,k));
    }

    long long I = 2, maxJ = 3, maxK = 4;
    FixedSizeMultiVector<long long> Js(maxJ, I);
    Js(0) = 3; Js(1) = 2;

    FixedSizeMultiVector<long long> Ks(0, maxJ);
    Ks(0) = 1; Ks(1) = 3; Ks(2) = 2;

    FixedSizeMultiVector<bool> filter(true, I, maxJ, maxK);
    math::makeFilter(filter, (long long)I, Js, Ks);

    for(int i = 0; i < I; i++)for(int j = 0; j < maxJ; j++)for(int k = 0; k < maxK; k++){
        EXPECT_EQ( i < I && j < Js(i) && k < Ks(j), filter(i,j,k));
    }

    vec = FixedSizeMultiVector<double>(val, 2, 3, 4);
    applied = math::applied(vec, std::exp, filter);
    for(int i = 0; i < I; i++)for(int j = 0; j < maxJ; j++)for(int k = 0; k < maxK; k++){
        if(i < I && j < Js(i) && k < Ks(j)){
            EXPECT_EQ(std::exp(val), applied(i,j,k));
        }else{
            EXPECT_EQ(val, applied(i,j,k));
        }
    }

    math::apply(vec, std::exp, filter);
    for(int i = 0; i < I; i++)for(int j = 0; j < maxJ; j++)for(int k = 0; k < maxK; k++){
        if(i < I && j < Js(i) && k < Ks(j)){
            EXPECT_EQ(std::exp(val), vec(i,j,k));
        }else{
            EXPECT_EQ(val, vec(i,j,k));
        }
    }
}

TEST_F(PMSwithTest, math_norm) {
    using namespace pmswitch;
    double eps = 1e-9;
    FixedSizeMultiVector<double> vec(0.0, 2, 3, 4);
    std::vector<double> v = {0.1,0.2,0.3,0.4};
    double factor = 3.5;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 4; k++){
                vec(i,j,k) = v[k] * factor;
            }
        }
    }
    math::norm<>(vec, (long long)1);
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 4; k++){
                EXPECT_TRUE( std::abs( vec(i,j,k) - v[k] ) < eps );
            }
        }
    }
}

TEST_F(PMSwithTest, math_subMax) {
    using namespace pmswitch;
    double eps = 1e-9;
    FixedSizeMultiVector<double> vec(0.0, 2, 3, 4);
    std::vector<double> v = {0.1,0.2,0.3,0.4};
    double factor = 3.5;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 4; k++){
                vec(i,j,k) = v[k] * factor;
            }
        }
    }
    math::norm<>(vec, (long long)1);
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 4; k++){
                EXPECT_TRUE( std::abs( vec(i,j,k) - v[k] ) < eps );
            }
        }
    }

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 4; k++){
                vec(i,j,k) = 1.0;
            }
        }
    }
    math::norm<>(vec, (long long)0);
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 4; k++){
                EXPECT_TRUE( std::abs( vec(i,j,k) - 1.0/12.0 ) < eps );
            }
        }
    }
}

TEST_F(PMSwithTest, filter) {
    using namespace pmswitch;
    long long I = 3; long long maxJ = 4; long long maxK = 5;
    double initVal = 3.0;
    FixedSizeMultiVector<double> vec(initVal, I, maxJ, maxK);
    FixedSizeMultiVector<bool> filter(1, I, maxJ, maxK);
    double defaultVal = -1.0;

    FixedSizeMultiVector<long long> Js(0, I);
    Js(0) = 2; Js(1) = 3; Js(2) = 4;

    FixedSizeMultiVector<long long> Ks(0, maxJ);
    Ks(0) = 2; Ks(1) = 3; Ks(2) = 4; Ks(3) = 5;

    math::makeFilter(filter, (long long)I, Js, Ks);
    // std::cerr << "print filter" << std::endl;
    // for(int i = 0; i < filter.size(); i++){std::cerr << "i: " << i  << ", fitler[i]: " << filter[i] << std::endl; }

    for(long long i = 0; i < I; i++)for(long long j = 0; j < maxJ; j++)for(long long k = 0; k < maxK; k++){
        bool filtVal = (i < I) && (j < Js(i)) && (k < Ks(j));
        EXPECT_EQ(filtVal, filter(i,j,k));
        if(filtVal != filter(i,j,k)) std::cerr << "i: " << i << ", j: " << j << ", k: " << k << ", filter(i,j,k): " <<  filter(i,j,k) << ", ex: " << filtVal << std::endl;
    }
    math::filter(vec, filter, defaultVal);
    for(long long i = 0; i < I; i++)for(long long j = 0; j < maxJ; j++)for(long long k = 0; k < maxK; k++){
        // std::cerr << "i: " << i << ", j: " << j << ", k: " << k << std::endl;
        bool filtVal = (i < I) && (j < Js(i)) && (k < Ks(j));
        if(filtVal) EXPECT_EQ(3.0, vec(i,j,k));
        else        EXPECT_EQ(defaultVal, vec(i,j,k));
    }
    std::vector<int> a(10,10);
}


TEST_F(PMSwithTest, math_calDirExp) {
    using namespace pmswitch;
    pmswitch::FixedSizeMultiVector<double> D1Vec1(2.0, 20);
    pmswitch::FixedSizeMultiVector<double> D1Vec2(2.0, 5);

    pmswitch::FixedSizeMultiVector<double> D2Vec(2.0, 3, 4, 5);

    pmswitch::FixedSizeMultiVector<double> D1Vec1Dig = math::calDirExp(D1Vec1);
    pmswitch::FixedSizeMultiVector<double> D1Vec2Dig = math::calDirExp(D1Vec2);

    pmswitch::FixedSizeMultiVector<double> D2VecDig1 = math::calDirExp(D2Vec, (long long)0);
    pmswitch::FixedSizeMultiVector<double> D2VecDig2 = math::calDirExp(D2Vec, (long long)1);
{/*
    std::cerr << "##############################################" << std::endl;
    for(long long i = 0; i < 3; i++)for(int j = 0; j < 4; j++)for(int k = 0; k < 5; k++){
        std::cerr << "i,j,k: " << i << "," << j << "," << k << ", D2Vec(i,j,k), D2VecDig1(i,j,k): " << D2Vec(i,j,k) << ", " << D2VecDig1(i,j,k) << std::endl;
    }
    std::cerr << "##############################################" << std::endl;
    for(long long i = 0; i < 3; i++)for(int j = 0; j < 4; j++)for(int k = 0; k < 5; k++){
        std::cerr << "i,j,k: " << i << "," << j << "," << k << ", D2Vec(i,j,k), D2VecDig2(i,j,k): " << D2Vec(i,j,k) << ", " << D2VecDig2(i,j,k) << std::endl;
    }
    std::cerr << "##############################################" << std::endl;
    std::cerr << "D1Vec1Dig(0) : " << D1Vec1Dig(0) << std::endl;
    std::cerr << "D1Vec2Dig(0) : " << D1Vec2Dig(0) << std::endl;
*/}

    EXPECT_EQ(D1Vec1Dig(0), D2VecDig1(0,0));
    EXPECT_EQ(D1Vec2Dig(0), D2VecDig2(0,0));
}

TEST_F(PMSwithTest, math_calDirExpFilter) {
    using namespace pmswitch;

    // pmswitch::FixedSizeMultiVector<double> D1Vec1(2.0, 20);
    pmswitch::FixedSizeMultiVector<double> D1Vec(2.0, 4);
    pmswitch::FixedSizeMultiVector<double> D2Vec(2.0, 3, 4, 5);

    long long I = 3, maxJ = 4, maxK = 5;
    FixedSizeMultiVector<long long> Js(4, I);

    FixedSizeMultiVector<long long> Ks(0, maxJ);
    Ks(0) = 4; Ks(1) = 4; Ks(2) = 4; Ks(3) = 4;

    FixedSizeMultiVector<bool> filter(1, I, maxJ, maxK);
    {
        math::makeFilter(filter, (long long)I, Js, Ks);
    }

    pmswitch::FixedSizeMultiVector<double> D1VecDig1 = math::calDirExp(D1Vec);
    pmswitch::FixedSizeMultiVector<double> D2VecDig1 = math::calDirExp(D2Vec, (long long)1, filter);


    for(long long i = 0; i < I; i++)for(long long j = 0; j < maxJ; j++)for(long long k = 0; k < maxK; k++){
        // std::cerr << "i,j,k: " << i << "," << j << "," << k << ", D2VecDig1(i,j,k), D1VecDig1(k) :" << D2VecDig1(i,j,k) << ", " << D1VecDig1(k) << std::endl;
        if(k < 4)EXPECT_NEAR(D2VecDig1(i,j,k), D1VecDig1(k), 1e-9 );
    }
}

TEST_F(PMSwithTest, math_calELogDir) {
    using namespace pmswitch;

    pmswitch::FixedSizeMultiVector<double> D1Vec(0.01, 5);
    pmswitch::FixedSizeMultiVector<double> D2Vec(0.01, 3, 4, 5);

    pmswitch::FixedSizeMultiVector<double> D1VecParam(2.5, 5);
    pmswitch::FixedSizeMultiVector<double> D2VecParam(2.5, 3, 4, 5);

    double elogdirD1 = math::calELogDir(D1Vec, D1VecParam);
    double elogdirD2 = math::calELogDir(D2Vec, D2VecParam, (long long)1);

    EXPECT_NEAR(elogdirD1*12, elogdirD2, 1e-9);
}


TEST_F(PMSwithTest, math_calELogDirFilter) {
    using namespace pmswitch;

    pmswitch::FixedSizeMultiVector<double> D1Vec(0.01, 4);
    pmswitch::FixedSizeMultiVector<double> D2Vec(0.01, 3, 4, 5);

    pmswitch::FixedSizeMultiVector<double> D1VecParam(2.5, 4);
    pmswitch::FixedSizeMultiVector<double> D2VecParam(2.5, 3, 4, 5);

    long long I = 3, maxJ = 4, maxK = 5;
    FixedSizeMultiVector<long long> Js(4, I);

    FixedSizeMultiVector<long long> Ks(0, maxJ);
    Ks(0) = 4; Ks(1) = 4; Ks(2) = 4; Ks(3) = 4;

    FixedSizeMultiVector<bool> filter(1, I, maxJ, maxK);
    {
        math::makeFilter(filter, (long long)I, Js, Ks);
    }

    double elogdirD1 = math::calELogDir(D1Vec, D1VecParam);
    double elogdirD2 = math::calELogDir(D2Vec, D2VecParam, (long long)1, filter);

    EXPECT_NEAR(elogdirD1*12, elogdirD2, 1e-9);
}












