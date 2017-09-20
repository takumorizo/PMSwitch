#ifndef ReadTest_Parameters_h
#define ReadTest_Parameters_h

#include <cstdlib>
#include <string>
#include <vector>
#include "cmdline.h"
#include <iostream>

namespace pmswitch{
template<typename Int = long long, typename Real = double>
class Parameters {
public:
    Parameters(int argc, const char *argv[]);
    void getFromCommandLineArguments(int argc, const char *argv[]);

    typedef enum { VB,ONLINEVB,MCMC } ALGO;
    ALGO method;
    std::string pathFV;
    std::string pathDB;
    std::string pathOutputDir;
    bool isDBUpdate;
    Int T;
    Real alpha;
    Real beta0, beta1;
    Real gamma;
    Real eta;
    Real s0, s1;
    Real u0, u1;

    std::string alphaPath;
    std::string betaPath;
    std::string gammaPath;
    std::string etaPath;
    std::string sPath;
    std::string uPath;
};

template<typename Int, typename Real>
Parameters<Int, Real>::Parameters(int argc, const char *argv[]){

}

template<typename Int, typename Real>
void Parameters<Int, Real>::getFromCommandLineArguments(int argc, const char *argv[]){
    cmdline::parser a;

    a.add<std::string>("input",     'i',  ".fv file for mutation feature",                         true,  "");
    a.add<std::string>("outDir",    'o',  "result directory for mutation feature",                 true,  "");
    a.add<std::string>("db" ,       'd',  "database for mutation feature",                         true,  "");
    a.add<std::string>("algorithm", '\0', "algorithm type used in inference",                      false, "");
    a.add("isDBUpdate",'\0', "if true, update existing database frequencies");
    a.add<Int>("truncation",        'T',  "truncation number used in dirichlet process inference", false, 10);

    a.add<Real>("alpha", '\0', "alpha value",    false, 1e-5);
    a.add<Real>("beta0", '\0', "beta  value @ 0", false, 1);
    a.add<Real>("beta1", '\0', "beta  value @ 1", false, 100);
    a.add<Real>("gamma", '\0', "gamma value",    false, 1);
    a.add<Real>("eta",   '\0', "gamma value",    false, 1);
    a.add<Real>("s0",    '\0', "s0    value",    false, 1.1);
    a.add<Real>("s1",    '\0', "s1    value",    false, 0.001);
    a.add<Real>("u0",    '\0', "u0    value",    false, 1.1);
    a.add<Real>("u1",    '\0', "u1    value",    false, 0.001);


    a.add<std::string>("alphaPath", '\0', "alpha vector file", false, "");
    a.add<std::string>("betaPath",  '\0', "beta vector file",  false, "");
    a.add<std::string>("gammaPath", '\0', "gamma vector file", false, "");
    a.add<std::string>("etaPath",   '\0', "gamma vector file", false, "");
    a.add<std::string>("sPath",     '\0', "s vector file",     false, "");
    a.add<std::string>("uPath",     '\0', "u vector file",     false, "");

    a.parse_check(argc, argv);

    std::string methodStr = a.get<std::string>("algorithm");
    if(methodStr == "VB"){
        method = Parameters<Int, Real>::VB;
    }else if(methodStr == "ONLINEVB"){
        method = Parameters<Int, Real>::ONLINEVB;
    }else if(methodStr == "MCMC"){
        method = Parameters<Int, Real>::MCMC;
    }else{
        std::cerr << "algorithm you specified is not in VB, ONLINEVB, MCMC. Specify only one of these." << std::endl; ;
        exit(1);
    }

    pathFV        = a.get<std::string>("input");
    pathDB        = a.get<std::string>("db");
    pathOutputDir = a.get<std::string>("outDir");
    isDBUpdate    = a.exist("isDBUpdate") ? true : false;
    T             = a.get<Int>("truncation");

    alpha = a.get<Real>("alpha");
    beta0 = a.get<Real>("beta0");
    beta1 = a.get<Real>("beta1");
    gamma = a.get<Real>("gamma");
    eta   = a.get<Real>("eta");
    s0    = a.get<Real>("s0");
    s1    = a.get<Real>("s1");
    u0    = a.get<Real>("u0");
    u1    = a.get<Real>("u1");


    alphaPath = a.get<std::string>("alphaPath");
    betaPath  = a.get<std::string>("betaPath");
    gammaPath = a.get<std::string>("gammaPath");
    etaPath   = a.get<std::string>("etaPath");
    sPath     = a.get<std::string>("sPath");
    uPath     = a.get<std::string>("uPath");
}

}

#endif
