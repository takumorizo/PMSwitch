#ifndef ReadTest_Parameters_h
#define ReadTest_Parameters_h

#include <cstdlib>
#include <string>
#include <vector>
#include "cmdline.h"
#include <iostream>

template<typename Int = int, typename Real = double>
class Parameters {
public:
    Parameters(int argc, const char *argv[]);
    void getFromCommandLineArguments(int argc, const char *argv[]);

    typedef enum { VB,ONLINEVB,MCMC } ALGO;
    ALGO method;
    std::string pathFV;
    std::string pathDB;
    std::string pathOutput;
    bool isDBUpdate;
    Int T;
    Real alpha;
};

template<typename Int, typename Real>
Parameters<Int, Real>::Parameters(int argc, const char *argv[]){

}

template<typename Int, typename Real>
void Parameters<Int, Real>::getFromCommandLineArguments(int argc, const char *argv[]){
    cmdline::parser a;

    a.add<std::string>("input",     'i',  ".fv file for mutation feature",    true, "");
    a.add<std::string>("out",       'o',  "result for mutation feature",      true, "");
    a.add<std::string>("db" ,       'd',  "database for mutation feature",    true, "");
    a.add<std::string>("algorithm", '\0', "algorithm type used in inference", false, "");
    a.add("isDBUpdate", '\0', "if true, update existing database frequencies");
    a.add<Int>("Truncation", 'T', "Truncation number used in dirichlet process inference", false, 10);

    a.parse_check(argc, argv);

    pathFV = a.get<std::string>("input");
    pathDB = a.get<std::string>("db");
    pathOutput = a.get<std::string>("out");
    isDBUpdate = a.exist("isDBUpdate") ? true : false;
    std::string methodStr = a.get<std::string>("algorithm");

    if(methodStr == "VB"){
        method = Parameters<Int>::VB;
    }else if(methodStr == "ONLINEVB"){
        method = Parameters<Int>::ONLINEVB;
    }else if(methodStr == "MCMC"){
        method = Parameters<Int>::MCMC;
    }else{
        std::cerr << "algorithm you specified is not in VB, ONLINEVB, MCMC. Specify only one of these." << std::endl; ;
        exit(1);
    }
}


#endif
