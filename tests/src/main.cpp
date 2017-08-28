#include "gtest/gtest.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/foreach.hpp>


#include <iostream>
#include <list>
#include <string>

#include "sim_test.hpp"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        std::cerr << "please specify the path of data folder" << std::endl;
        exit(1);
    }
    data_path = std::string(argv[1]);
    return RUN_ALL_TESTS();
    // return RUN_ALL_TESTS();
    std::cout << "aaa" << std::endl;
}

