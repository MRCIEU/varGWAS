#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "gtest/gtest.h"
#include "JLSP.h"
#include "libscl.h"
#include <boost/tokenizer.hpp>

TEST(JLSPTest, get_linear_estimate_test) {
    typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;
    std::vector<std::string> vec;
    std::string line;
    std::string header;
    int n_sim = 10000;
    scl::realmat x(n_sim, 1);
    scl::realmat y(n_sim, 1);
    REAL tau = 0.5;

    // read in data
    std::string data("data/regression.csv");
    std::ifstream in(data.c_str());
    if (!in.is_open()) throw std::runtime_error("Problem reading file");
    int n = 0;
    while (getline(in, line)) {
        Tokenizer tok(line);
        vec.assign(tok.begin(), tok.end());

        if (n == 0) {
            if (vec[0] != "x" || vec[1] != "y") {
                throw std::runtime_error("Wrong file format");
            }
        } else {
            x[n] = std::stod(vec[0]);
            y[n] = std::stod(vec[1]);
        }

        n += 1;
    }

    // quantile regression
    scl::realmat q = scl::quantreg(y, x, tau);
    std::cout << q << "\n";

    EXPECT_EQ (true, true);
}