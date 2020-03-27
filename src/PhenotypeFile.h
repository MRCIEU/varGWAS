//
// Created by Matthew Lyon on 27/03/2020.
//

#ifndef JLST_CPP_PHENOTYPEFILE_H
#define JLST_CPP_PHENOTYPEFILE_H

#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include <cassert>
#include <stdexcept>


class PhenotypeFile {
public:
    PhenotypeFile(std::string &filename);
    void readfile();


private:
    std::string filename;
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
};


#endif //JLST_CPP_PHENOTYPEFILE_H
