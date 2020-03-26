#include "csv.h"
#include "gtest/gtest.h"
#include <Eigen/Core>
#include <Eigen/SVD>

TEST(LinRegTest, slope_residual) {
    Eigen::MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;


    /*int x_f;
    int c1_f;
    double c2_f;
    double y_f;
    double d_f;

    // get data (see data/regression.R)
    io::CSVReader<5> in("data/regression.csv");
    in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y", "d");
    int t = 0;
    while (in.read_row(x_f, c1_f, c2_f, y_f, d_f)) {
        t++;
    }

    // linear regression
    Matrix615<double> tmpy(argv[1]); // read n * 1 matrix y
    Matrix615<double> tmpX(argv[2]); // read n * p marrix X
    int n = tmpX.numRows();
    int p = tmpX.numCols();
    MatrixXd y, X; // copy the matrices into Eigen::Matrix objects
    tmpy.copyTo(y);
    tmpX.copyTo(X);

    // assertions*/
}