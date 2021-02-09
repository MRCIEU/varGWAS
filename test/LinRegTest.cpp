#include "csv.h"
#include "gtest/gtest.h"
#include <Eigen/Core>
#include <Eigen/SVD>

/*
 * Test for performing linear regression model
 * */

TEST(LinRegTest, slope_residual) {
    std::cout << "Hi 23" << std::endl;
    const double intercept = 1.0;
    double x_f;
    double c1_f;
    double c2_f;
    double y_f;
    double d_f;
    int n = 50000;
    int p = 3;

    Eigen::MatrixXd X = Eigen::MatrixXd(n, p + 1);
    Eigen::VectorXd y = Eigen::VectorXd(n);

    /*// get data (see data/regression.R)
    std::cout << "Reading phenotype file..." << std::endl;
    io::CSVReader<5> in("data/regression.csv");
    in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y", "d");
    int t = 0;
    while (in.read_row(x_f, c1_f, c2_f, y_f, d_f)) {
        X(t, 0) = intercept;
        X(t, 1) = x_f;
        X(t, 2) = c1_f;
        X(t, 3) = c2_f;
        y(t, 0) = y_f;
        t++;
    }
    std::cout << "Done!" << std::endl;

    // linear regression using SVD
    // adapted from: https://genome.sph.umich.edu/w/images/2/2c/Biostat615-lecture14-presentation.pdf
    Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd betasSvd = svd.solve(y);

    // calculate VDˆ{-1}
    Eigen::MatrixXd ViD = svd.matrixV() * svd.singularValues().asDiagonal().inverse();
    double sigmaSvd = (y - X * betasSvd).squaredNorm() / (n - p); // compute \sigmaˆ2
    Eigen::MatrixXd varBetasSvd = sigmaSvd * ViD * ViD.transpose(); // Cov(\hat{beta})

    // assertions
    std::cout << "Int: " << betasSvd(0, 0) << std::endl;
    std::cout << "X: " << betasSvd(1, 0) << std::endl;
    std::cout << "C1: " << betasSvd(2, 0) << std::endl;
    std::cout << "C2: " << betasSvd(3, 0) << std::endl;

    ASSERT_NEAR(betasSvd(0, 0), 25, 0.1);
    ASSERT_NEAR(betasSvd(1, 0), 0.6, 0.1);
    ASSERT_NEAR(betasSvd(2, 0), 2, 0.1);
    ASSERT_NEAR(betasSvd(3, 0), 0.05, 0.002);*/
}