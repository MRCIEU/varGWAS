// Adapted from https://raw.githubusercontent.com/cran/cqrReg/master/src/QRMM.cpp
// Pietrosanu, M., Gao, J., Kong, L., Jiang, B., and Niu, D. (2020). Advanced algorithms for penalized quantile and composite quantile regression. Comput. Stat. 2020 361 36, 333–346.

#include <Eigen/Core>
#include <Eigen/QR>
#include "QRMM.h"

namespace CqrReg {
Eigen::VectorXd QRMM::fit(
    const Eigen::MatrixXd &X,
    const Eigen::VectorXd &y,
    const Eigen::VectorXd &init,
    double toler,
    int maxit,
    double tau) {
  Eigen::VectorXd beta = init;
  Eigen::MatrixXd xt = X.transpose();

  double error = 10000, epsilon = 0.9999;
  int iteration = 1;

  while (iteration <= maxit && error > toler) {
    Eigen::VectorXd fitted = X * beta;
    Eigen::VectorXd r = y - fitted;
    Eigen::VectorXd v = 1 - 2 * tau - r.array() / (r.array().abs() + epsilon);

    Eigen::VectorXd W = 1 / (epsilon + r.array().abs());
    Eigen::MatrixXd product = xt.array().rowwise() * W.transpose().array();
    Eigen::MatrixXd pX = product * X;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(pX);
    if (qr.rank() < pX.cols()) {
      throw std::runtime_error("rank-deficient matrix");
    }
    Eigen::VectorXd delta = qr.solve(xt * v);
    beta = beta - delta;

    error = delta.cwiseAbs().sum();
    iteration++;
  }

  return (beta);
}
}