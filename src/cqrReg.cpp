//
// Created by Matt Lyon on 20/10/2021.
//

#include <armadillo>
#include <Eigen/Core>
#include "cqrReg.h"

namespace cqrReg {
Eigen::VectorXd cqrReg::qrmm(Eigen::MatrixXd &X,
                             Eigen::VectorXd &Y,
                             Eigen::MatrixXd &ols_beta,
                             double toler,
                             int maxit,
                             double tau) {
  int n = X.rows();
  int p = X.cols();
  arma::mat x = arma::mat(X.data(), X.rows(), X.cols(), false, false);
  arma::vec y = arma::vec(Y.data(), Y.rows(), false, false);
  arma::vec beta = arma::vec(ols_beta.data(), ols_beta.rows(), false, false);
  arma::mat product, xt;
  arma::vec W, newX, z, signw, v, r;
  arma::vec delta;
  arma::uvec order, index;

  double error = 10000, epsilon = 0.9999;
  int iteration = 1;
  product.ones(p, n);
  xt = x.t();

  while (iteration <= maxit && error > toler) {
    r = y - x * beta;
    v = 1 - 2 * tau - r / (arma::abs(r) + epsilon);

    W = 1 / (epsilon + arma::abs(r));

    for (int i = 0; i < n; i++) { product.col(i) = xt.col(i) * W(i); }

    delta = arma::solve(product * x, xt * v);
    beta = beta - delta;

    error = sum(abs(delta));
    iteration++;
  }

  Eigen::VectorXd res = Eigen::Map<Eigen::VectorXd>(beta.memptr(), beta.n_rows, beta.n_cols);

  return(res);
}
}