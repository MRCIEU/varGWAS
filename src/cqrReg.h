//
// Created by Matt Lyon on 20/10/2021.
//

#include <armadillo>
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Dense>

#ifndef VARGWAS_SRC_CQRREG_H_
#define VARGWAS_SRC_CQRREG_H_

namespace cqrReg {
class cqrReg {
 public:
  static Eigen::VectorXd qrmm(Eigen::MatrixXd &X,
                              Eigen::VectorXd &Y,
                              Eigen::MatrixXd &ols_beta,
                              double toler,
                              int maxit,
                              double tau
  );
};
}

#endif //VARGWAS_SRC_CQRREG_H_
