// Taken from https://raw.githubusercontent.com/cran/cqrReg/master/src/QRMM.cpp
// Pietrosanu, M., Gao, J., Kong, L., Jiang, B., and Niu, D. (2020). Advanced algorithms for penalized quantile and composite quantile regression. Comput. Stat. 2020 361 36, 333–346.

#include <Eigen/Core>
#include <Eigen/QR>
#include <iostream>
#include <sstream>

#ifndef VARGWAS_SRC_CQRREG_H_
#define VARGWAS_SRC_CQRREG_H_

namespace CqrReg {
class QRMM {
 public:
  static Eigen::VectorXd fit(const Eigen::MatrixXd &X,
                             const Eigen::VectorXd &Y,
                             const Eigen::VectorXd &init,
                             double toler,
                             int maxit,
                             double tau
  );
};
}

#endif //VARGWAS_SRC_CQRREG_H_
