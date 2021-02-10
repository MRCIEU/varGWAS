//
// Created by Matt Lyon on 10/02/2021.
//

#ifndef JLST_CPP_SRC_PHENOTYPEFILEEXCEPTION_H_
#define JLST_CPP_SRC_PHENOTYPEFILEEXCEPTION_H_

#include <exception>
#include <string>

namespace jlst {
class PhenotypeFileException : public std::exception {
 public:
  explicit PhenotypeFileException(const std::string &msg);
  const char *what() const noexcept override;
 private:
  std::string msg;
};
}

#endif //JLST_CPP_SRC_PHENOTYPEFILEEXCEPTION_H_
