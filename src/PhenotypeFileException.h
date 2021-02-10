//
// Created by Matt Lyon on 10/02/2021.
//

#ifndef JLST_CPP_SRC_PHENOTYPEFILEEXCEPTION_H_
#define JLST_CPP_SRC_PHENOTYPEFILEEXCEPTION_H_

namespace jlst {
class PhenotypeFileException : public exception {

 public:
  PhenotypeFileException(const std::string &ex);
  virtual const char *what() const throw();
 private:
  std::string ex;

};
}

#endif //JLST_CPP_SRC_PHENOTYPEFILEEXCEPTION_H_
