//
// Created by Matt Lyon on 10/02/2021.
//

#include "PhenotypeFileException.h"

// TODO convert to parent class and subclass
namespace jlst {
PhenotypeFileException::PhenotypeFileException(const std::string &msg) {
  this->msg = msg;
};
const char *PhenotypeFileException::what() const noexcept {
  return msg.c_str();
};
}