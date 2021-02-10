//
// Created by Matt Lyon on 10/02/2021.
//

#include "PhenotypeFileException.h"

namespace jlst {
PhenotypeFileException::PhenotypeFileException(const std::string &ex) : ex{

}
virtual const char *PhenotypeFileException::what() const throw() {
  return ex.c_str();
}
}
