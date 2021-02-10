//
// Created by Matt Lyon on 10/02/2021.
//
#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include "PhenotypeFile.h"
#include "BgenParser.h"

#ifndef JLST_CPP_SRC_ASSOC_H_
#define JLST_CPP_SRC_ASSOC_H_

namespace jlst {
class Model {
 public:
  static void run(jlst::PhenotypeFile &phenotype_file, genfile::bgen::BgenParser &bgen_parser);
 private:
};
}

#endif //JLST_CPP_SRC_ASSOC_H_
