//
// Created by Matt Lyon on 18/02/2021.
//

#ifndef JLST_CPP_SRC_SYNCHRONIZEDFILE_H_
#define JLST_CPP_SRC_SYNCHRONIZEDFILE_H_

#include <string>
#include <iostream>
#include <utility>
#include <fstream>
#include "Result.h"

namespace jlst {
class SynchronizedFile {
 public:
  explicit SynchronizedFile(const std::string &path) : path(path) {
    file.open(path);
    file << "CHR\tPOS\tRSID\tOA\tEA\tBETA\tSE\tP\tN\tEAF\n";
  }
  void write(const jlst::Result &result);
  void close();
 private:
  std::string path;
  std::ofstream file;
  std::mutex writer_mutex;
};
}

#endif //JLST_CPP_SRC_SYNCHRONIZEDFILE_H_
