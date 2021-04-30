//
// Created by Matt Lyon on 30/04/2021.
//

#ifndef JLST_CPP_SRC_PLINKPARSER_H_
#define JLST_CPP_SRC_PLINKPARSER_H_

namespace jlst {
class PlinkParser {
 public:
  explicit PlinkParser(std::string &file_path) : _file_path(file_path) {}
  std::vector<std::string> get_samples();
 private:
  std::string &_file_path;
};
}

#endif //JLST_CPP_SRC_PLINKPARSER_H_