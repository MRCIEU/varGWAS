//
// Created by Matthew Lyon on 19/03/2020.
//

#include "gtest/gtest.h"
#include <glog/logging.h>

int main(int argc, char **argv) {
  // Initialize Google's logging library.
  google::InitGoogleLogging(argv[0]);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
