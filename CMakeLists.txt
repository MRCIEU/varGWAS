# adapted from: https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html

cmake_minimum_required(VERSION 3.15)
project(jlst_cpp)

set(CMAKE_CXX_STANDARD 14)

# include
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/genfile/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/3rd_party/boost_1_55_0)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/3rd_party/zstd-1.1.0/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/3rd_party/sqlite3/sqlite3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/db/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/fast-cpp-csv-parser)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen-3.3.7)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/cxxopts-2.2.0/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/ThreadPool)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# libraries
set(BGEN_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/libbgen.a)
set(BOOST_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/boost_1_55_0/libboost.a)
set(ZSTD_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/zstd-1.1.0/libzstd.a)
set(SQLITE_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/sqlite3/libsqlite3.a)
set(DB_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/db/libdb.a)
set(SPDLOG_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/build/libglog.a)
set(LIBSCL_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp/libscl.a)

# src and test cmakelist
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0)