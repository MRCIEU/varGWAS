# adapted from: https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html

cmake_minimum_required(VERSION 3.15)
project(jlst_cpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Threads)

# include
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/genfile/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/3rd_party/boost_1_55_0)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd-1.1.0/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/3rd_party/sqlite3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/db/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/fast-cpp-csv-parser)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen-3.3.9)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/cxxopts-2.2.0/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/ThreadPool)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# libraries
find_library(BGEN_LIB_PATH NAMES libbgen.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build)
find_library(BOOST_LIB_PATH NAMES libboost.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/boost_1_55_0)
find_library(ZSTD_LIB_PATH NAMES libzstd.so libzstd.dylib HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd-1.1.0/build/cmake/build/lib)
find_library(SQLITE_LIB_PATH NAMES libsqlite3.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/sqlite3)
find_library(DB_LIB_PATH NAMES libdb.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/db)
find_library(SPDLOG_LIB_PATH NAMES libglog.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/build)
find_library(LIBSCL_LIB_PATH NAMES libscl.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp)

#set(BGEN_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/libbgen.a)
#set(BOOST_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/boost_1_55_0/libboost.a)
#set(ZSTD_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd-1.1.0/build/cmake/build/lib/libzstd.so)
#set(SQLITE_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/sqlite3/libsqlite3.a)
#set(DB_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/db/libdb.a)
#set(SPDLOG_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/build/libglog.a)
#set(LIBSCL_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp/libscl.a)

# copy dynamic libs to lib folder
#file(COPY ${ZSTD_LIB_PATH} DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#get_filename_component(zstdlib ${ZSTD_LIB_PATH} NAME)

# src
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

# build tests if in debug mode
if (${CMAKE_BUILD_TYPE} STREQUAL Debug)
    message("Building tests")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0)
endif ()