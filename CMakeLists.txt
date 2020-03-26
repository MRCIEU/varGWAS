# adapted from: https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html

cmake_minimum_required(VERSION 3.15)
project(jlst_cpp)

set(CMAKE_CXX_STANDARD 14)

# include
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/genfile/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/3rd_party/boost_1_55_0)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/3rd_party/zstd-1.1.0/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/3rd_party/sqlite3/sqlite3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/db/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp)
include_directories(src)

# libraries
set(BGEN_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/build/libbgen.a)
set(BOOST_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/build/3rd_party/boost_1_55_0/libboost.a)
set(ZSTD_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/build/3rd_party/zstd-1.1.0/libzstd.a)
set(SQLITE_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/build/3rd_party/sqlite3/libsqlite3.a)
set(DB_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/gavinband-bgen-44fcabbc5c38/build/db/libdb.a)
set(LIBSCL_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp/libscl.a)

# src and test cmakelist
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(lib/googletest-release-1.10.0)