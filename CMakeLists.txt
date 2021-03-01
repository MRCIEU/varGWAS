# adapted from: https://raymii.org/s/tutorials/Cpp_project_setup_with_cmake_and_unit_tests.html

cmake_minimum_required(VERSION 3.15)
project(jlst_cpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Threads)

# include
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/genfile/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/boost_1_75_0)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd-1.1.0/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/3rd_party/sqlite3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/db/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/fast-cpp-csv-parser)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen-master)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/cxxopts-2.2.0/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/ThreadPool)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# libraries
find_library(BGEN_LIB_PATH NAMES libbgen.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build)
# TODO upgrade
find_library(ZSTD_LIB_PATH NAMES libzstd.so libzstd.dylib HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd-1.1.0/build/cmake/build/lib)
# TODO upgrade
find_library(SQLITE_LIB_PATH NAMES libsqlite3.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/3rd_party/sqlite3)
find_library(DB_LIB_PATH NAMES libdb.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/bgen/build/db)
find_library(SPDLOG_LIB_PATH NAMES libglog.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/build)
find_library(LIBSCL_LIB_PATH NAMES libscl.a HINTS ${CMAKE_CURRENT_SOURCE_DIR}/lib/libscl/gpp)

# src
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

# build tests & configure for performance profiling if in debug mode
if (${CMAKE_BUILD_TYPE} STREQUAL Debug)
    message("Building in DEBUG mode")

    # adding test libs
    message("Building tests")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest-release-1.10.0)

    # To make sure none of the frames are missing due to inlining
    message("Setting flags for performance profiling")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

    # compilers can use the frame pointer register as a general-purpose register for optimization purposes
    # which may lead to broken stack traces
    if (APPLE)
        if (${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
            message("Adding flags for Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer")
            set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer")
        elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
            message("Adding flags for GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
            set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer")
        endif ()
    endif ()

endif ()