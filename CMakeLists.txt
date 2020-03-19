cmake_minimum_required(VERSION 3.15)
project(jlst_cpp)

set(CMAKE_CXX_STANDARD 14)

set(BOOST_SOURCE_DIR "includes/bgen/3rd_party/boost_1_55_0")
set(ZSTD_SOURCE_DIR "includes/bgen/3rd_party/zstd-1.1.0")
set(SQLITE_SOURCE_DIR "includes/bgen/3rd_party/sqlite3")
set(DB_SOURCE_DIR "includes/bgen/db")
set(BGEN_SOURCE_DIR "includes/bgen/src")
set(GENFILE_SOURCE_DIR "includes/bgen/genfile")
set(LIBSCL_SOURCE_DIR "includes/libscl_float/src")

file(GLOB BOOST_SRC_FILES
        ${BOOST_SOURCE_DIR}/libs/system/src/*.cpp
        ${BOOST_SOURCE_DIR}/libs/thread/src/*.cpp
        ${BOOST_SOURCE_DIR}/libs/thread/src/pthread/once_atomic.cpp
        ${BOOST_SOURCE_DIR}/libs/thread/src/pthread/thread.cpp
        ${BOOST_SOURCE_DIR}/libs/filesystem/src/*.cpp
        ${BOOST_SOURCE_DIR}/libs/date_time/src/posix_time/*.cpp
        ${BOOST_SOURCE_DIR}/libs/timer/src/*.cpp
        ${BOOST_SOURCE_DIR}/libs/chrono/src/*.cpp
        )

file(GLOB DB_SRC_FILES ${DB_SOURCE_DIR}/src/*.cpp)
file(GLOB ZSTD_SRC_FILES
        ${ZSTD_SOURCE_DIR}/lib/common/*.c
        ${ZSTD_SOURCE_DIR}/lib/decompress/*.c
        ${ZSTD_SOURCE_DIR}/lib/compress/*.c
        )
file(GLOB BGEN_SRC_FILES ${BGEN_SOURCE_DIR}/*.cpp)

include_directories(${GENFILE_SOURCE_DIR}/include)
include_directories(${DB_SOURCE_DIR}/include)
include_directories(${BOOST_SOURCE_DIR})
include_directories(${ZSTD_SOURCE_DIR})
include_directories(${ZSTD_SOURCE_DIR}/lib)
include_directories(${ZSTD_SOURCE_DIR}/lib/common)
include_directories(${ZSTD_SOURCE_DIR}/lib/compress)
include_directories(${ZSTD_SOURCE_DIR}/lib/decompress)
include_directories(${SQLITE_SOURCE_DIR})
include_directories(${LIBSCL_SOURCE_DIR})

add_definitions(-DSQLITE_ENABLE_COLUMN_METADATA -DSQLITE_ENABLE_STAT4 -DSQLITE_MAX_EXPR_DEPTH=10000 -DSQLITE_USE_URI=1)

set(SOURCE
        ${BOOST_SRC_FILES}
        ${DB_SRC_FILES}
        ${ZSTD_SRC_FILES}
        ${SQLITE_SOURCE_DIR}/sqlite3/sqlite3.c
        ${BGEN_SRC_FILES}
        ${LIBSCL_SOURCE_DIR}/quantreg.cpp
        ${LIBSCL_SOURCE_DIR}/unsk.cpp
        ${LIBSCL_SOURCE_DIR}/ran.cpp
        ${LIBSCL_SOURCE_DIR}/intvec.cpp
        ${LIBSCL_SOURCE_DIR}/realmat.cpp
        ${LIBSCL_SOURCE_DIR}/fmt.cpp
        ${LIBSCL_SOURCE_DIR}/dcnd.cpp
        ${LIBSCL_SOURCE_DIR}/dginv.cpp
        ${LIBSCL_SOURCE_DIR}/dgmpnt.cpp
        ${LIBSCL_SOURCE_DIR}/dsolve.cpp
        ${LIBSCL_SOURCE_DIR}/dsweep.cpp
        ${LIBSCL_SOURCE_DIR}/dsvd.cpp
        ${LIBSCL_SOURCE_DIR}/sclerror.cpp
        ${LIBSCL_SOURCE_DIR}/starbox.cpp
        main.cpp
        DosageSetter.cpp
        DosageSetter.h
        ProbSetter.cpp
        ProbSetter.h
        BgenParser.cpp
        BgenParser.h
        JLSP.cpp
        JLSP.h
        )

add_executable(jlst_cpp ${SOURCE})
add_library(jlst_lib STATIC ${SOURCE})
target_link_libraries(jlst_cpp z)

#add_subdirectory(tests)