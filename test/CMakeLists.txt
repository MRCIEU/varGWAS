set(BINARY ${CMAKE_PROJECT_NAME}_test)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES})
add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(
        ${BINARY}
        PUBLIC
        ${CMAKE_PROJECT_NAME}_lib
        gtest
        z
        ${ZSTD_LIB_PATH}
        ${SQLITE_LIB_PATH}
        ${DB_LIB_PATH}
        ${BGEN_LIB_PATH}
        ${LIBSCL_LIB_PATH}
        ${GLOG_LIB_PATH}
        ${LIBPLINKIO_LIB_PATH}
        ${LIBCSV_LIB_PATH}
        OpenMP::OpenMP_CXX
        ${CMAKE_THREAD_LIBS_INIT}
)

add_subdirectory(data)