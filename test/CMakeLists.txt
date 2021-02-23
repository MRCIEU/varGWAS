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
        ${BGEN_LIB_PATH}
        ${BOOST_LIB_PATH}
        ${ZSTD_LIB_PATH}
        ${SQLITE_LIB_PATH}
        ${DB_LIB_PATH}
        ${LIBSCL_LIB_PATH}
        ${GLOG_LIB_PATH}
        ${CMAKE_THREAD_LIBS_INIT}
        ${SQLITE3_LIBRARIES}
)

add_subdirectory(data)