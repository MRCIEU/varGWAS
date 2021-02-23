set(BINARY ${CMAKE_PROJECT_NAME}_test)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES})
add_test(NAME ${BINARY} COMMAND ${BINARY})

#${SQLITE_LIB_PATH}
target_link_libraries(
        ${BINARY}
        PUBLIC
        ${CMAKE_PROJECT_NAME}_lib
        gtest
        z
        ${BOOST_LIB_PATH}
        ${ZSTD_LIB_PATH}
        ${SQLITE3_LIBRARIES}
        ${DB_LIB_PATH}
        ${BGEN_LIB_PATH}
        ${LIBSCL_LIB_PATH}
        ${GLOG_LIB_PATH}
        ${CMAKE_THREAD_LIBS_INIT}
)

add_subdirectory(data)