set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)

set(SOURCES ${SOURCES})

add_executable(${BINARY}_run ${SOURCES})
target_link_libraries(
        ${BINARY}_run
        z
        ${BGEN_LIB_PATH}
        ${BOOST_LIB_PATH}
        ${ZSTD_LIB_PATH}
        ${SQLITE_LIB_PATH}
        ${DB_LIB_PATH}
        ${LIBSCL_LIB_PATH}
        ${GLOG_LIB_PATH}
)

add_library(${BINARY}_lib STATIC ${SOURCES})