set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)

set(SOURCES ${SOURCES})

#${SQLITE_LIB_PATH}
add_executable(${BINARY} ${SOURCES})
target_link_libraries(
        ${BINARY}
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

add_library(${BINARY}_lib STATIC ${SOURCES})