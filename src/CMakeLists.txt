set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)

set(SOURCES ${SOURCES})

add_executable(${BINARY} ${SOURCES})

target_link_libraries(
        ${BINARY}
        ${BGEN_LIB_PATH}
        ZLIB::ZLIB
        OpenMP::OpenMP_CXX
        ${ZSTD_LIB_PATH}
        ${SQLITE_LIB_PATH}
        ${DB_LIB_PATH}
        ${SPDLOG_LIB_PATH}
        ${ARMADILLO_LIBRARIES}
)

add_library(${BINARY}_lib STATIC ${SOURCES})