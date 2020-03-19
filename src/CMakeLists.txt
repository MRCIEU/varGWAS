set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)

set(SOURCES ${SOURCES})

add_executable(${BINARY}_run ${SOURCES} ${LIB_SRC_FILES})
target_link_libraries(${BINARY}_run z)

add_library(${BINARY}_lib STATIC ${SOURCES} ${LIB_SRC_FILES})