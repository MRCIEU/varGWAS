set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)

set(SOURCES ${SOURCES})

add_executable(${BINARY} ${SOURCES})

target_link_libraries(
        ${BINARY}
        ${LIBS}
        ${LIBS}
)

add_library(${BINARY}_lib STATIC ${SOURCES})