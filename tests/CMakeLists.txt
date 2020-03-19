project(tests)
set(GTEST_SOURCE_DIR "lib/googletest-release-1.10.0")
add_subdirectory(${GTEST_SOURCE_DIR})
include_directories(${GTEST_SOURCE_DIR}/include ${GTEST_SOURCE_DIR})

# adding the Google_Tests_run target
add_executable(Google_Tests_run JLSPTest.cpp)

# linking Google_Tests_run with JLSP which will be tested
target_link_libraries(Google_Tests_run jlst_lib)

target_link_libraries(Google_Tests_run gtest gtest_main)