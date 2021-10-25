
include_directories(${CMAKE_SOURCE_DIR}/adf/inc)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(adf_test main.cpp adftest.hpp)
target_link_libraries(adf_test ${GTEST_BOTH_LIBRARIES} pthread adf)
#add_test(outtest outtest)
#add_dependencies(adf_tests outtest)
