set(TESTS
        adftest.cpp
        main.cpp)
#include_directories(${CMAKE_SOURCE_DIR}/adf/inc ${CMAKE_SOURCE_DIR}/adf/src)
#include_directories(${GTEST_INCLUDE_DIRS})

add_executable(adf_test ${TESTS})
target_link_libraries(adf_test adf ${GTEST_BOTH_LIBRARIES} pthread)
#add_test(outtest outtest)
#add_dependencies(adf_tests outtest)
