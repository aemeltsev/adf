# Testing
##################################################
set(TESTS
    adftest.hpp
    polytest.hpp
    comptest.hpp
    main.cpp)

add_executable(adf_test ${TESTS})
target_compile_features(adf_test PRIVATE cxx_auto_type)
target_link_libraries(adf_test adf GTest::gtest GTest::gtest_main pthread)
add_test(NAME adf_test COMMAND adf_test)
