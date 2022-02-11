find_package(GTest QUIET)

if(NOT GTEST_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/cmake/GoogleTest-CMakeLists.txt.in
        ${CMAKE_BINARY_DIR}/googletest-download/CMakeLists.txt)
    execute_process(COMMAND &{CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)
endif()

set(TESTS
        adftest.cpp
        main.cpp)

add_executable(adf_test ${TESTS})
target_compile_features(adf_test PRIVATE cxx_auto_type)
target_link_libraries(adf_test adf gtest_main pthread)
add_test(NAME adf_test COMMAND adf_test)
