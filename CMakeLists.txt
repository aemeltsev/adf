# Project name, version & main languages
##################################################
cmake_minimum_required(VERSION 3.5)
set                   (BIN adf)
project               (${BIN} VERSION 0.0.1 LANGUAGES CXX)
set                   (CMAKE_CXX_STANDARD 17)

if(MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()

# Targets
##################################################
add_subdirectory(adf)

# Options
##################################################
option(BUILD_TESTS "Build tests." ON)

if(BUILD_TESTS)
    enable_testing()
    configure_file(CMakeLists.txt.in
        googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}googletest-build)

    add_subdirectory(test)
endif()

# Installation TODO corrected
##################################################

install(
    TARGETS ${BIN}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
 )
