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

# Options
##################################################
option(BUILD_TESTS "Build tests." ON)

# Targets
##################################################
include(GNUInstallDirs)

add_subdirectory(adf)

# Installation TODO corrected
##################################################

install(
    TARGETS ${BIN}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
 )
