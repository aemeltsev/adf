cmake_minimum_required(VERSION 3.5)

# Project Name Version & main Languages

project(adf VERSION 0.0.1 LANGUAGES CXX)
include(GNUInstallDirs)
include_directories(adf/inc adf/src)

# define the library

add_library(adf STATIC
    adf/inc/advmath.hpp
    adf/inc/complex.hpp
    adf/inc/errornum.hpp
    adf/inc/adfilter.hpp
    adf/inc/genfilter.hpp
    adf/src/adfilter.cpp
    adf/src/advmath.cpp
    adf/src/complex.cpp
    adf/src/errornum.cpp
    adf/src/genfilter.cpp
)

set_target_properties(adf PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(adf PROPERTIES SOVERSION 1)

set_target_properties(adf PROPERTIES PUBLIC_HEADER adf/inc/adfilter.hpp)

configure_file(adfilter.pc.in adfilter.pc @ONLY)

target_include_directories(adf PRIVATE include)
target_include_directories(adf PRIVATE src)

# install and export the library

install(TARGETS adf
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/adfilter.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
