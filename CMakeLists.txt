cmake_minimum_required(VERSION 3.5)

############################################################
# Project Name Version & main Languages
############################################################

set(BIN adf)
project(${BIN} VERSION 0.0.1 LANGUAGES CXX)
include(GNUInstallDirs)

############################################################
# Create a library
############################################################

add_library(adf STATIC
    adf/inc/advmath.hpp
    adf/inc/complex.hpp
    adf/inc/base.hpp
    adf/inc/adfilter.hpp
    adf/inc/genfilter.hpp
    adf/src/adfilter.cpp
    adf/src/complex.cpp
    adf/src/genfilter.cpp
)
############################################################
# Add directories
############################################################
target_include_directories(adf PRIVATE adf/inc)
target_include_directories(adf PRIVATE adf/src)

############################################################
# Version
############################################################
set_target_properties(adf PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(adf PROPERTIES SOVERSION 1)

############################################################
# Public heders
############################################################
set(MYLIB_HEADERS adf/inc/adfilter.hpp)
set_target_properties(adf PROPERTIES PUBLIC_HEADER ${MYLIB_HEADERS})

############################################################
# Public heders
############################################################
configure_file(adfilter.pc.in adfilter.pc @ONLY)

target_compile_features(adf PUBLIC cxx_std_17)

# install and export the library

install(TARGETS adf
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/adfilter.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
