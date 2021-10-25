cmake_minimum_required(VERSION 3.5)

# Project name, version & main languages
set(BIN adf)
project(${BIN} VERSION 0.0.1 LANGUAGES CXX)
include(GNUInstallDirs)

add_subdirectory(adf)

find_package(GTest)
if(GTEST_FOUND)
    enable_testing()
    add_subdirectory(test)
    #add_test()
    #add_dependencies()
endif()

# Add configs
#configure_file(adfilter.pc.in adfilter.pc @ONLY)

# Compiler settings
target_compile_features(${BIN} PUBLIC cxx_std_11)

# Install and export the library
#install(TARGETS adf
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#install(FILES ${CMAKE_BINARY_DIR}/adfilter.pc
#    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
