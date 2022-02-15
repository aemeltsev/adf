set(SOURCES
    inc/advmath.hpp
    inc/complex.hpp
    inc/base.hpp
    inc/adfilter.hpp
    inc/genfilter.hpp
    inc/anfilter.hpp
    inc/iirfilter.hpp
    inc/fresp.hpp
    inc/polynomial.hpp
    inc/fft.hpp
    src/complex.cpp
    src/polynomial.cpp
    )

#Declare library
add_library(adf STATIC
           ${SOURCES})

# Specify here the include directories exported
# by this library
target_include_directories(adf PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    #adf/inc
    )
#target_include_directories(adf PUBLIC
    #${CMAKE_CURRENT_SOURCE_DIR}
    #adf/src
    #)
target_compile_features(adf PUBLIC cxx_std_11)

install(
    TARGETS adf
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
 )



# Add configs
#configure_file(adfilter.pc.in adfilter.pc @ONLY)

# Compiler settings
#target_compile_features(${BIN} PUBLIC cxx_std_11)

# Install and export the library
#install(TARGETS adf
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#install(FILES ${CMAKE_BINARY_DIR}/adfilter.pc
#    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

##############################################
## Add test
enable_testing()
add_subdirectory(test)
