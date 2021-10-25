set(adf_sources
    inc/advmath.hpp
    inc/complex.hpp
    inc/base.hpp
    inc/adfilter.hpp
    inc/genfilter.hpp
    src/adfilter.cpp
    src/complex.cpp
    src/genfilter.cpp
    )

#Declare library
add_library(adf STATIC
           ${adf_sources})

# Specify here the include directories exported
# by this library
target_include_directories(adf PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    #adf/inc
    )
target_include_directories(adf PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    #adf/src
    )
