set(
    SOURCES
        inc/advmath.hpp
        inc/complex.hpp
        inc/base.hpp
        inc/adfilter.hpp
        inc/genfilter.hpp
        src/complex.cpp
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

 add_subdirectory(test)
