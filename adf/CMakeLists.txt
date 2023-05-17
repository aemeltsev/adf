# Sources
##################################################
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
    inc/firfilter.hpp
    src/complex.cpp
    src/polynomial.cpp
    )

# Targets
##################################################
add_library(${BIN} STATIC ${SOURCES})
target_include_directories(${BIN} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${BIN} PUBLIC ${PROJECT_INCLUDE_DIRS})
target_compile_features(${BIN} PUBLIC ${PROJECT_COMPILE_DEFINITIONS})
set_target_properties(${BIN} PROPERTIES LINKER_LANGUAGE CXX)
