add_library(dcdflib STATIC
    dcdflib.cpp
    dcdflib.h
    cdflib.h
    ipmpar.h)

target_include_directories(dcdflib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
