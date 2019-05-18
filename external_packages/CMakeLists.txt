cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (libJetScapeThird)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../lib/")

set (SOURCES
    fjcore.cc
    tinyxml2.cc
    cornelius.cpp
)

add_subdirectory( hydro_from_external_file )

if (${ZLIB_FOUND})
    set (SOURCES ${SOURCES}  gzstream.cc )
endif()

add_library(JetScapeThird SHARED ${SOURCES})
if (${ZLIB_FOUND})
    target_link_libraries(JetScapeThird ${ZLIB_LIBRARIES} )
endif()

