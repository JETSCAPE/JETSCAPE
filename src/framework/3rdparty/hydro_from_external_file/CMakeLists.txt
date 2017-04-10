cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (lib_Hydro_from_file)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../../lib/")

set (SOURCES
     main.cpp
     FluidcellStatistic.cpp
     ParameterReader.cpp
     arsenal.cpp
     SurfaceFinder.cpp
     cornelius.cpp
     Hydroinfo_MUSIC.cpp)

add_library(hydroFromFile SHARED ${SOURCES})

