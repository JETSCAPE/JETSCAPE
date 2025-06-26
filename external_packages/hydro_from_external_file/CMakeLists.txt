cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project (lib_Hydro_from_file)

# for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "../../lib/")

if(${HDF5_FOUND})
  message(STATUS "hydro from file is compiled with hdf5 library ...")
  set(SOURCES
      main.cpp
      FluidcellStatistic.cpp
      ParameterReader.cpp
      arsenal.cpp
      SurfaceFinder.cpp
      cornelius.cpp
      Hydroinfo_MUSIC.cpp
      Hydroinfo_h5.cpp)
  add_library(hydroFromFile SHARED ${SOURCES})
  set_target_properties(hydroFromFile PROPERTIES COMPILE_FLAGS "-DUSE_HDF5")
  target_link_libraries(hydroFromFile ${_hdf5_libs})
else(${HDF5_FOUND})
  message(STATUS "hydro from file is compiled without hdf5 library ...")
  set(SOURCES
      main.cpp
      FluidcellStatistic.cpp
      ParameterReader.cpp
      arsenal.cpp
      SurfaceFinder.cpp
      cornelius.cpp
      Hydroinfo_MUSIC.cpp)
  add_library(hydroFromFile SHARED ${SOURCES})
endif(${HDF5_FOUND})
